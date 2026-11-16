"""Run a configured proteomics/phosphoproteomics analysis, without shared imputation state."""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
import proteomics as p


def scenario_name(options):
    name = options.get('name', options.get('method', 'knn') + '_' + options.get('scope', 'within_group'))
    if not re.fullmatch(r'[A-Za-z0-9_-]+', name):
        raise ValueError('Imputation scenario names may use letters, digits, _ and - only')
    return name


def run(config):
    out = Path(config['output_dir'])
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError('output_dir must be empty: use a fresh directory to prevent stale results')
    p.write_json(config, out / 'analysis_config.json')
    summary = {'status': 'running', 'versions': p.dependency_versions(), 'source': 'Trulsson/biopsy_vs_resection', 'contrasts': {}, 'scopes': {}}
    p.write_json(summary, out / 'run_summary.json')
    try:
        # Fail before producing apparently complete results when the statistical backend is missing.
        from inmoose.limma import lmFit
        raw, ann, md, report = p.read_inputs(config)
        if not config.get('contrasts'):
            raise ValueError('At least one explicit contrast is required')
        alpha = config.get('alpha', .05)
        if not 0 < alpha < 1:
            raise ValueError('alpha must lie between zero and one')
        formats = config.get('formats', ['pdf', 'png'])
        if not formats or set(formats) - {'pdf', 'png', 'svg'}:
            raise ValueError('formats must contain pdf, png and/or svg')
        raw.to_csv(out / 'log2_observed_input.csv')
        ann.to_csv(out / 'feature_annotations.csv')
        md.to_csv(out / 'sample_metadata.csv')
        summary['input_report'] = report
        summary['input_sha256'] = {key: p.file_sha256(path) for key, path in {'matrix': config['input']['path'], 'metadata': config['metadata_path']}.items()}
        summary['pairing'] = {c['name']: c.get('paired', False) for c in config['contrasts']}
        scenario_options = config.get('imputations', [{'method': 'knn', 'scope': 'within_group', 'neighbors': 'features', 'k': 3, 'weights': 'uniform'}])
        if not scenario_options:
            raise ValueError('At least one imputation scenario is required for PCA and notebook t-tests')
        names = [scenario_name(o) for o in scenario_options]
        if len(set(names)) != len(names):
            raise ValueError('Duplicate imputation scenario names')
        contrast_names = [c['name'] for c in config['contrasts']]
        if len(set(contrast_names)) != len(contrast_names):
            raise ValueError('Contrast names must be unique')
        gene_sets, summary['enrichment'] = p.get_gene_sets(config, out)
        assays = {'abundance': (raw, ann)}
        if config.get('protein_adjustment'):
            if config.get('assay', 'proteomics') != 'phosphoproteomics':
                raise ValueError('protein_adjustment requires assay=phosphoproteomics')
            adjusted, adjusted_ann, adjust_report = p.adjust_phosphosites(raw, ann, config)
            assays['protein_adjusted'] = (adjusted, adjusted_ann)
            summary['protein_adjustment'] = adjust_report
            summary['input_sha256']['matched_protein'] = p.file_sha256(config['protein_adjustment']['path'])
            adjusted.to_csv(out / 'protein_adjusted_log2_observed.csv')
        scope_filter = config.get('filtering', {})
        for assay, (assay_raw, assay_ann) in assays.items():
            # Global scope starts from the original matrix, never a union of contrast-filtered matrices.
            global_groups = config.get('global_group_columns', [config.get('condition_column', 'condition')])
            scopes = [('all_samples', assay_raw, md, global_groups, None)]
            for contrast in config['contrasts']:
                if not re.fullmatch(r'[A-Za-z0-9_-]+', contrast['name']):
                    raise ValueError('Contrast names may use letters, digits, _ and - only')
                if 'paired' in contrast and not isinstance(contrast['paired'], bool):
                    raise ValueError('paired must be a JSON boolean')
                subset, subset_md, excluded = p.select_comparison(assay_raw, md, contrast)
                summary['contrasts'][contrast['name']] = {'paired': contrast.get('paired', False), 'excluded_incomplete_pair_samples': excluded, 'test': contrast['test'], 'reference': contrast['reference'], 'subset': contrast.get('subset', {})}
                scopes.append((contrast['name'], subset, subset_md, contrast.get('imputation_group_columns', [contrast['factor']]), contrast))
            for scope, original, scope_md, grouping, contrast in scopes:
                dest = out / assay / scope
                dest.mkdir(parents=True)
                groups = p.groups_for(scope_md, grouping)
                normalized, offsets = p.normalize(original, config.get('normalization', 'none'))
                offsets.to_csv(dest / 'normalization_offsets.csv')
                group_map = pd.DataFrame({'group_id': groups.map(str)}, index=scope_md.index)
                group_map.to_csv(dest / 'imputation_groups.csv')
                opts = dict(scope_filter)
                if contrast:
                    opts.update(contrast.get('filtering', {}))
                retained, filters = p.filter_observed(normalized, groups, opts.get('min_total', 4), opts.get('min_per_group', 1))
                filters.join(assay_ann, rsuffix='_annotation').to_csv(dest / 'feature_filtering.csv')
                retained.to_csv(dest / 'observed_log2_filtered.csv')
                if retained.empty:
                    raise ValueError(f'No features retained for {assay}/{scope}')
                scope_key = f'{assay}/{scope}'
                summary['scopes'][scope_key] = {'sample_ids': scope_md.index.tolist(), 'features_input': len(original), 'features_retained': len(retained), 'imputation_group_columns': grouping, 'scenarios': {}}
                lres = None
                if contrast:
                    lres, design = p.limma(retained, scope_md, contrast, config.get('limma', {}))
                    design.to_csv(dest / 'limma_design.csv')
                    annotated = assay_ann.reindex(lres.index).join(lres)
                    if gene_sets is not None:
                        membership = p.enrichment(lres, assay_ann, gene_sets, dest, alpha, formats)
                        annotated = annotated.join(membership)
                    annotated.to_csv(dest / 'limma_nonimputed_results.csv')
                    p.volcano(lres, assay_ann, dest / 'volcano_limma_nonimputed', f'{scope}: non-imputed limma', alpha, formats)
                    summary['scopes'][scope_key]['limma_testable_features'] = int(lres.testable.sum())
                for method_options in scenario_options:
                    method_name = scenario_name(method_options)
                    method_dest = dest / method_name
                    method_dest.mkdir()
                    filled, mask = p.impute(retained, groups, method_options, config.get('seed', 42))
                    filled.to_csv(method_dest / 'imputed_log2.csv')
                    mask.astype(int).to_csv(method_dest / 'imputation_mask.csv')
                    per_feature = p.plot_qc(retained, filled, mask, groups, method_dest, formats)
                    color = contrast['factor'] if contrast else config.get('pca_color', config.get('condition_column', 'condition'))
                    if scope_md[color].isna().any():
                        raise ValueError('Missing PCA color metadata')
                    used = p.plot_pca(filled, scope_md, color, method_dest / 'pca', formats, config.get('pca_top_features', 500), config.get('pca_labels', False), config.get('palette'))
                    method_summary = {'pca_features': used, 'imputed_entries': int(mask.sum().sum()), 'imputed_fraction': float(mask.to_numpy().mean()), 'options': method_options}
                    if contrast:
                        tres = p.ttest(filled, scope_md, contrast)
                        annotated = assay_ann.reindex(tres.index).join(tres).join(per_feature)
                        if gene_sets is not None:
                            membership = p.enrichment(tres, assay_ann, gene_sets, method_dest, alpha, formats)
                            annotated = annotated.join(membership)
                        annotated.to_csv(method_dest / 'ttest_imputed_results.csv')
                        p.volcano(tres, assay_ann, method_dest / 'volcano_ttest_imputed', f'{scope}: imputed t-test', alpha, formats)
                        method_summary['impact'] = p.imputation_audit(tres, lres, per_feature, method_dest, alpha, formats)
                    summary['scopes'][scope_key]['scenarios'][method_name] = method_summary
        summary['status'] = 'complete'
        p.write_json(summary, out / 'run_summary.json')
        methods = (
            'Proteomics/phosphoproteomics analysis derived from Trulsson/biopsy_vs_resection.\n'
            'Source blobs: proteomics 3f82fbdabc039c0e3926010005a582dde27743d2; phosphoproteomics 667c1425fdb714e8337d10998da3ee565775290f.\n'
            f"Input format and feature level: {report['input_format']}, {report['feature_level']}; sheet {report['sheet_name']}. No automatic peptide-to-protein summarization.\n"
            'Features are rows; samples are columns. Linear zero intensities become missing before log2; log2 zeros remain valid.\n'
            'Each all-sample and comparison scope starts from the original observed matrix and independently normalizes, filters, imputes, selects variable features, standardizes and fits PCA.\n'
            'PCA uses the most variable retained features after imputation, with per-feature centering/unit variance.\n'
            'Default imputation: within-condition feature-neighbor KNN, k=3, uniform weights, nan-euclidean distance. Gaussian, when requested, uses each sample observed mean minus shift*SD and width*SD in log2 space.\n'
            'Imputed testing: Student equal-variance two-sided t-test for unpaired data; paired two-sided t-test only for explicitly paired contrasts.\n'
            'Non-imputed InMoose limma: paired differences with >=4 complete pairs by default; unpaired group/covariate model with >=2 observations per group and positive residual df. Abundance-dependent empirical Bayes, robust=False.\n'
            'Limma uses InMoose lmFit and squeezeVar with the identical moderated-t/df/p formulas as its eBayes; B-statistics are omitted to avoid the upstream scalar-infinite-prior edge case.\n'
            'BH adjustment is performed separately for each test family over finite p-values. Untestable features remain explicit, with missing statistics.\n'
            'Imputation audit compares two different statistical tests; agreement does not isolate a causal imputation effect or establish imputed-value accuracy. Limma does not correct non-random missingness.\n'
            'Enrichment: separate up/down hypergeometric ORA with BH across all supplied pathways, unique unambiguous mapped genes in each method testable background. Top ten are exported and plotted even when not significant.\n'
            'Phosphosite Hallmark results describe the genes carrying changing sites, not kinase activity; repeated sites are counted once per gene and direction. A gene with sites changing in both directions can occur in both directional analyses.\n'
            'Protein adjustment, when requested, subtracts matched observed protein log2 abundance per sample before independent processing, and is not absolute occupancy.\n'
            'Volcano y-axis: nominal p-value; significance/coloring: BH-adjusted p-value. PDF fonts are editable.\n'
            'Exact parameters, subset/pair exclusions, input/GMT SHA256 hashes, and dependency versions are in analysis_config.json and run_summary.json.\n'
            'Configuration:\n' + json.dumps(config, indent=2) + '\nVersions:\n' + json.dumps(summary['versions'], indent=2) + '\n'
        )
        (out / 'methods.txt').write_text(methods)
        return summary
    except Exception as exc:
        summary['status'], summary['error'] = 'failed', f'{type(exc).__name__}: {exc}'
        p.write_json(summary, out / 'run_summary.json')
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--config')
    selection.add_argument('--inspect', metavar='INPUT_PATH', help='List worksheets/columns before configuring input')
    parser.add_argument('--header-row', type=int, default=0, help='Zero-based header row for --inspect')
    args = parser.parse_args()
    if args.inspect:
        print(json.dumps(p.inspect_input(args.inspect, args.header_row), indent=2))
        return
    config_path = Path(args.config).resolve()
    config = json.loads(config_path.read_text())
    # Relative paths are relative to the configuration, independent of invocation cwd.
    def resolve(path):
        path = Path(path).expanduser()
        return str(path if path.is_absolute() else config_path.parent / path)
    config['input']['path'] = resolve(config['input']['path'])
    config['metadata_path'] = resolve(config['metadata_path'])
    config['output_dir'] = resolve(config['output_dir'])
    if config.get('enrichment', {}).get('gmt_path'):
        config['enrichment']['gmt_path'] = resolve(config['enrichment']['gmt_path'])
    if config.get('protein_adjustment'):
        config['protein_adjustment']['path'] = resolve(config['protein_adjustment']['path'])
    summary = run(config)
    print(json.dumps({'status': summary['status'], 'output_dir': config['output_dir'], 'scopes': len(summary['scopes'])}))


if __name__ == '__main__':
    main()
