"""Feature x sample proteomics analysis. Observed and imputed matrices stay separate."""
from __future__ import annotations

import hashlib
import json
import re
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
from urllib.request import urlopen

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer
from sklearn.preprocessing import StandardScaler

plt.style.use('seaborn-v0_8-whitegrid')
matplotlib.rc('pdf', fonttype=42)
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False


FACILITY_ANNOTATIONS = ['PEPTIDE_SEQ', 'PROTEIN_UNIPROT_AC', 'FASTA_HEADERS',
                        'GENE_NAMES', 'PROTEIN_DESCRIPTIONS']
FACILITY_PROTEIN_ANNOTATIONS = ['UNIPROT_AC', 'GENE_NAMES', 'FASTA_HEADERS',
                                'PROTEIN_DESCRIPTIONS', 'NUM_PRECURSORS',
                                'PROTEIN_WEIGHT', 'NUM_UNIQUE_PEPTIDES']


def is_excel(path):
    return Path(path).suffix.lower() in {'.xlsx', '.xlsm'}


def select_excel_sheet(book, requested=None, default=None):
    if requested is not None:
        if requested not in book.sheet_names:
            raise ValueError(f'Worksheet {requested!r} absent; available sheets: {book.sheet_names}')
        return requested
    if default is not None and default in book.sheet_names:
        return default
    if len(book.sheet_names) == 1 and default is None:
        return book.sheet_names[0]
    raise ValueError(f'Specify sheet_name explicitly; available sheets: {book.sheet_names}')


def read_table(path, *, sheet_name=None, header_row=0, **kwargs):
    suffix = Path(path).suffix.lower()
    if is_excel(path):
        with pd.ExcelFile(path, engine='openpyxl') as book:
            selected = select_excel_sheet(book, sheet_name)
            return pd.read_excel(book, sheet_name=selected, header=header_row, **kwargs)
    if sheet_name is not None:
        raise ValueError('sheet_name is only valid for Excel workbooks')
    if suffix == '.parquet':
        return pd.read_parquet(path, **kwargs)
    return pd.read_csv(path, sep='\t' if suffix in {'.tsv', '.txt'} else ',', header=header_row, **kwargs)


def inspect_input(path, header_row=0):
    """Inventory schemas without selecting a measurement sheet or loading a full report."""
    def describe(columns):
        columns = [str(c) for c in columns]
        facility = set(FACILITY_ANNOTATIONS).issubset(columns)
        protein = set(FACILITY_PROTEIN_ANNOTATIONS).issubset(columns) and any(c.startswith('Raw_Intensities') for c in columns)
        diann = {'Run', 'Genes', 'Genes.MaxLFQ'}.issubset(columns)
        return {'columns': columns, 'suggested_format': 'facility_peptides' if facility else 'facility_proteins' if protein else 'diann_protein' if diann else 'matrix',
                'candidate_sample_columns': [c for c in columns if c not in FACILITY_ANNOTATIONS] if facility else [c for c in columns if c.startswith('Raw_Intensities')] if protein else [],
                'note': 'Match sample columns to metadata; numeric annotation columns are not evidence of sample identity.'}
    if is_excel(path):
        with pd.ExcelFile(path, engine='openpyxl') as book:
            return {'path': str(path), 'header_row': header_row, 'sheets': {
                name: describe(pd.read_excel(book, sheet_name=name, header=header_row, nrows=0).columns)
                for name in book.sheet_names}}
    if Path(path).suffix.lower() == '.parquet':
        import pyarrow.parquet as pq
        return {'path': str(path), **describe(pq.read_schema(path).names)}
    return {'path': str(path), **describe(read_table(path, header_row=header_row, nrows=0).columns)}


def write_json(data, path):
    def clean(x):
        if isinstance(x, dict):
            return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        if isinstance(x, np.generic):
            return clean(x.item())
        if isinstance(x, float) and not np.isfinite(x):
            return None
        return x
    Path(path).write_text(json.dumps(clean(data), indent=2, default=str, allow_nan=False) + '\n')


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def bh(pvalues):
    p = np.asarray(pvalues, float)
    q = np.full(p.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid], kind='stable')]
    if len(order):
        q[order] = np.minimum.accumulate((p[order] * len(order) / np.arange(1, len(order) + 1))[::-1])[::-1].clip(0, 1)
    return q


def validate_ids(ids, label):
    if pd.Index(ids).isna().any() or (pd.Index(ids).astype(str).str.strip() == '').any():
        raise ValueError(f'Empty {label}')
    if pd.Index(ids).duplicated().any():
        raise ValueError(f'Duplicate {label}; resolve technical replicates explicitly')


def read_inputs(config):
    inp = config['input']
    md = read_table(config['metadata_path'], sheet_name=config.get('metadata_sheet_name'), dtype=str)
    sid = config.get('sample_id_column', 'sample_id')
    validate_ids(md[sid], 'metadata sample IDs')
    md = md.set_index(sid)
    fmt = inp.get('format', 'matrix')
    idcol = inp.get('feature_id_column', 'feature_id')
    gene = inp.get('gene_column', 'gene_symbol')
    sheet = inp.get('sheet_name')
    workbook_sheets = None
    if is_excel(inp['path']):
        with pd.ExcelFile(inp['path'], engine='openpyxl') as book:
            workbook_sheets = book.sheet_names
            sheet = select_excel_sheet(book, sheet, {'facility_peptides':'Peptides (raw)', 'facility_proteins':'Summary'}.get(fmt))
    read_options = {'sheet_name': sheet, 'header_row': inp.get('header_row', 0)}
    if fmt == 'diann_protein' and Path(inp['path']).suffix.lower() == '.parquet':
        # DIA-NN reports may have millions of precursor rows and many unrelated columns.
        import pyarrow.parquet as pq
        names = pq.read_schema(inp['path']).names
        required = [inp.get('run_column', 'Run'), inp.get('gene_column', 'Genes'), inp.get('value_column', 'Genes.MaxLFQ')]
        protein_col = inp.get('protein_column', 'Protein.Ids')
        read_options['columns'] = required + ([protein_col] if protein_col in names else [])
    tab = read_table(inp['path'], **read_options)
    input_report = {'input_format': fmt, 'sheet_name': sheet, 'available_sheets': workbook_sheets,
                    'source_rows': len(tab), 'feature_level': 'phosphosite' if fmt == 'diann_phosphosite' else 'gene' if fmt == 'diann_protein' else 'peptide' if fmt == 'facility_peptides' else 'protein' if fmt == 'facility_proteins' else inp.get('feature_level', 'unspecified')}
    if fmt == 'diann_protein':
        if inp['scale'] != 'linear':
            raise ValueError('DIA-NN MaxLFQ adapter requires linear intensities')
        run, value = inp.get('run_column', 'Run'), inp.get('value_column', 'Genes.MaxLFQ')
        gene = inp.get('gene_column', 'Genes')
        protein_col = inp.get('protein_column', 'Protein.Ids')
        if tab[run].isna().any():
            raise ValueError('DIA-NN Run IDs cannot be missing')
        original_runs = tab[run].astype(str)
        tab[run] = original_runs.replace(inp.get('sample_mapping', {}))
        if tab[run].nunique() != original_runs.nunique():
            raise ValueError('sample_mapping merges distinct DIA-NN runs; resolve technical replicates explicitly')
        if tab[gene].isna().any() or tab[gene].astype(str).str.strip().eq('').any():
            raise ValueError('DIA-NN rows without gene IDs must be resolved before aggregation')
        agg = inp.get('aggregation')
        if agg not in ('mean', 'median'):
            raise ValueError('DIA-NN long input requires explicit aggregation: mean or median')
        tab[value] = pd.to_numeric(tab[value], errors='raise')
        if (tab[value].dropna() < 0).any() or np.isinf(tab[value].to_numpy()).any():
            raise ValueError('Invalid linear DIA-NN intensity')
        tab[value] = tab[value].replace(0, np.nan)
        if protein_col in tab:
            if tab[protein_col].isna().any() or tab[protein_col].astype(str).str.strip().eq('').any():
                raise ValueError('Resolve missing Protein.Ids before notebook-style aggregation')
            # Reproduce notebook weighting: collapse repeated precursor rows for each
            # gene/protein/run first, then equally aggregate protein groups per gene/run.
            grouped = tab.groupby([gene, protein_col, run], sort=False)[value].agg(agg).reset_index()
            values = grouped.groupby([gene, run], sort=False)[value].agg(agg).reset_index()
            protein_ids = tab.groupby(gene)[protein_col].agg(lambda x: ';'.join(sorted(set(x.astype(str)))))
            input_report['aggregation'] = f'{agg} within gene/protein/run, then {agg} across proteins per gene/run'
        else:
            values = tab.groupby([gene, run], sort=False)[value].agg(agg).reset_index()
            protein_ids = None
            input_report['aggregation'] = f'{agg} within gene/run; no protein identifier column supplied'
        matrix = values.pivot(index=gene, columns=run, values=value)
        ann = pd.DataFrame({'gene_symbol': matrix.index}, index=matrix.index)
        if protein_ids is not None:
            ann['protein_id'] = protein_ids.reindex(matrix.index)
        ann.index.name = 'feature_id'
    elif fmt in ('matrix', 'diann_phosphosite', 'facility_peptides', 'facility_proteins'):
        source_raw_columns = None
        if fmt == 'facility_proteins':
            missing = set(FACILITY_PROTEIN_ANNOTATIONS) - set(tab.columns)
            if missing:
                raise ValueError(f'Facility protein sheet lacks columns: {sorted(missing)}')
            prefix = inp.get('intensity_prefix', 'Raw_Intensities')
            if not prefix:
                raise ValueError('A non-empty protein intensity_prefix is required')
            source_raw_columns = [c for c in tab.columns if str(c).startswith(prefix)]
            if not source_raw_columns or any(str(c) == prefix for c in source_raw_columns):
                raise ValueError('Protein raw intensity columns must have the configured prefix followed by a sample ID')
            gene = inp.get('gene_column', 'GENE_NAMES')
            if idcol not in tab:
                tab[idcol] = tab['UNIPROT_AC']
            tab['protein_id'] = tab['UNIPROT_AC']
        if fmt == 'facility_peptides':
            required = set(FACILITY_ANNOTATIONS)
            if not required.issubset(tab.columns):
                raise ValueError(f'Facility peptide sheet lacks columns: {sorted(required - set(tab.columns))}')
            gene = inp.get('gene_column', 'GENE_NAMES')
            if idcol not in tab:
                for col in ['PEPTIDE_SEQ', 'PROTEIN_UNIPROT_AC']:
                    if tab[col].isna().any() or tab[col].astype(str).str.strip().eq('').any():
                        raise ValueError(f'Facility peptide IDs require non-empty {col}')
                tab[idcol] = tab['PEPTIDE_SEQ'].astype(str) + '|' + tab['PROTEIN_UNIPROT_AC'].astype(str)
            tab['protein_id'] = tab['PROTEIN_UNIPROT_AC']
        if fmt == 'diann_phosphosite':
            gene = inp.get('gene_column', 'Gene.Names')
            cols = [gene, inp.get('residue_column', 'Residue'), inp.get('position_column', 'Site')]
            if tab[cols].isna().any().any():
                raise ValueError('Phosphosite identifiers require gene, residue and site position')
            tab[idcol] = tab[gene].astype(str) + '_' + tab[cols[1]].astype(str) + tab[cols[2]].astype(str)
        validate_ids(tab[idcol], 'feature IDs')
        tab[idcol] = tab[idcol].astype(str)
        tab = tab.set_index(idcol)
        mapping = inp.get('sample_mapping', {})
        if source_raw_columns is not None:
            if set(mapping) - set(source_raw_columns):
                raise ValueError('Protein sample mappings may only reference raw-intensity columns')
            mapping = {c: mapping.get(c, str(c)[len(prefix):]) for c in source_raw_columns}
        input_report['effective_sample_mapping'] = mapping
        unknown = set(mapping) - set(tab.columns)
        if unknown:
            raise ValueError(f'sample_mapping source columns absent from sheet: {sorted(unknown)}')
        if fmt == 'facility_peptides' and set(mapping).intersection(FACILITY_ANNOTATIONS + ['protein_id']):
            raise ValueError('Sample mappings cannot rename peptide annotations')
        tab = tab.rename(columns=mapping)
        validate_ids(tab.columns, 'columns after sample mapping')
        sample_columns = inp.get('sample_columns', md.index.tolist())
        absent = set(md.index) - set(sample_columns)
        if absent:
            raise ValueError(f'Metadata samples omitted by sample_columns: {sorted(absent)}')
        missing_columns = set(sample_columns) - set(tab.columns)
        if missing_columns:
            raise ValueError(f'Sample columns absent from selected sheet: {sorted(missing_columns)}')
        if source_raw_columns is not None:
            raw_names = [mapping[c] for c in source_raw_columns]
            if set(sample_columns) - set(raw_names):
                raise ValueError('Protein sample_columns must select Raw_Intensities measurements, not annotation columns')
            input_report['excluded_raw_intensity_columns'] = [c for c in raw_names if c not in sample_columns]
        matrix = tab.loc[:, sample_columns].copy()
        ann = tab.drop(columns=sample_columns).copy()
        if gene in ann:
            if fmt in ('facility_peptides', 'facility_proteins'):
                ann['gene_symbol'] = ann[gene]  # Preserve the original facility annotation columns.
            else:
                ann = ann.rename(columns={gene: 'gene_symbol'})
        if 'gene_symbol' not in ann:
            raise ValueError('Supply a gene-symbol annotation column; IDs must remain separate')
        if fmt == 'facility_peptides':
            input_report['unselected_sheet_columns'] = [c for c in tab.columns if c not in sample_columns and c not in FACILITY_ANNOTATIONS + ['protein_id']]
    else:
        raise ValueError('input.format must be matrix, diann_protein, diann_phosphosite, facility_peptides, or facility_proteins')
    matrix.columns = matrix.columns.astype(str)
    validate_ids(matrix.columns, 'matrix sample IDs')
    missing = md.index.difference(matrix.columns).tolist()
    if missing:
        raise ValueError(f'Metadata samples absent from matrix: {missing}')
    excluded = matrix.columns.difference(md.index).tolist()
    matrix = matrix.loc[:, md.index].apply(pd.to_numeric, errors='raise').astype(float)
    if np.isinf(matrix.to_numpy()).any():
        raise ValueError('Infinite intensities are invalid; resolve them explicitly')
    scale = inp['scale']
    if scale == 'linear':
        if matrix.lt(0).any().any():
            raise ValueError('Negative linear intensity')
        matrix = np.log2(matrix.replace(0, np.nan))
    elif scale != 'log2':
        raise ValueError('input.scale must explicitly be linear or log2')
    empty_samples = matrix.columns[matrix.isna().all()].tolist()
    if empty_samples:
        raise ValueError(f'Entirely missing samples cannot be imputed: {empty_samples}')
    loc = inp.get('localization_column')
    report = {**input_report, 'sample_mapping': inp.get('sample_mapping', {}), 'excluded_matrix_samples': excluded, 'localization_filter': 'not applied; no column configured'}
    if loc:
        cutoff = inp.get('localization_threshold', 0.9)
        if not 0 <= cutoff <= 1:
            raise ValueError('Localization threshold must be between zero and one')
        probabilities = pd.to_numeric(ann[loc], errors='raise')
        if not probabilities.dropna().between(0, 1).all():
            raise ValueError('Localization probabilities must use the 0..1 scale')
        keep = probabilities.ge(cutoff)
        matrix, ann = matrix.loc[keep], ann.loc[keep]
        report['localization_filter'] = {'threshold': cutoff, 'excluded_features': int((~keep).sum())}
    # Preserve ambiguous protein groups, but do not expand them into invented independent genes.
    symbols = ann['gene_symbol'].fillna('').astype(str).str.strip()
    ambiguous = symbols.str.contains(r'[;,|]', regex=True) | symbols.eq('')
    ann['gene_symbol'] = symbols.mask(ambiguous)
    report['unmapped_or_ambiguous_gene_features'] = int(ambiguous.sum())
    ann.index.name = matrix.index.name = 'feature_id'
    return matrix, ann.reindex(matrix.index), md, report


def normalize(matrix, method='none'):
    if method == 'none':
        return matrix.copy(), pd.Series(0., index=matrix.columns, name='normalization_offset')
    if method != 'median':
        raise ValueError('normalization must be none or median')
    med = matrix.median()
    offsets = med - med.median()
    return matrix.subtract(offsets, axis=1), offsets.rename('normalization_offset')


def groups_for(metadata, columns):
    if not columns:
        return pd.Series('all', index=metadata.index)
    if metadata[columns].isna().any().any():
        raise ValueError('Missing imputation-group metadata')
    # Tuples avoid collisions between category labels containing separators.
    return pd.Series(list(map(tuple, metadata[columns].astype(str).to_numpy())), index=metadata.index)


def filter_observed(matrix, groups, min_total=4, min_per_group=1):
    if min_total < 1 or min_per_group < 0:
        raise ValueError('Invalid observation thresholds')
    diagnostics = pd.DataFrame({'n_observed': matrix.notna().sum(axis=1)})
    keep = diagnostics.n_observed.ge(min_total)
    for i, group in enumerate(pd.unique(groups)):
        cols = groups.index[groups.map(lambda x: x == group)]
        n = matrix[cols].notna().sum(axis=1)
        diagnostics[f'n_observed_group_{i}'] = n
        if min_per_group:
            keep &= n.ge(min_per_group)
    diagnostics['retained'] = keep
    return matrix.loc[keep].copy(), diagnostics


def impute(matrix, groups, options, seed=42):
    method = options.get('method', 'knn')
    scope = options.get('scope', 'within_group')
    if scope not in ('within_group', 'pooled'):
        raise ValueError('Imputation scope must be within_group or pooled')
    if matrix.empty or matrix.isna().all(axis=0).any():
        raise ValueError('No features or an entirely missing sample after filtering')
    filled = matrix.copy()
    rng = np.random.default_rng(seed)
    blocks = [matrix.columns] if scope == 'pooled' else [groups.index[groups.map(lambda x: x == g)] for g in pd.unique(groups)]
    for cols in blocks:
        block = matrix.loc[:, cols]
        if method == 'knn':
            orientation = options.get('neighbors', 'features')
            if orientation not in ('features', 'samples'):
                raise ValueError('knn.neighbors must be features or samples')
            if int(options.get('k', 3)) < 1:
                raise ValueError('KNN k must be positive')
            data = block if orientation == 'features' else block.T
            if data.isna().all(axis=0).any():
                raise ValueError('KNN input has an all-missing column; increase group filtering or change scope')
            imp = KNNImputer(n_neighbors=int(options.get('k', 3)), weights=options.get('weights', 'uniform'))
            arr = imp.fit_transform(data)
            filled.loc[:, cols] = arr if orientation == 'features' else arr.T
        elif method == 'gaussian':
            # Perseus-style heuristic: each sample's observed log2 mean/SD.
            shift, width = float(options.get('shift', 1.8)), float(options.get('width', 0.3))
            if shift < 0 or width <= 0:
                raise ValueError('Gaussian shift must be >=0 and width >0')
            for col in cols:
                observed = block[col].dropna()
                mask = block[col].isna()
                if not mask.any():
                    continue
                sd = observed.std(ddof=1)
                if len(observed) < 2 or not np.isfinite(sd) or sd <= 0:
                    raise ValueError(f'Cannot estimate Gaussian imputation distribution for {col}')
                filled.loc[mask, col] = rng.normal(observed.mean() - shift * sd, width * sd, int(mask.sum()))
        else:
            raise ValueError('Imputation method must be knn or gaussian')
    mask = matrix.isna()
    if filled.isna().any().any() or not np.isfinite(filled.to_numpy()).all():
        raise ValueError('Imputation left unresolved or invalid values')
    if not np.array_equal(matrix.to_numpy()[~mask], filled.to_numpy()[~mask]):
        raise AssertionError('Observed measurements changed during imputation')
    return filled, mask


def select_comparison(matrix, metadata, contrast):
    factor, test, reference = contrast['factor'], str(contrast['test']), str(contrast['reference'])
    if test == reference:
        raise ValueError('Test and reference levels must differ')
    md = metadata.copy()
    for col, values in contrast.get('subset', {}).items():
        md = md[md[col].isin([str(x) for x in values])]
    md = md[md[factor].isin([test, reference])]
    if set(md[factor]) != {test, reference}:
        raise ValueError('Both comparison levels must be present')
    excluded = []
    if contrast.get('paired', False):
        subject = contrast['subject_column']
        if md[subject].isna().any() or md.duplicated([subject, factor]).any():
            raise ValueError('Paired analysis requires exactly one sample per subject/condition; resolve duplicates')
        pairs = md.reset_index().pivot(index=subject, columns=factor, values=md.index.name)
        complete = pairs.dropna(subset=[test, reference])
        excluded = md.index[~md[subject].isin(complete.index)].tolist()
        md = md[md[subject].isin(complete.index)]
        if len(complete) < 2:
            raise ValueError('Paired analysis needs at least two metadata-complete pairs')
    if (md[factor].value_counts() < 2).any():
        raise ValueError('Each comparison group needs at least two biological samples')
    return matrix.loc[:, md.index].copy(), md, excluded


def comparison_arrays(matrix, metadata, contrast):
    factor, test, reference = contrast['factor'], str(contrast['test']), str(contrast['reference'])
    if contrast.get('paired', False):
        subject = contrast['subject_column']
        pairs = metadata.reset_index().pivot(index=subject, columns=factor, values=metadata.index.name).dropna(subset=[test, reference])
        return matrix[pairs[test]].to_numpy(), matrix[pairs[reference]].to_numpy(), pairs
    return matrix.loc[:, metadata[factor].eq(test)].to_numpy(), matrix.loc[:, metadata[factor].eq(reference)].to_numpy(), None


def ttest(matrix, metadata, contrast):
    a, b, _ = comparison_arrays(matrix, metadata, contrast)
    if np.isnan(a).any() or np.isnan(b).any():
        raise ValueError('The notebook t-test requires the imputed matrix')
    if contrast.get('paired', False):
        res = stats.ttest_rel(a, b, axis=1, nan_policy='propagate')
    else:
        res = stats.ttest_ind(a, b, axis=1, equal_var=True, nan_policy='propagate')
    out = pd.DataFrame({'log2FoldChange': a.mean(axis=1) - b.mean(axis=1), 'stat': res.statistic, 'pvalue': res.pvalue}, index=matrix.index)
    out['padj'] = bh(out.pvalue)
    out['testable'] = np.isfinite(out.pvalue)
    out['testability_reason'] = np.where(out.testable, 'tested', 'undefined_t_statistic')
    return out


def limma(matrix, metadata, contrast, options=None):
    """Notebook paired-difference model or unpaired condition/covariate model, with NAs."""
    from inmoose.limma import lmFit, squeezeVar
    from patsy import dmatrix
    opts = options or {}
    if opts.get('robust', False):
        raise ValueError('robust empirical Bayes is not implemented by the notebook InMoose backend')
    a, b, pairs = comparison_arrays(matrix, metadata, contrast)
    covariates = contrast.get('covariates', [])
    diagnostics = pd.DataFrame(index=matrix.index)
    diagnostics['n_observed_test'] = np.isfinite(a).sum(axis=1)
    diagnostics['n_observed_reference'] = np.isfinite(b).sum(axis=1)
    diagnostics['raw_missing_fraction'] = matrix.isna().mean(axis=1)
    diagnostics['mean_log2_intensity'] = matrix.mean(axis=1)
    if pairs is not None:
        if covariates:
            raise ValueError('Notebook paired-difference limma does not support covariates; request a separate design')
        y = pd.DataFrame(a - b, index=matrix.index, columns=pairs.index)
        diagnostics['n_complete_pairs'] = y.notna().sum(axis=1)
        keep = diagnostics.n_complete_pairs.ge(int(opts.get('min_complete_pairs', 4)))
        design = dmatrix('0 + paired_difference', {'paired_difference': np.ones(y.shape[1])})
        coef = 'paired_difference'
    else:
        factor, test = contrast['factor'], str(contrast['test'])
        dm = pd.DataFrame({'condition_test': metadata[factor].eq(test).astype(float)}, index=metadata.index)
        terms = ['condition_test']
        for name in covariates:
            if name == factor:
                raise ValueError('Condition cannot also be a covariate')
            if not re.fullmatch(r'[A-Za-z_]\w*', name):
                raise ValueError('Covariate names must be simple identifiers')
            if metadata[name].isna().any():
                raise ValueError(f'Missing covariate: {name}')
            if name in contrast.get('numeric_covariates', []):
                dm[name] = pd.to_numeric(metadata[name], errors='raise')
                terms.append(name)
            else:
                dm[name] = metadata[name].astype(str)
                terms.append(f'C({name})')
        design = dmatrix('1 + ' + ' + '.join(terms), dm, NA_action='raise')
        y, coef = matrix.copy(), 'condition_test'
        keep = diagnostics.n_observed_test.ge(int(opts.get('min_per_group', 2))) & diagnostics.n_observed_reference.ge(int(opts.get('min_per_group', 2)))
    x = np.asarray(design)
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError('Confounded limma design; do not drop covariates silently')
    estimable = []
    dfs = []
    for values in y.to_numpy():
        observed = np.isfinite(values)
        rank = np.linalg.matrix_rank(x[observed]) if observed.any() else 0
        dfs.append(int(observed.sum() - rank))
        estimable.append(rank == x.shape[1] and observed.sum() > rank)
    keep &= np.asarray(estimable)
    diagnostics['residual_df'] = dfs
    diagnostics['testable'] = keep
    diagnostics['testability_reason'] = np.where(keep, 'tested', 'insufficient_observations_pairs_or_residual_df')
    if keep.sum() < 2:
        raise ValueError('Fewer than two features are testable by limma; revise the design/threshold explicitly')
    fit = lmFit(y.loc[keep], design=design)
    # For differences Amean would be mean difference; use original abundance for the trend.
    trend = diagnostics.loc[keep, 'mean_log2_intensity'].to_numpy() if opts.get('trend', True) else False
    # Use InMoose's empirical-Bayes moderation directly. Its eBayes wrapper's
    # unrelated B-statistic branch raises KeyError when df_prior is scalar infinity
    # in 0.9.1. The moderated t/df/p formulas below are exactly its _ebayes formulas.
    moderation = squeezeVar(fit.sigma ** 2, fit.df_residual,
                            covariate=trend if not isinstance(trend, bool) else None,
                            robust=False)
    moderated_t = fit.coefficients[coef].to_numpy() / (
        fit.stdev_unscaled[coef].to_numpy() * np.sqrt(moderation['var_post']))
    total_df = np.minimum(fit.df_residual + moderation['df_prior'], np.nansum(fit.df_residual))
    pvalues = 2 * stats.t.cdf(-np.abs(moderated_t), df=total_df)
    res = pd.DataFrame({'log2FoldChange': fit.coefficients[coef], 'stat': moderated_t,
                        'pvalue': pvalues, 'padj': bh(pvalues),
                        'moderated_df': total_df,
                        'posterior_variance': moderation['var_post']}, index=y.index[keep])
    out = diagnostics.join(res)
    undefined = out.testable & ~np.isfinite(out.pvalue)
    out.loc[undefined, 'testable'] = False
    out.loc[undefined, 'testability_reason'] = 'undefined_moderated_statistic'
    return out, pd.DataFrame(x, index=y.columns, columns=design.design_info.column_names)


def savefig(fig, stem, formats=('pdf', 'png')):
    for fmt in formats:
        fig.savefig(f'{stem}.{fmt}', bbox_inches='tight', dpi=180)
    plt.close(fig)


def placeholder(stem, text, formats):
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.axis('off')
    ax.text(0.5, 0.5, text, ha='center', va='center', wrap=True)
    savefig(fig, stem, formats)


def plot_pca(matrix, metadata, color, stem, formats=('pdf', 'png'), top_features=500, labels=False, palette=None):
    variance = matrix.std(axis=1)
    features = variance[variance.gt(0)].sort_values(ascending=False, kind='stable').head(top_features).index
    if min(len(features), matrix.shape[1] - 1) < 2:
        raise ValueError('PCA needs >=2 variable features and >=3 samples')
    scaled = StandardScaler().fit_transform(matrix.loc[features].T)
    pca = PCA(n_components=2, svd_solver='full')
    scores = pd.DataFrame(pca.fit_transform(scaled), index=matrix.columns, columns=['PC1', 'PC2'])
    scores.join(metadata).to_csv(f'{stem}_scores.csv')
    pd.DataFrame(pca.components_.T, index=features, columns=['PC1', 'PC2']).to_csv(f'{stem}_loadings.csv')
    pd.DataFrame({'PC': ['PC1', 'PC2'], 'variance_explained': pca.explained_variance_ratio_}).to_csv(f'{stem}_variance.csv', index=False)
    fig, ax = plt.subplots(figsize=(7, 6))
    for level in metadata[color].unique():
        idx = metadata.index[metadata[color].eq(level)]
        ax.scatter(scores.loc[idx, 'PC1'], scores.loc[idx, 'PC2'], s=50, alpha=.8, label=level, color=(palette or {}).get(level))
    if labels:
        for sample, row in scores.iterrows():
            ax.annotate(sample, (row.PC1, row.PC2), fontsize=7)
    ax.set(xlabel=f'PC1 ({100*pca.explained_variance_ratio_[0]:.2f}%)', ylabel=f'PC2 ({100*pca.explained_variance_ratio_[1]:.2f}%)')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.grid(alpha=.2)
    savefig(fig, stem, formats)
    return features.tolist()


def plot_qc(raw, filled, mask, groups, out, formats=('pdf', 'png')):
    per_sample = pd.DataFrame({'n_observed': raw.notna().sum(), 'n_imputed': mask.sum(), 'imputed_fraction': mask.mean(), 'group': groups.map(str)})
    per_feature = pd.DataFrame({'n_observed': raw.notna().sum(axis=1), 'n_imputed': mask.sum(axis=1), 'imputed_fraction': mask.mean(axis=1)})
    per_sample.to_csv(out / 'imputation_per_sample.csv')
    per_feature.to_csv(out / 'imputation_per_feature.csv')
    per_sample.groupby('group')[['n_observed', 'n_imputed']].sum().assign(imputed_fraction=lambda d: d.n_imputed / (d.n_observed + d.n_imputed)).to_csv(out / 'imputation_per_group.csv')
    fig, ax = plt.subplots(figsize=(max(7, .35*len(per_sample)), 4))
    ax.bar(per_sample.index, per_sample.n_observed, label='Observed')
    ax.bar(per_sample.index, per_sample.n_imputed, bottom=per_sample.n_observed, label='Imputed', color='#D55E00')
    ax.tick_params(axis='x', labelrotation=90)
    ax.set_ylabel('Features')
    ax.legend()
    savefig(fig, out / 'observed_imputed_counts', formats)
    # Paginate so hundreds of samples never produce a single unreadable figure.
    for page, start in enumerate(range(0, raw.shape[1], 12), 1):
        cols = raw.columns[start:start+12]
        fig, axes = plt.subplots(int(np.ceil(len(cols)/3)), 3, figsize=(12, 3*np.ceil(len(cols)/3)), squeeze=False)
        for ax, col in zip(axes.flat, cols):
            ax.hist([raw[col].dropna(), filled.loc[mask[col], col]], bins=20, stacked=True, label=['Observed', 'Imputed'], color=['#0072B2', '#D55E00'])
            ax.set(title=str(col), xlabel='Log2 intensity', ylabel='Count')
        for ax in list(axes.flat)[len(cols):]:
            ax.axis('off')
        axes.flat[0].legend(fontsize=8)
        fig.tight_layout()
        savefig(fig, out / f'observed_imputed_distributions_{page}', formats)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, data, title in zip(axes, [raw, filled], ['Observed (pairwise complete)', 'After imputation']):
        corr = data.corr(min_periods=3)
        corr.to_csv(out / ('correlations_observed.csv' if title.startswith('Observed') else 'correlations_imputed.csv'))
        img = ax.imshow(corr, vmin=-1, vmax=1, cmap='coolwarm')
        ax.set(title=title, xticks=range(len(corr)), yticks=range(len(corr)), xticklabels=corr.columns, yticklabels=corr.index)
        ax.tick_params(axis='x', labelrotation=90, labelsize=6)
        ax.tick_params(axis='y', labelsize=6)
    fig.colorbar(img, ax=axes, shrink=.6)
    savefig(fig, out / 'correlation_comparison', formats)
    return per_feature


def volcano(res, ann, stem, title, alpha, formats):
    fig, ax = plt.subplots(figsize=(7, 6))
    y = -np.log10(res.pvalue.clip(lower=1e-300))
    sig = res.padj.lt(alpha)
    ax.scatter(res.log2FoldChange, y, c=np.where(sig, '#D55E00', '#aaaaaa'), s=16, alpha=.7)
    for feature in res.loc[sig].sort_values('padj').head(15).index:
        label = ann.at[feature, 'gene_symbol']
        ax.annotate(str(label) if pd.notna(label) else str(feature), (res.at[feature, 'log2FoldChange'], y.at[feature]), fontsize=7)
    ax.set(title=title, xlabel='Log2 fold change (test - reference)', ylabel='-log10 nominal p-value')
    savefig(fig, stem, formats)


def imputation_audit(tres, lres, per_feature, out, alpha, formats):
    audit = per_feature.join(tres[['log2FoldChange', 'padj']].add_prefix('ttest_')).join(lres[['log2FoldChange', 'padj', 'testable', 'testability_reason']].add_prefix('limma_'))
    hit = audit.ttest_padj.lt(alpha)
    testable = audit.limma_testable.fillna(False)
    both = hit & testable
    retained = both & audit.limma_padj.lt(alpha)
    same = np.sign(audit.ttest_log2FoldChange) == np.sign(audit.limma_log2FoldChange)
    summary = {
        'matrix_entries': int(per_feature.n_observed.sum() + per_feature.n_imputed.sum()),
        'imputed_entries': int(per_feature.n_imputed.sum()),
        'imputed_fraction': float(per_feature.n_imputed.sum() / (per_feature.n_observed.sum() + per_feature.n_imputed.sum())),
        'ttest_significant_hits': int(hit.sum()),
        'hits_with_no_imputation': int((hit & audit.n_imputed.eq(0)).sum()),
        'hits_with_imputation': int((hit & audit.n_imputed.gt(0)).sum()),
        'fraction_hits_with_no_imputation': float((hit & audit.n_imputed.eq(0)).sum() / hit.sum()) if hit.any() else None,
        'median_imputed_values_among_imputed_hits': float(audit.loc[hit & audit.n_imputed.gt(0), 'n_imputed'].median()),
        'hits_testable_by_limma': int(both.sum()),
        'hits_untestable_by_limma': int((hit & ~testable).sum()),
        'hits_significant_by_limma': int(retained.sum()),
        'fraction_testable_hits_significant_by_limma': float(retained.sum() / both.sum()) if both.any() else None,
        'fraction_testable_hits_same_direction': float((both & same).sum() / both.sum()) if both.any() else None,
        'interpretation': 'Differences combine imputation, feature eligibility, and t-test versus empirical-Bayes modeling; not a causal estimate of imputation alone.'
    }
    audit['ttest_significant'] = hit
    audit['limma_significant'] = audit.limma_padj.lt(alpha)
    audit.to_csv(out / 'imputation_method_audit.csv')
    audit.loc[hit].to_csv(out / 'significant_hits_imputation_audit.csv')
    write_json(summary, out / 'imputation_impact_summary.json')
    fig, ax = plt.subplots(figsize=(6, 6))
    ok = audit.ttest_log2FoldChange.notna() & audit.limma_log2FoldChange.notna()
    points = ax.scatter(audit.loc[ok, 'ttest_log2FoldChange'], audit.loc[ok, 'limma_log2FoldChange'], c=audit.loc[ok, 'imputed_fraction'], cmap='viridis', s=15)
    fig.colorbar(points, ax=ax, label='Imputed fraction')
    bounds = np.r_[audit.loc[ok, 'ttest_log2FoldChange'], audit.loc[ok, 'limma_log2FoldChange']]
    if len(bounds):
        ax.plot([bounds.min(), bounds.max()], [bounds.min(), bounds.max()], '--', c='gray')
    ax.set(xlabel='Imputed t-test log2FC', ylabel='Non-imputed limma log2FC')
    savefig(fig, out / 'effect_size_agreement', formats)
    fig, ax = plt.subplots(figsize=(7, 4))
    categories = {'Both significant': int(retained.sum()), 'T-test only (limma testable)': int((both & ~audit.limma_padj.lt(alpha)).sum()), 'T-test hit (limma untestable)': int((hit & ~testable).sum()), 'Limma only': int((~hit & audit.limma_padj.lt(alpha)).sum())}
    ax.barh(list(categories), list(categories.values()), color='#0072B2')
    ax.set_xlabel('Features')
    savefig(fig, out / 'significant_hit_agreement', formats)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist([audit.loc[hit, 'imputed_fraction'], audit.loc[~hit, 'imputed_fraction']], bins=np.linspace(0, 1, 11), label=['T-test hits', 'Other features'])
    ax.set(xlabel='Fraction imputed', ylabel='Features')
    ax.legend()
    savefig(fig, out / 'imputation_burden_among_hits', formats)
    return summary


def read_gmt(path):
    sets = {}
    for line in Path(path).read_text().splitlines():
        fields = line.rstrip().split('\t')
        if len(fields) >= 3:
            if fields[0] in sets:
                raise ValueError(f'Duplicate GMT term: {fields[0]}')
            sets[fields[0]] = set(fields[2:]) - {''}
    if not sets:
        raise ValueError('No gene sets in GMT')
    return sets


def get_gene_sets(config, output):
    opts = config.get('enrichment', {})
    if not opts.get('enabled', True):
        return None, {'status': 'disabled explicitly'}
    path = opts.get('gmt_path')
    release = opts.get('release', '2024.1.Hs')
    if not path:
        if config.get('organism', 'human') != 'human':
            raise ValueError('Non-human analysis requires an explicitly matched GMT; never use human sets silently')
        if not re.fullmatch(r'\d{4}\.\d+\.Hs', release):
            raise ValueError('Invalid MSigDB human release name')
        path = output / f'h.all.v{release}.symbols.gmt'
        url = f'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{release}/h.all.v{release}.symbols.gmt'
        try:
            with urlopen(url, timeout=45) as response:
                Path(path).write_bytes(response.read())
        except Exception as exc:
            raise RuntimeError('Hallmark retrieval failed; supply enrichment.gmt_path. Do not report a complete enrichment run.') from exc
    return read_gmt(path), {'gmt_path': str(path), 'sha256': file_sha256(path), 'release': release if not opts.get('gmt_path') else opts.get('release', 'user supplied')}


def enrichment(res, ann, gene_sets, out, alpha, formats):
    # Unique mapped genes among testable features; repeated phosphosites count once.
    named = res.join(ann[['gene_symbol']])
    background = set(named.loc[named.pvalue.notna(), 'gene_symbol'].dropna())
    membership = {feature: [] for feature in named.index}
    for direction, sign in [('up', 1), ('down', -1)]:
        significant = named.padj.lt(alpha) & (named.log2FoldChange * sign).gt(0)
        hits = set(named.loc[significant, 'gene_symbol'].dropna()) & background
        rows = []
        for term, genes in gene_sets.items():
            eligible = genes & background
            overlap = eligible & hits
            p = stats.hypergeom.sf(len(overlap)-1, len(background), len(eligible), len(hits)) if background else 1.
            rows.append({'Term': term, 'background_genes': len(background), 'set_genes_in_background': len(eligible), 'directional_hit_genes': len(hits), 'overlap': len(overlap), 'pvalue': float(p), 'Genes': ';'.join(sorted(overlap))})
        table = pd.DataFrame(rows)
        table['padj'] = bh(table.pvalue)
        table['fdr_threshold'], table['significant'] = alpha, table.padj.lt(alpha)
        table = table.sort_values(['padj', 'pvalue', 'Term'])
        table.to_csv(out / f'hallmark_{direction}_all.csv', index=False)
        top = table.head(10)
        top.to_csv(out / f'hallmark_{direction}_top10.csv', index=False)
        fig, ax = plt.subplots(figsize=(9, max(4, .4*len(top))))
        ax.barh(top.Term[::-1], -np.log10(top.padj[::-1].clip(lower=1e-300)), color=np.where(top.significant[::-1], '#D55E00', '#aaaaaa'))
        ax.axvline(-np.log10(alpha), linestyle='--', color='black', label=f'FDR {alpha:g}')
        ax.set(title=f'Hallmark {direction}; {len(hits)} unique hit genes', xlabel='-log10 adjusted p-value')
        ax.legend()
        savefig(fig, out / f'hallmark_{direction}_top10', formats)
        contributors = []
        for row in top.itertuples():
            genes = set(row.Genes.split(';')) - {''}
            for feature, values in named.loc[significant & named.gene_symbol.isin(genes)].iterrows():
                contributors.append({'Term': row.Term, 'feature_id': feature, 'gene_symbol': values.gene_symbol, 'log2FoldChange': values.log2FoldChange})
        contributors = pd.DataFrame(contributors, columns=['Term', 'feature_id', 'gene_symbol', 'log2FoldChange'])
        contributors.to_csv(out / f'hallmark_{direction}_contributors.csv', index=False)
        if contributors.empty:
            placeholder(out / f'hallmark_{direction}_contributors', 'No significant directional contributors', formats)
        else:
            fig, ax = plt.subplots(figsize=(9, max(4, .4*contributors.Term.nunique())))
            terms = list(dict.fromkeys(contributors.Term))
            ax.scatter(contributors.log2FoldChange, contributors.Term.map({v: i for i, v in enumerate(terms)}), s=20, alpha=.6)
            ax.set(yticks=range(len(terms)), yticklabels=terms, xlabel='Feature log2 fold change')
            savefig(fig, out / f'hallmark_{direction}_contributors', formats)
        for row in table.loc[table.significant].itertuples():
            genes = set(row.Genes.split(';')) - {''}
            for feature in named.index[significant & named.gene_symbol.isin(genes)]:
                membership[feature].append(f'{direction}:{row.Term}')
    return pd.Series({k: ';'.join(v) for k, v in membership.items()}, name='Hallmark')


def adjust_phosphosites(raw, ann, config):
    opts = config['protein_adjustment']
    protein = read_table(opts['path'], sheet_name=opts.get('sheet_name'))
    idcol = opts.get('feature_id_column', 'feature_id')
    validate_ids(protein[idcol], 'matched protein IDs')
    protein[idcol] = protein[idcol].astype(str)
    protein = protein.set_index(idcol).loc[:, raw.columns].apply(pd.to_numeric, errors='raise').astype(float)
    if np.isinf(protein.to_numpy()).any():
        raise ValueError('Infinite matched protein values')
    if opts['scale'] == 'linear':
        if protein.lt(0).any().any():
            raise ValueError('Negative linear protein intensities')
        protein = np.log2(protein.replace(0, np.nan))
    elif opts['scale'] != 'log2':
        raise ValueError('Matched protein scale must be linear or log2')
    mapping = ann[opts.get('mapping_column', 'protein_id')]
    mapped = mapping.notna() & mapping.astype(str).isin(protein.index)
    matched = protein.reindex(mapping.loc[mapped].astype(str))
    matched.index = raw.index[mapped]
    adjusted = raw.loc[mapped] - matched
    report = {'mapped_sites': int(mapped.sum()), 'unmapped_sites': int((~mapped).sum()), 'additional_missing_entries': int((adjusted.isna() & raw.loc[mapped].notna()).sum().sum()), 'interpretation': 'Log2 phosphosite minus matched log2 protein abundance; relative adjustment, not absolute occupancy.'}
    return adjusted, ann.loc[mapped].copy(), report


def dependency_versions():
    versions = {}
    for package in ['numpy', 'pandas', 'scipy', 'scikit-learn', 'matplotlib', 'inmoose', 'patsy', 'pyarrow', 'openpyxl']:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = 'unavailable'
    return versions
