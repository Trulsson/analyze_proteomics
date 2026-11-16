"""Synthetic regression check: run directly; creates temporary analysis outputs only."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import numpy as np
import pandas as pd
from scipy import stats
import proteomics as p
from run_proteomics import run


def check():
    rng = np.random.default_rng(71)
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        # Interleave and reverse condition order to detect accidental positional pairing.
        md = pd.DataFrame({'sample_id': [f's{i}' for i in range(16)], 'condition': ['B', 'A']*8, 'subject': np.repeat([f'p{i}' for i in range(8)], 2), 'batch': ['x', 'x', 'y', 'y']*4}).set_index('sample_id')
        data = rng.normal(20, 1, (90, 16))
        data[:10, md.condition.eq('B')] += 3
        data[10:20, md.condition.eq('B')] -= 3
        raw = pd.DataFrame(data, index=[f'f{i}' for i in range(90)], columns=md.index)
        raw = raw.mask(rng.random(raw.shape) < .08)
        # Present in enough samples overall but too few complete pairs for paired limma.
        raw.loc['f89'] = np.nan
        raw.loc['f89', ['s0', 's2', 's4', 's7', 's9', 's11']] = [20., 21., 19., 20., 22., 18.]
        groups = p.groups_for(md, ['condition'])
        retained, filtering = p.filter_observed(raw, groups)
        filled, mask = p.impute(retained, groups, {})
        assert np.array_equal(filled.values[~mask], retained.values[~mask])
        assert not filled.isna().any().any()
        gaussian, _ = p.impute(retained, groups, {'method': 'gaussian'})
        gaussian2, _ = p.impute(retained, groups, {'method': 'gaussian'})
        assert np.array_equal(gaussian, gaussian2)
        assert np.array_equal(gaussian.values[~mask], retained.values[~mask])
        np.testing.assert_allclose(p.bh([.01, .03, np.nan, .2]), [.03, .045, np.nan, .2], equal_nan=True)
        unpaired = {'name': 'unpaired', 'factor': 'condition', 'test': 'B', 'reference': 'A', 'paired': False}
        paired = dict(unpaired, name='paired', paired=True, subject_column='subject')
        tres = p.ttest(filled, md, unpaired)
        np.testing.assert_allclose(tres.pvalue, stats.ttest_ind(filled.loc[:, md.condition.eq('B')], filled.loc[:, md.condition.eq('A')], axis=1).pvalue)
        # The result must be invariant to metadata order when pairing is by subject ID.
        pairres = p.ttest(filled, md, paired)
        shuffled = md.sample(frac=1, random_state=2)
        np.testing.assert_allclose(pairres.pvalue, p.ttest(filled, shuffled, paired).pvalue)
        for contrast in [unpaired, paired]:
            res, design = p.limma(retained, md, contrast)
            assert res.loc['f0', 'log2FoldChange'] > 0
            assert res.loc['f10', 'log2FoldChange'] < 0
            assert res.index.equals(retained.index)
            if contrast['paired']:
                assert not res.loc['f89', 'testable'] and pd.isna(res.loc['f89', 'padj'])
        confounded = md.assign(confounded=md.condition)
        try:
            p.limma(retained, confounded, dict(unpaired, covariates=['confounded']))
        except ValueError as exc:
            assert 'Confounded' in str(exc)
        else:
            raise AssertionError('Confounding was not rejected')
        # Match the standard eBayes/topTable path when its prior df is finite.
        from inmoose.limma import lmFit, eBayes, topTable
        from patsy import dmatrix
        parity_rng = np.random.default_rng(328)
        scale = np.exp(parity_rng.normal(0, 1, 160))
        parity = pd.DataFrame(parity_rng.normal(size=(160, 16))*scale[:, None]+20,
                              index=[f'F{i}' for i in range(160)], columns=md.index)
        parity_md = md.assign(condition=['A']*8+['B']*8)
        check_res, _ = p.limma(parity, parity_md, unpaired)
        dm = dmatrix('1 + condition_test', {'condition_test': np.r_[np.zeros(8), np.ones(8)]})
        standard = eBayes(lmFit(parity, design=dm), trend=parity.mean(axis=1).to_numpy(), robust=False)
        standard = pd.DataFrame(topTable(standard, coef='condition_test', number=160, sort_by='P', adjust_method='fdr_bh')).reindex(parity.index)
        np.testing.assert_allclose(check_res.pvalue, standard.pvalue, rtol=1e-12)
        np.testing.assert_allclose(check_res.padj, standard.adj_pvalue, rtol=1e-12)
        annotations = pd.DataFrame({'gene_symbol': [f'G{i}' for i in range(90)], 'protein_id': [f'P{i}' for i in range(90)]}, index=raw.index)
        # A feature visible only in the first comparison should survive that comparison's
        # own filtering even if filtering on the global grouping would remove it.
        comparison_groups = p.groups_for(md.iloc[:8], ['condition'])
        scoped, _ = p.filter_observed(raw.iloc[:, :8], comparison_groups)
        assert set(scoped.index) <= set(raw.index)
        annotations.join(raw).reset_index(names='feature_id').to_csv(root / 'matrix.csv', index=False)
        md.to_csv(root / 'metadata.csv')
        (root / 'sets.gmt').write_text('HALLMARK_UP\tna\t'+'\t'.join(f'G{i}' for i in range(10))+'\nHALLMARK_DOWN\tna\t'+'\t'.join(f'G{i}' for i in range(10,20))+'\nHALLMARK_OTHER\tna\tG50\tG51\n')
        ingestion_cfg = {'input': {'path': str(root/'matrix.csv'), 'scale':'log2'}, 'metadata_path':str(root/'metadata.csv')}
        read_raw, read_ann, read_md, report = p.read_inputs(ingestion_cfg)
        np.testing.assert_allclose(read_raw, raw, equal_nan=True)
        for scenario in [{'method':'knn','scope':'pooled'}, {'method':'knn','neighbors':'samples'}]:
            complete, indicator = p.impute(retained, groups, scenario)
            assert np.array_equal(complete.values[~indicator], retained.values[~indicator])
        linear = annotations.join(np.exp2(raw)).reset_index(names='feature_id')
        linear.loc[0,'s0'] = 0
        linear.to_csv(root/'linear.csv', index=False)
        ingestion_cfg['input'] = {'path':str(root/'linear.csv'),'scale':'linear'}
        lin, _, _, _ = p.read_inputs(ingestion_cfg)
        assert pd.isna(lin.loc['f0','s0'])
        long = np.exp2(raw).reset_index(names='Genes').melt(id_vars='Genes', var_name='Run', value_name='Genes.MaxLFQ')
        long.to_csv(root/'diann.csv',index=False)
        ingestion_cfg['input'] = {'path':str(root/'diann.csv'),'format':'diann_protein','scale':'linear','aggregation':'mean'}
        dia, _, _, _ = p.read_inputs(ingestion_cfg)
        np.testing.assert_allclose(dia.reindex(raw.index), raw, equal_nan=True)
        cfg = {'input': {'path': str(root/'matrix.csv'), 'format': 'matrix', 'scale': 'log2'}, 'metadata_path': str(root/'metadata.csv'), 'output_dir': str(root/'results'), 'contrasts': [unpaired, paired], 'formats': ['png'], 'imputations': [{'name':'knn', 'method':'knn'}, {'name':'gaussian','method':'gaussian'}], 'enrichment': {'gmt_path':str(root/'sets.gmt')}}
        summary = run(cfg)
        assert summary['status'] == 'complete'
        assert len(summary['scopes']) == 3
        hit_audit = pd.read_csv(root/'results/abundance/paired/knn/imputation_method_audit.csv', index_col=0)
        assert not hit_audit.loc['f89','limma_testable']
        for scope in ['all_samples', 'paired', 'unpaired']:
            for scenario in ['knn', 'gaussian']:
                dest = root / 'results/abundance' / scope / scenario
                assert (dest/'pca.png').exists() and (dest/'imputation_mask.csv').exists()
        assert (root/'results/abundance/paired/hallmark_up_top10.png').exists()
        # Known ORA with duplicate sites: background contains unique genes, not site count.
        oraout = root/'ora'
        oraout.mkdir()
        ora_res = pd.DataFrame({'pvalue':[.001,.001,.5], 'padj':[.002,.002,.5], 'log2FoldChange':[2,3,0]}, index=['s1','s2','s3'])
        ora_ann = pd.DataFrame({'gene_symbol':['ONE','ONE','TWO']}, index=ora_res.index)
        p.enrichment(ora_res, ora_ann, {'TERM':{'ONE'}}, oraout, .05, ['png'])
        ora = pd.read_csv(oraout/'hallmark_up_all.csv')
        assert ora.loc[0,'background_genes'] == 2 and ora.loc[0,'directional_hit_genes'] == 1
        # Protein subtraction uses original NAs; no imputed value becomes an observation.
        proteins = pd.DataFrame(19., index=[f'P{i}' for i in range(90)], columns=md.index)
        proteins.loc['P0', 's0'] = np.nan
        proteins.reset_index(names='protein_id').to_csv(root/'protein.csv', index=False)
        adjusted, adjusted_ann, report = p.adjust_phosphosites(raw, annotations, {'protein_adjustment':{'path':str(root/'protein.csv'),'feature_id_column':'protein_id','scale':'log2'}})
        assert pd.isna(adjusted.loc['f0','s0'])
        np.testing.assert_allclose(adjusted.loc['f1'], raw.loc['f1']-19, equal_nan=True)
        print('PASS: end-to-end paired/unpaired KNN/Gaussian, real limma, PCA/ORA/audits, NA handling, protein adjustment and confounding checks')


if __name__ == '__main__':
    check()
