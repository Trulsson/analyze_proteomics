"""Regression check for facility Excel and notebook-style DIA-NN parquet input adapters."""
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
import pandas as pd
import proteomics as p


def check():
    fixture = Path(__file__).resolve().parents[1] / 'assets/test_inputs/facility.xlsx'
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        md = pd.DataFrame({'sample_id': [f's{i}' for i in range(1,9)], 'condition':['A']*4+['B']*4})
        md.to_csv(root/'md.csv', index=False)
        sample_mapping = {f'Run {i}': f's{i}' for i in range(1,9)}
        config = {'input': {'path':str(fixture), 'format':'facility_peptides', 'scale':'linear', 'sample_mapping':sample_mapping}, 'metadata_path':str(root/'md.csv')}
        inventory = p.inspect_input(fixture)
        assert len(inventory['sheets']) == 7
        assert inventory['sheets']['Peptides (raw)']['suggested_format'] == 'facility_peptides'
        assert inventory['sheets']['Summary']['suggested_format'] == 'facility_proteins'
        assert len(inventory['sheets']['Summary']['candidate_sample_columns']) == 8
        protein_md = md.copy()
        protein_md['sample_id'] = [f'{i:03d}_MP002_B26' for i in range(11,19)]
        protein_md.to_csv(root/'protein_md.csv', index=False)
        protein_config = {'input':{'path':str(fixture), 'format':'facility_proteins', 'scale':'linear'}, 'metadata_path':str(root/'protein_md.csv')}
        protein_raw, protein_ann, _, protein_report = p.read_inputs(protein_config)
        assert protein_raw.shape == (32,8) and protein_report['feature_level']=='protein'
        assert protein_report['sheet_name']=='Summary'
        assert set(p.FACILITY_PROTEIN_ANNOTATIONS).issubset(protein_ann.columns)
        assert list(protein_raw.columns)==protein_md.sample_id.tolist()
        assert protein_report['effective_sample_mapping']['Raw_Intensities011_MP002_B26']=='011_MP002_B26'
        assert protein_ann.at['P0','NUM_PRECURSORS']==12
        assert 'NUM_PRECURSORS' not in protein_raw.columns
        assert pd.isna(protein_raw.at['P0','011_MP002_B26'])
        np.testing.assert_allclose(protein_raw.at['P0','012_MP002_B26'],np.log2(2340))
        renamed_md = protein_md.copy()
        renamed_md.loc[0,'sample_id'] = 'control_01'
        renamed_md.to_csv(root/'renamed_md.csv', index=False)
        override = dict(protein_config, metadata_path=str(root/'renamed_md.csv'), input=dict(protein_config['input'], sample_mapping={'Raw_Intensities011_MP002_B26':'control_01'}))
        renamed_raw, _, _, _ = p.read_inputs(override)
        assert list(renamed_raw.columns) == renamed_md.sample_id.tolist()
        try:
            p.read_inputs(dict(protein_config, input=dict(protein_config['input'], sample_columns=['NUM_PRECURSORS',*protein_md.sample_id])))
        except ValueError as exc:
            assert 'not annotation' in str(exc)
        else:
            raise AssertionError('A numeric protein annotation became a sample intensity')
        raw, ann, metadata, report = p.read_inputs(config)
        assert raw.shape == (32,8) and report['feature_level'] == 'peptide'
        assert report['sheet_name'] == 'Peptides (raw)'
        assert set(p.FACILITY_ANNOTATIONS).issubset(ann.columns)
        assert ann.loc['PEP0|P0','GENE_NAMES'] == 'G0'
        assert set(['PEP0|P0','PEP0|P1']) <= set(raw.index)
        assert pd.isna(raw.loc['PEP0|P0','s1'])
        np.testing.assert_allclose(raw.loc['PEP0|P0','s2'],np.log2(2340))
        assert pd.isna(ann.loc['PEP31|P31','gene_symbol'])
        assert ann.loc['PEP31|P31','GENE_NAMES'] == 'G30;G31'
        table = p.read_table(fixture, sheet_name='Peptides (raw)')
        table.to_csv(root/'peptides.csv', index=False)
        csv_config = dict(config, input=dict(config['input'], path=str(root/'peptides.csv')))
        csvraw, csvann, _, _ = p.read_inputs(csv_config)
        pd.testing.assert_frame_equal(raw, csvraw)
        # Generic multi-sheet reads never silently choose README/the first sheet.
        try:
            p.read_table(fixture)
        except ValueError as exc:
            assert 'sheet_name' in str(exc) and 'Peptides (raw)' in str(exc)
        else:
            raise AssertionError('Ambiguous worksheet was selected silently')
        wrong = dict(config, input=dict(config['input'],sheet_name='IBAQ Values'))
        try:
            p.read_inputs(wrong)
        except ValueError as exc:
            assert 'lacks columns' in str(exc)
        else:
            raise AssertionError('Wrong sheet was interpreted as peptide measurements')
        bad = pd.concat([table,table.iloc[:1]], ignore_index=True)
        bad.to_csv(root/'duplicates.csv',index=False)
        try:
            p.read_inputs(dict(config,input=dict(config['input'],path=str(root/'duplicates.csv'))))
        except ValueError as exc:
            assert 'Duplicate feature' in str(exc)
        else:
            raise AssertionError('Duplicate peptide/protein rows were silently aggregated')
        # Unequal precursor row counts must not weight one protein more strongly.
        dia = pd.DataFrame([
            {'Run': f'run{i}', 'Genes':'ONE', 'Protein.Ids':protein, 'Genes.MaxLFQ':value, 'Unused':'ignore'}
            for i in range(1,9)
            for protein,value in [('P1',4.),('P1',4.),('P1',4.),('P2',16.)]
        ])
        dia.to_parquet(root/'report.parquet',index=False)
        assert p.inspect_input(root/'report.parquet')['suggested_format']=='diann_protein'
        diaconfig = {'input':{'path':str(root/'report.parquet'),'format':'diann_protein','scale':'linear','aggregation':'mean','sample_mapping':{f'run{i}':f's{i}' for i in range(1,9)}},'metadata_path':str(root/'md.csv')}
        values, annotations, _, report = p.read_inputs(diaconfig)
        np.testing.assert_allclose(values.loc['ONE'],np.log2(10.))
        assert annotations.at['ONE','protein_id']=='P1;P2'
        assert 'across proteins' in report['aggregation']
        assert 'Unused' not in annotations.columns
        collision = dict(diaconfig['input'], sample_mapping={f'run{i}': ('s1' if i==2 else f's{i}') for i in range(1,9)})
        try:
            p.read_inputs(dict(diaconfig, input=collision))
        except ValueError as exc:
            assert 'merges distinct' in str(exc)
        else:
            raise AssertionError('Many-to-one run mapping silently merged biological samples')
        print('PASS: multi-sheet Excel selection, peptide annotations/identity, zero handling, sample mapping, wrong-sheet/duplicate rejection, and two-stage DIA-NN parquet aggregation, Summary proteins and Raw_Intensities prefix mapping')


if __name__ == '__main__':
    check()
