# Inputs and configuration

Contents: [Inputs](#required-data), [formats](#formats), [Excel facility workflow](#facility-workbooks), [scopes](#preprocessing-and-scopes), [imputation](#imputation-scenarios), [comparisons](#comparisons-pairing-and-limma), [protein adjustment](#optional-matched-protein-adjustment), [enrichment](#enrichment-and-visualization), [execution](#execution-and-deliverables).

## Required data

Supply a feature-by-sample CSV/TSV/Excel sheet with unique `feature_id`, a `gene_symbol` annotation and sample columns matching metadata. A phosphosite ID must preserve its unique site identity (e.g. `AKT1_S473`); retain protein accession, residue, position and localization columns as annotations when available. Annotation names must not collide with exported statistical fields. The default is a human analysis; supply species-matched symbols and a local GMT for other species.

Metadata needs unique `sample_id` and the comparison factor (default `condition`); a `subject_id` is required only for paired analysis. Additional categorical columns can specify tissue, batch, etc. Use one biological sample per row. Pairing cannot be inferred reliably from intensity values. Matrix columns absent from metadata are reported as excluded when explicitly included in `input.sample_columns`; metadata samples absent from the matrix cause an error. By default only exact metadata sample columns are read and other columns remain annotations.

Set `input.scale` explicitly to `linear` or `log2`. Zeros become NaN only on the linear scale. Invalid numeric text, negative linear intensities, infinities, duplicate IDs and wholly absent samples raise errors. Do not encode missing log2 values with zero. Supply empty cells/NaN instead.

## Formats

- `matrix`: CSV/TSV/Excel sheet with feature ID, gene annotation and sample columns. Configure `feature_id_column`, `gene_column`, and optional `sample_columns`. Default names are `feature_id` and `gene_symbol`.
- `diann_protein`: long CSV/TSV/parquet with `Run`, `Genes`, `Genes.MaxLFQ`. Configure `run_column`, `gene_column`, `value_column` as necessary. Set `scale: linear`, and explicitly choose `aggregation: mean` (notebook behavior) or `median`. For standard parquet reports, read only run/gene/intensity/protein-ID columns. When `Protein.Ids` is available, first aggregate repeated precursor rows per gene/protein/run, then equally aggregate protein groups per gene/run, matching notebook weighting. Without that column aggregate per gene/run and record this difference. Retain joined protein accessions as annotations. Configure `protein_column` when needed. Use optional `sample_mapping` from original Run strings to metadata IDs. No implicit run-path or extension stripping. Exclude contaminants/decoys using documented upstream filters and record them; the adapter does not infer such labels. Parquet requires pyarrow.
- `diann_phosphosite`: wide DIA-NN phosphosite CSV/TSV with `Gene.Names`, `Residue`, `Site` and exact sample columns. Creates gene_residue-position IDs; rejects duplicates. Configure column names when needed.

- `facility_proteins`: `.xlsx`/`.xlsm` or CSV/TSV containing the facility's protein-level `Summary` table. Preserve `UNIPROT_AC`, `GENE_NAMES`, `FASTA_HEADERS`, `PROTEIN_DESCRIPTIONS`, `NUM_PRECURSORS`, `PROTEIN_WEIGHT`, and `NUM_UNIQUE_PEPTIDES`. Use `UNIPROT_AC` as the unique feature ID unless `feature_id_column` is supplied. Default sheet is exactly `Summary`; select measurements only from columns starting with `Raw_Intensities`. Remove that prefix to obtain metadata sample IDs, e.g. `Raw_Intensities011_MP002_B26` maps to `011_MP002_B26`. Preserve leading zeros. Counts/weights remain annotations. Configure `intensity_prefix` only for a deliberately different naming convention. Explicit `sample_mapping` entries override prefix removal and must reference full original raw-intensity column names. Set `scale: linear`. Prefer this existing protein-level summary for protein-level analyses.
- `facility_peptides`: multi-sheet `.xlsx`/`.xlsm`, or a CSV/TSV export, with `PEPTIDE_SEQ`, `PROTEIN_UNIPROT_AC`, `FASTA_HEADERS`, `GENE_NAMES`, `PROTEIN_DESCRIPTIONS`, and raw-intensity sample columns. Excel defaults to the exact `Peptides (raw)` sheet; if absent, require an explicit `sheet_name`. Preserve all five source annotation columns plus canonical `gene_symbol` and `protein_id`. By default form feature IDs as `PEPTIDE_SEQ|PROTEIN_UNIPROT_AC`; alternatively supply an existing unique `feature_id_column`. Repeated sequences assigned to different protein groups remain separate source rows; duplicate full feature IDs cause an error rather than averaging. Preserve shared/multi-gene mappings, with unresolved genes omitted from ORA. This adapter produces peptide-level data and performs no automatic protein summarization.

## Facility workbooks

Inspect a workbook before preparing its configuration:

```bash
python /absolute/skill/path/scripts/run_proteomics.py --inspect facility.xlsx
```

This inventories every sheet and its column names, identifies both facility protein and peptide schemas, and lists candidate sample columns without treating them as confirmed sample identities. The same command inspects a parquet schema without loading the full report. For a workbook with a title block, add `--header-row 2` (zero-based) and configure the same `input.header_row` for analysis. Workbooks may contain README, Summary, Identification Qvalue/PEP, Spectral Counts and IBAQ Values tabs. `Summary` contains protein raw intensities when its `Raw_Intensities` columns and protein annotations match the supplied schema. Use that table for protein analyses. Other score/count/IBAQ tabs are not substitutes for the requested raw-intensity measurements. Do not import a score/count table as intensities or combine sheets automatically.

For protein-level analysis, use this block with the ordinary metadata/contrast configuration:

```json
"input": {
  "path": "facility.xlsx",
  "format": "facility_proteins",
  "sheet_name": "Summary",
  "scale": "linear"
}
```

Metadata sample IDs should be `011_MP002_B26`, etc. To map a different metadata convention, add e.g. `"sample_mapping": {"Raw_Intensities011_MP002_B26": "control_01"}`; other raw columns still use prefix removal. Only raw-intensity columns can be selected or mapped as protein samples. Unused raw-intensity columns are reported explicitly. Duplicate protein IDs require an explicit unique feature ID or a documented resolution; they are never silently combined.

For peptide-level analysis, use this block:

```json
"input": {
  "path": "facility.xlsx",
  "format": "facility_peptides",
  "sheet_name": "Peptides (raw)",
  "scale": "linear",
  "sample_mapping": {"Facility run 01": "sample_01", "Facility run 02": "sample_02"}
}
```

Omit `sample_mapping` when workbook sample names match metadata exactly. The mapping keys are original measurement-column names, and values are metadata sample IDs; an explicit `sample_columns` list uses mapped names. Reject missing source columns, many-to-one mappings and mappings that rename annotation columns. All measurements used in the analysis must have metadata. Other unselected columns are preserved as annotations and listed in the input report; inspect them for omitted samples. Never select every numeric column automatically.

`matrix` and `diann_phosphosite` also support explicit Excel `sheet_name`; a multi-sheet generic workbook without one errors rather than reading the first tab. Metadata can itself be Excel with `metadata_sheet_name`. Optional matched-protein Excel input uses `protein_adjustment.sheet_name` and must contain a protein-level abundance matrix, not this peptide table. Formula cells require valid cached numeric values: do not fill missing formula caches with zero. Raw facility peptide intensities may require normalization; preserve scale and review the facility README/methods before choosing `normalization: none` or `median`.

Ask for the intended analysis level when unclear. Use the existing `Summary` protein raw intensities for protein-level work; ask before converting peptide measurements to protein abundance. No peptide-to-protein aggregation or MaxLFQ estimation is inferred. Do not label peptide-level tests as protein-level changes or infer localized phosphosites from peptide annotations alone. Other protein/site sheets can be read as `matrix` once the appropriate identifiers, gene annotations, measurement meaning and sample columns are explicitly configured.

For phosphosite localization, set `input.localization_column`; `localization_threshold` defaults to 0.9 and expects probabilities on 0..1. Rows lacking a probability fail the filter. Without that column no localization filtering is claimed. A prefiltered input can be used, but document its upstream threshold rather than guessing from its filename.

## Preprocessing and scopes

`normalization: none` preserves already normalized LFQ values (notebook behavior). Optional `median` subtracts each sample's observed log2 median relative to the scope-wide median of medians. This is performed independently within each scope; do not use it without considering global biological shifts.

`filtering.min_total: 4` and `min_per_group: 1` count observed values before any imputation. `global_group_columns` defaults to `[condition_column]`; use `["tissue", "condition"]` to reproduce separate tissue/condition fits in a combined dataset. This grouping determines both the group presence filter and within-group KNN partitions. If a feature should be retained despite being entirely absent from a group, explicitly change `min_per_group`; do not silently change it. Feature-neighbor KNN can borrow from other proteins/sites in that group's sample columns; sample-neighbor KNN cannot infer an entirely missing feature within a group and will error. Entirely missing sample columns after filtering also error.

Global and contrast processing each start from observed input, so a feature may be retained for one comparison but absent globally. Comparison PCA uses only that comparison's samples and its own retained features, imputation and feature ranking. It is not a view of the global fitted PCA. Coordinate axes from separate PCAs are not directly comparable.

## Imputation scenarios

`imputations` is a list; each entry gets its own output subdirectory. Give every entry a distinct filesystem-safe `name`. Defaults:

```json
{"name":"knn_within_condition","method":"knn","scope":"within_group","neighbors":"features","k":3,"weights":"uniform"}
```

`scope: pooled` fits across all samples of the current scope. `neighbors: features` reproduces notebook matrix orientation; `samples` transposes before fitting, for a deliberately different sensitivity analysis. KNN uses sklearn's nan-euclidean distance and excludes missingness indicators from the expression matrix.

Optional Gaussian scenario:

```json
{"name":"gaussian","method":"gaussian","scope":"within_group","shift":1.8,"width":0.3}
```

For every missing entry in each sample, draw from N(observed sample log2 mean - shift*sample SD, (width*sample SD)^2). `seed` defaults to 42; sample SD uses ddof=1. This sample-wise model gives the same distributions for pooled/within-group settings; the KNN grouping distinction does not create a new Gaussian algorithm. Require at least two observed values and a positive SD when filling is needed. This is a left-censoring heuristic; it is not a missingness mechanism inferred from the data. Each scenario is deterministic for the same ordered inputs/configuration/seed and preserves observed measurements exactly.

## Comparisons, pairing and limma

An unpaired contrast:

```json
{"name":"drug_liver","factor":"condition","test":"drug","reference":"vehicle","paired":false,"subset":{"tissue":["liver"]}}
```

A paired contrast:

```json
{"name":"resection_vs_biopsy","factor":"collection","test":"Resection","reference":"Biopsy","paired":true,"subject_column":"patient_id","subset":{"tissue":["tumor"]}}
```

`paired` defaults to false in code; explicitly set it based on study design when preparing configurations. Paired comparisons exclude metadata-incomplete subjects and report those sample IDs, then use that same sample set for comparison PCA, t-test and limma. Duplicate subject/condition measurements error. This is a balanced two-condition paired workflow, not a mixed model.

The imputed t-test is paired only on request; the unpaired default is equal-variance Student's, not Welch's. Non-imputed limma always runs. Paired limma works on original within-subject differences; default `limma.min_complete_pairs: 4` applies per feature. Unpaired limma requires `limma.min_per_group: 2`, feature-wise full-rank observed design and positive residual df. Default `trend: true`, `robust: false` reproduce the executed notebook. Robust InMoose moderation is not supported; requesting it errors rather than silently disabling it.

For unpaired limma, add `covariates: ["batch", "age"]` and `numeric_covariates: ["age"]` to a contrast. Other covariates are categorical. The t-test remains unadjusted, and the audit must acknowledge this additional model difference. Stop if variables are confounded. Covariates in the paired-difference workflow require a separate model design and are rejected, not silently ignored.

## Optional matched-protein adjustment

Set `assay: phosphoproteomics` and add:

```json
"protein_adjustment": {
  "path":"matched_protein.csv",
  "feature_id_column":"protein_id",
  "mapping_column":"protein_id",
  "scale":"linear"
}
```

The phosphosite annotations must contain an explicit, unambiguous `protein_id` matching a unique row in the protein matrix. Sample columns must match exactly. Supply consistently normalized input modalities; no separate matched-protein median normalization is inferred. The code subtracts protein log2 intensities from site log2 intensities using observed values only. Missing proteins create missing adjusted sites. Unmapped sites are counted and omitted from the adjusted analysis. The unadjusted abundance analysis is also retained. Never select one member of an ambiguous protein group automatically.

## Enrichment and visualization

`enrichment` defaults to human Hallmark `2024.1.Hs`, matching the transcriptomics reference. By default download the official GMT and record its SHA256. For offline use set `gmt_path` and its release; use `enabled: false` only when explicitly opting out. Retrieval errors fail the run with a diagnostic; enrichment is never silently omitted. Do not bundle downloaded MSigDB material inside the skill.

`alpha: 0.05`, `pca_top_features: 500`, `pca_labels: false`, `formats: ["pdf", "png"]`, optional `palette` category-to-color mapping. Set matching condition colors across assays/scopes; categories missing from the mapping use Matplotlib's color cycle. Export feature scores/loadings/variance for every PCA. Statistical results include all retained features, with untestable limma rows retained and explained.

## Execution and deliverables

Install `scripts/requirements.txt`; for source builds without clang, use available compilers, e.g. `CC=gcc CXX=g++ python -m pip install -r scripts/requirements.txt`. Excel reading requires openpyxl>=3.1 and parquet reading requires pyarrow>=14 (both included in requirements.txt). No R installation is needed. Require the real InMoose backend; do not emulate empirical-Bayes variance moderation.

Run `python /absolute/skill/path/scripts/run_proteomics.py --config analysis.json`. Paths resolve relative to that JSON. Use an empty `output_dir` to avoid mixing runs.

The root contains observed input, feature annotations, metadata, configuration, methods, hashes/versions and summary. Each `assay/scope` has filtering/normalization logs and observed non-imputed limma results. Each scenario has its own imputed matrix, binary mask, missingness summaries, histograms, correlation comparison, PCA and (for contrast scopes) t-test, volcano, Hallmark and imputation-impact audit. `status: failed` means outputs are partial and must not be described as a complete run. When only a partial stage is requested, import its function instead of invoking the whole runner.
