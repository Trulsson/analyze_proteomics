# Analyze Proteomics

**A reusable analysis skill for quantitative bulk proteomics and phosphoproteomics.**

`analyze-proteomics` provides configurable Python workflows for preprocessing, exploratory analysis, differential abundance testing, pathway enrichment, and missing-value diagnostics. It can be used as an AI-assistant skill through its [`SKILL.md`](SKILL.md) instructions or run directly with the bundled command-line tools.

The statistical approach is adapted from the proteomics and phosphoproteomics notebooks in [Trulsson/biopsy_vs_resection](https://github.com/Trulsson/biopsy_vs_resection), with extensions for general experimental designs and additional input formats.

## Features

| Area | Capabilities |
| --- | --- |
| **Input formats** | Feature-by-sample CSV/TSV/Excel matrices; DIA-NN protein parquet/CSV/TSV reports and phosphosite tables; facility Excel `Summary` and `Peptides (raw)` sheets |
| **Preprocessing** | Explicit linear/log2 intensity handling, optional median normalization, observed-value filtering, and sample metadata validation |
| **Missing values** | Within-condition KNN imputation (default, *k* = 3); optional pooled KNN or Gaussian left-shifted imputation; retained observed-value and imputation masks |
| **Exploratory analysis** | All-sample and contrast-specific PCA, each independently preprocessed from the observed matrix; missingness and intensity QC |
| **Differential abundance** | Paired or unpaired t-tests on imputed intensities; non-imputed empirical-Bayes limma using InMoose; Benjamini–Hochberg FDR correction |
| **Functional enrichment** | Separate up- and down-regulated MSigDB Hallmark over-representation analyses with contributor tables and figures |
| **Sensitivity and diagnostics** | Comparison of t-test and limma results, hit-level imputation burden, effect-size agreement, and missingness visualizations |
| **Phosphoproteomics** | Site-level statistics with preserved site identities; optional analysis of phosphosite intensities adjusted by matched protein abundances |
| **Reporting** | Annotated CSVs, PDF/PNG figures, configuration and method records, software versions, and a machine-readable run summary |

## Repository layout

```text
analyze_proteomics/
├── SKILL.md                       # Skill instructions and analytical guardrails
├── agents/openai.yaml             # Assistant-facing skill metadata
├── assets/test_inputs/            # Example facility input for testing
├── references/
│   ├── example-config.json        # Starting configuration
│   ├── input-and-options.md       # Input schemas and configuration reference
│   └── methods.md                 # Statistical methods and interpretation
└── scripts/
    ├── proteomics.py              # Reusable analysis functions
    ├── run_proteomics.py          # Configuration-driven command-line runner
    ├── requirements.txt           # Python dependencies
    ├── test_inputs.py             # Input adapter checks
    └── test_proteomics.py         # Synthetic workflow regression tests
```

## Getting started

### 1. Install the dependencies

Clone the repository and install the bundled requirements into a dedicated Python environment:

```bash
git clone https://github.com/Trulsson/analyze_proteomics.git
cd analyze_proteomics
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
python -m pip install -r scripts/requirements.txt
```

The pipeline uses NumPy, pandas, SciPy, scikit-learn, matplotlib, InMoose, openpyxl, and pyarrow. Installing InMoose may require a C/C++ compiler. An R installation is not required.

### 2. Inspect your input data

```bash
python scripts/run_proteomics.py --inspect path/to/measurements.xlsx
```

The inspection command reports worksheets, columns, and candidate intensity fields so that the appropriate input format and biological samples can be selected explicitly. It also supports supported CSV/TSV and parquet inputs.

### 3. Prepare an analysis configuration

Start with [`references/example-config.json`](references/example-config.json):

```bash
cp references/example-config.json analysis.json
```

Edit the paths, input format, intensity scale, sample metadata, contrasts, and output directory. For example, a facility workbook with protein intensities in its `Summary` sheet uses:

```json
"input": {
  "path": "facility.xlsx",
  "format": "facility_proteins",
  "sheet_name": "Summary",
  "scale": "linear"
}
```

A minimal sample metadata file includes one row per biological sample:

```csv
sample_id,condition,subject_id
sample_01,control,patient_01
sample_02,treated,patient_01
sample_03,control,patient_02
sample_04,treated,patient_02
```

These rows illustrate the schema only; real analyses require enough replicates for the configured filtering and statistical models. For matched samples, explicitly set `"paired": true` and `"subject_column": "subject_id"` in the contrast. For independent samples, set `"paired": false`. Pairing must not be inferred from row order.

### 4. Run the analysis

```bash
python scripts/run_proteomics.py --config analysis.json
```

Input and output paths in the JSON file are resolved relative to the configuration file. Use a **new, empty output directory** for every run.

After completion, check `results/run_summary.json`. Only `"status": "complete"` indicates that the full configured workflow finished successfully.

## Outputs

The results are organized by assay, analysis scope, and imputation scenario. Depending on the configuration, output files include:

- **PCA:** scores, loadings, explained variance, and publication-ready figures for all samples and each requested contrast.
- **Differential abundance:** annotated t-test and non-imputed limma results, log2 fold changes, raw and adjusted p-values, and volcano plots.
- **Pathway enrichment:** directional Hallmark over-representation results, top-pathway plots, and contributing-gene tables.
- **Imputation diagnostics:** observed/imputed intensity distributions, missingness summaries, imputation masks, and agreement between statistical methods.
- **Reproducibility records:** `analysis_config.json`, `methods.txt`, `run_summary.json`, input checksums, and dependency versions.

The all-sample PCA and each contrast-specific PCA are fitted **independently from the original observed measurements**, rather than subsetting a globally imputed matrix.

## Methodological considerations

- **Missing-data assumptions matter.** Imputation is a modeling choice, not a measurement of the missing abundance. The workflow preserves observed values and reports the impact of imputation on reported results.
- **Paired and unpaired designs are different.** Matching requires explicit subject identifiers. Unpaired t-tests use the equal-variance Student test; limma uses the relevant design and observed data without imputation.
- **Statistical agreement is not validation of imputation.** The imputed t-test and non-imputed limma also differ in feature eligibility and statistical modeling.
- **Phosphosite enrichment is gene-level context.** Hallmark enrichment does not establish kinase activity; protein-adjusted phosphosite values are not absolute phosphorylation occupancy.
- **Input processing is explicit.** Peptide measurements are not automatically summarized to protein abundance, and precursor counts or protein weights are not treated as intensity columns.

See [`references/methods.md`](references/methods.md) for implementation details and [`references/input-and-options.md`](references/input-and-options.md) for supported formats, configuration parameters, and statistical assumptions.

## Scope and limitations

This skill is intended for **quantitative bulk proteomics and phosphoproteomics analyses of processed intensity data**. It does not perform raw mass-spectrometry processing, RNA-seq analysis, single-cell proteomics, PCA distance metrics, permutation analyses, or kinase-activity inference. The bundled tests use synthetic/example data and do not substitute for biological validation.

## Origin and citation

The original analysis methods were developed for:

Heusler M-A, Nuciforo S, Trulsson F, *et al.* (2026). **Biopsy and Surgical Resection Yield Distinct Multi-Omic Profiles in Hepatocellular Carcinoma.** *The Journal of Pathology: Clinical Research*. [https://doi.org/10.1002/2056-4538.70121](https://doi.org/10.1002/2056-4538.70121)

This skill generalizes parts of that analytical workflow; it is **not** the archived code snapshot for the publication. The publication-specific notebooks are preserved at [Trulsson/biopsy_vs_resection](https://github.com/Trulsson/biopsy_vs_resection).

## Maintainer

[Trulsson](https://github.com/Trulsson)
