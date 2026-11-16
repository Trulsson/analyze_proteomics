# Methods and interpretation

## Provenance and intentional changes

Procedures were extracted from `Trulsson/biopsy_vs_resection` on 2026-11-02:

- Proteomics notebook `02_proteomics_biopsy_resection_clean.ipynb`, Git blob `3f82fbdabc039c0e3926010005a582dde27743d2`.
- Phosphoproteomics notebook `03_phosphoproteomics_biopsy_resection_clean.ipynb`, Git blob `667c1425fdb714e8337d10998da3ee565775290f`.
- Directional Hallmark output conventions follow the installed `analyze-bulk-transcriptomics` skill.

Preserve log2 intensities, minimum four observed measurements overall/one per condition, within-condition three-neighbor KNN with features as neighbors, top-500 variable-feature scaled PCA, imputed Student's t-tests, and non-imputed paired-difference limma with at least four complete pairs. The source helper declares `robust=True` but executes `eBayes(... robust=False)`; use the executed behavior and reject unsupported robust requests.

Intentional changes requested by the user: pairing is data-dependent and defaults to unpaired in code; reprocess all-sample and comparison-only scopes independently from observed input; include optional Gaussian imputation, unpaired limma and parallel protein-adjusted phosphosite analysis; use directional Hallmark outputs for both assays. Do not copy the source's cohort exclusions. Do not include transcriptomics, PCA distances, permutations or PTMsigDB/kinase analyses.

The notebook's combined PCA concatenates separately filtered tissue matrices and performs another KNN pass. This skill deliberately does not reproduce that step: global PCA starts from the original full dataset, filters/imputes on that scope, and fits its own PCA. Compare subset structure visually, without treating independent PCA axes as the same coordinate system.

## Input formats and feature identity

DIA-NN parquet import reproduces the notebook's two-stage aggregation when `Protein.Ids` is available: collapse duplicate precursor-level gene/protein/run intensities, then equally average unique protein groups within each gene/run. This avoids weighting one protein group by its number of report rows. Read only the required parquet columns; preserve joined protein IDs and report the exact aggregation. Without a protein-ID column, record the single-stage gene/run aggregation explicitly.

Facility Excel import selects a named measurement sheet. For protein-level `Summary`, use `UNIPROT_AC` as the feature identifier and only columns with the `Raw_Intensities` prefix as sample measurements; remove that prefix to match metadata, preserving leading zeros. Preserve the seven protein annotation columns, including numeric precursor counts, molecular weight and unique-peptide counts, outside the intensity matrix. Do not derive protein abundance from peptides when the facility has supplied the requested protein-level intensities. Peptide-level import selects `Peptides (raw)` and reads only metadata-matched raw sample columns into the numerical matrix. Preserve original peptide sequences, accessions, FASTA headers, gene names and descriptions; keep composite peptide/protein feature IDs separate from gene symbols. Preserve multi-gene mappings in original annotations without expanding them into independent enrichment genes. Peptide-level measurements remain peptide-level; protein summarization and phosphosite localization cannot be inferred from these columns. Record the selected sheet, feature level, sample-name mapping and available sheets in the input report.

## Missing values and PCA

Intensity zero in linear-scale input means missing; log2 zero means an observed intensity of one and must remain observed. Keep observed log2 values and missingness indicators separate. Optional median normalization operates on observed values before filtering in the current scope.

KNN takes a feature x sample matrix exactly as in the source. Thus its neighbors are other proteins/sites with similar profiles in the selected sample columns. A within-condition fit excludes the opposite condition's columns when finding neighbors, but still uses condition labels in preprocessing. It can alter group separation, so PCA is descriptive rather than independent evidence of treatment effects. Pooled KNN and sample-neighbor KNN are explicit alternatives. No imputation indicator columns enter PCA or tests.

Optional Gaussian imputation samples each missing measurement from a lower-intensity distribution derived from observed values in that sample, with mean shifted down 1.8 SD and width 0.3 SD by default. It assumes missing low intensities rather than demonstrating this mechanism. Record its seed and assess whether it creates condition-specific effects. All methods preserve observed measurements exactly.

Select up to 500 features with greatest log2 SD after scope-specific imputation; remove zero-variance features, standardize each feature across the scope's samples, then compute PC1/PC2 with explained variance. Export scores, loadings and selected feature IDs. Each scenario may select different variable features; this is part of its full preprocessing effect, not a comparison holding the feature set fixed.

## Differential abundance

For each requested comparison, use only its subset's samples. Positive log2FC always means test minus reference. For unpaired data, run a two-sided equal-variance Student t-test on imputed log2 values, matching the source helper; do not silently substitute Welch's test. For paired data, align subject IDs and run a paired t-test on imputed differences. Biological replicates, not precursor/technical repetitions, define sample size.

Run non-imputed InMoose limma separately. In paired analyses, compute an observed test-minus-reference difference per subject, with NaN whenever either member is missing. Fit an intercept-only empirical-Bayes model for features with >=4 complete pairs and positive residual df. Supply the original mean log2 abundance explicitly to the variance-trend fit: mean difference is not an abundance covariate. Metadata-incomplete subjects are excluded and reported before comparison processing.

For unpaired analyses, fit intercept + test indicator + requested covariates, require >=2 observations in each group by default, a full-rank feature-wise observed design and positive residual df. Categorical and numeric covariates must be specified explicitly. A covariate-adjusted limma comparison versus an unadjusted t-test includes a model difference; the paired-difference implementation rejects covariates rather than silently omitting them.

Use `inmoose==0.9.1` model fitting and empirical-Bayes `squeezeVar` with observed mean abundance as the trend covariate and `robust=False`; compute moderated t-statistics, degrees of freedom and two-sided p-values with the same formulas as InMoose `_ebayes`, then BH-adjust finite p-values. B-statistics are not computed. This avoids an upstream InMoose 0.9.1 B-statistic failure when the estimated prior df is scalar infinity; the variance moderation and t-test remain unchanged. Keep non-testable features in the output with NA statistics and a reason; do not count them as evidence of no effect. Adjust finite t-test p-values using BH separately; do not let an undefined t-statistic propagate NA into every adjusted p-value, as can happen with the notebook's hand-written adjustment. This is a correctness fix for undefined/constant features, not a different multiple-testing method.

## Enrichment

Use MSigDB Hallmark `2024.1.Hs` by default and record the downloaded GMT hash; alternatively use a supplied species-matched GMT. Run hypergeometric over-representation separately on significant positive and negative log2FC mapped genes. Deduplicate gene symbols within each direction and restrict gene sets to the method-specific testable mapped background. BH-adjust across all supplied gene sets, including zero-overlap terms; this makes the tested hypothesis family explicit. This differs from tools that suppress zero-overlap terms before multiple-testing adjustment and is recorded rather than labeled numerically identical to their outputs.

For phosphosites, retain distinct sites in the statistical tables and contributor plots, but count their gene once per directional ORA. One gene can have sites changing in both directions and appear in both analyses; do not infer a single direction for its protein activity. Exclude unresolved multi-gene protein groups from ORA rather than assigning the same protein signal to multiple independent genes. Record unmapped/ambiguous counts. Hallmark describes gene-level biological context and is not a test of kinase activity or phosphosite-specific regulation.

Always export full and top-ten tables and plots, with an explicit FDR flag and threshold. Include nominal and adjusted p-values, unique background size, set size in background, hit count, overlap and contributing genes. Add significant directional Hallmark membership to result tables. Use the same principle for both non-imputed limma and each imputed t-test scenario. Empty hit lists produce p=1 tables and an explanatory contributor placeholder.

## Protein-adjusted phosphosites

Retain the unadjusted analysis. When an explicit unique protein mapping and matched samples are available, subtract observed matched-protein log2 intensity from observed phosphosite log2 intensity, then independently process that adjusted matrix through the full workflow. Missing protein measurements remain missing and contribute to the adjusted-data imputation audit. Do not impute protein values first or interpret the resulting ratio as absolute occupancy. Differences in modality normalization, ambiguous protein groups and isoform/site mappings constrain interpretation.

## Imputation-impact reporting

Export sample/group/feature observed and imputed counts and fractions, pagination of observed-versus-imputed intensity histograms, pairwise-complete observed and imputed correlation matrices, and the imputation mask. For t-test hits, summarize the number with zero imputation, number with imputation, median imputed entries among imputed hits, number testable by limma, number significant by limma and direction agreement among testable hits. Keep untestable limma hits separate from testable non-significant hits. Plot log2FC agreement colored by imputed fraction, hit overlap categories, and imputation burden among hits versus other features.

These diagnostics quantify missingness and robustness to a different analysis; they do not estimate the true unobserved intensities. The t-test/limma comparison changes both imputation and modeling/eligibility, so it does not isolate a causal effect of imputation. Missing-not-at-random mechanisms, small sample size, normalization shifts, localization uncertainty and covariate confounding remain limitations. Do not call a plausible imputed distribution proof of accuracy.

## Primary technical references

- Source: https://github.com/Trulsson/biopsy_vs_resection
- KNN row-neighbor semantics: https://scikit-learn.org/stable/modules/generated/sklearn.impute.KNNImputer.html
- InMoose lmFit (missing observations): https://inmoose.readthedocs.io/en/stable/generated/inmoose.limma.lmFit.html
- InMoose empirical Bayes: https://inmoose.readthedocs.io/en/stable/generated/inmoose.limma.eBayes.html
- InMoose topTable: https://inmoose.readthedocs.io/en/stable/generated/inmoose.limma.topTable.html
- Hallmark: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
