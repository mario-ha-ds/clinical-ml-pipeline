# Clinical ML pipeline: a comparative benchmark across two biomedical geometries

This project builds a robust, leak-free data mining pipeline to benchmark several machine learning algorithms across biomedical datasets with contrasting natures. The goal is to evaluate how different algorithms behave when faced with two opposing scenarios: a well-behaved dataset with clear morphological boundaries, and a highly complex, multifactorial dataset characterized by severe class overlap.

It is important to clarify that this is not an algorithm optimization project; the objective is not to squeeze out the highest possible evaluation metrics. The real purpose is to build a sound pipeline, benchmark diverse algorithmic families, and analyze their behavior across contrasting biomedical topologies, all while respecting the clinical asymmetry between a false negative and a false positive.

## The two scenarios

-   **Breast Cancer Wisconsin (Diagnostic):** 569 observations and 32 variables describing morphological properties of digitized cell nuclei, with a binary target (malignant/benign). Source: [UCI Machine Learning Repository](https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic).
-   **CAIR CVD 2025 (Cardiovascular Risk from Bangladesh):** 1,529 observations and 22 variables (demographic, lifestyle, and biochemical) describing cardiovascular risk, with a three-class target (low/intermediate/high). Source: [Kaggle](https://www.kaggle.com/datasets/jocelyndumlao/cair-cvd-2025-cardiovascular-risk-from-bangladesh/data).

These two datasets are deliberately chosen for their opposing geometries. The oncological dataset presents a largely linear, well-separated signal, while the cardiovascular dataset presents a dense, multifactorial continuum with heavy class overlap. Contrasting how each algorithmic family adapts to these two realities is the central thread running through every notebook.

## Model battery

We evaluate two main families of algorithms:

-   **Unsupervised (Notebook 02):** two centroid-based partitioning algorithms (K-Means under $L_2$ Euclidean distance and K-Medians under $L_1$ Manhattan distance), and a density-based algorithm (DBSCAN, guided by OPTICS reachability profiles). The objective here is to assess whether natural biological boundaries emerge prior to label exposure.
-   **Supervised (Notebook 03):** three distinct predictive paradigms tuned via 10-fold stratified cross-validation: a penalized (elastic net) Logistic Regression, a single cost-sensitive C5.0 decision tree, and a Random Forest ensemble. The objective here is to evaluate how each algorithmic family behaves when exposed to the target, and whether the underlying geometry of the dataset affects its ability to generalize.

## Project structure

The pipeline is organized into four sequential R Markdown notebooks (`.Rmd`) alongside their pre-rendered reports (`.html`), enabling direct examination of full analytical outputs without requiring a full execution. Each notebook represents a discrete phase of the project and pairs with a dedicated R utility script (`_utils.R` in the `R/` directory) containing custom auxiliary functions.

The `data/` directory manages data across multiple stages. Only `data/raw/` is committed to the repository (it holds the original datasets, kept for reproducibility). Three additional subfolders are generated during execution: `split/` (holding the quarantined train/test/folds partitions), `processed/` (holding preprocessed CSV files and recipe objects), and `models/` (holding trained, serialized models).

The repository also includes other key components, such as the isolated environment managed by `renv`, this documentation (`README.md`), and an end-to-end execution script (`run_all.R`).

``` text
clinical-ml-pipeline/
├── notebooks/
│   ├── 01_EDA_and_preprocessing.Rmd
│   ├── 01_EDA_and_preprocessing.html
│   ├── 02_unsupervised_modeling.Rmd
│   ├── 02_unsupervised_modeling.html
│   ├── 03_supervised_modeling.Rmd
│   ├── 03_supervised_modeling.html
│   ├── 04_clinical_synthesis.Rmd
│   ├── 04_clinical_synthesis.html
│   └── images/
├── R/
│   ├── 01_utils.R
│   ├── 02_utils.R
│   ├── 03_utils.R
│   └── 04_utils.R
├── data/
│   ├── raw/
│   ├── split/
│   ├── processed/
│   └── models/
├── README.md
├── run_all.R
├── renv/
├── renv.lock
└── clinical-ml-pipeline.Rproj
└── ...
```

***Note:** GitHub displays `.html` files as raw code instead of rendering them. To view the compiled notebooks properly, please download the files and open them locally in your browser.*

This project leans heavily on the `tidymodels` ecosystem, and two of its concepts show up repeatedly across the notebooks:

-   **Recipes (`recipes`):** define the preprocessing steps applied to the predictors (imputation, transformations, standardization, dummy encoding...), fit exclusively on training data to prevent leakage.
-   **Workflows (`workflows`):** bundle a given supervised model together with its recipe into a single object, so the entire preprocessing, training, and prediction process can be applied as one isolated unit, without ever manually re-running individual recipe steps outside of it.

Workflows are used specifically for the three supervised models, while the unsupervised route uses a classical implementation without `workflow` objects.

The data flow and primary tasks across the pipeline are structured as follows:

-   **Notebook 01 (`01_EDA_and_preprocessing.Rmd`):** ingests raw data, isolates the stratified train/test/folds partitions, audits statistical assumptions, and conducts exploratory data analysis. It builds three specialized `recipes` (`_cluster`, `_logistic`, `_trees`). The `_cluster` recipe is prepped and baked immediately to export processed CSVs (`data/processed/`), while the `_logistic` and `_trees` recipes are saved unbaked as serialized RDS objects to prevent data leakage.
-   **Notebook 02 (`02_unsupervised_modeling.Rmd`):** reads the processed cluster CSVs to conduct blind geometric exploration via K-Means, K-Medians, and DBSCAN/OPTICS. It determines optimal hyperparameters ($k$ and $\text{minPts}$/$ε$), evaluates internal clustering validity (Silhouette, Davies-Bouldin, Calinski-Harabasz), interprets the resulting topologies, and serializes the final models to `data/`
-   **Notebook 03 (`03_supervised_modeling.Rmd`):** imports the unbaked recipes, bundles them into `workflows` with three predictive architectures (Penalized Logistic Regression, C5.0 with a cost matrix, and Random Forest), and tunes hyperparameters via 10-fold cross-validation using a clinically prioritized metric hierarchy. It fits the final models, extracts coefficients and feature importances, performs internal validation on the training set, and serializes the workflows to `data/models/`.
-   **Notebook 04 (`04_clinical_synthesis.Rmd`):** does not train new models; it loads all serialized objects from Notebooks 02 and 03 to benchmark both unsupervised and supervised models against the quarantined 20% test set. It maps unsupervised clusters to clinical classes by majority vote, evaluates generalization metrics on blind data, and closes with the project's global synthesis, limitation analysis, and conclusions.

## Reproducibility

This project uses [`renv`](https://rstudio.github.io/renv/) to lock every package to the exact version used during development (tested with R 4.5; Windows users may need RTools to compile packages from source).

Opening `clinical-ml-pipeline.Rproj` in RStudio is strictly required: it triggers `.Rprofile` to bootstrap `renv` and activates the isolated project library without touching your global R installation.

The raw datasets are committed directly under `data/raw`, so the pipeline runs end to end without requiring live API access, and keeps working even if a future Kaggle update breaks the automated download. (Notebook 01 does include an automated Kaggle ingestion path, but it only triggers if a required file is missing locally, and expects valid `kaggle_user` and `kaggle_key` credentials to do so).

## Running the pipeline

For running the pipeline, you need R 4.5. and RStudio installed. Clone or download the repository, then open `clinical-ml-pipeline.Rproj` in RStudio to load the project environment.

***Note:** If you are on Windows, ensure you add the project's root folder to your Windows Defender exclusions before running the code. This prevents Windows Defender or Smart App Control from falsely blocking project files (such as `.dll` libraries) and halting execution.*

You can then run the pipeline using either of the following workflows:

#### Option A: automated end-to-end execution (recommended)

Run the root script `run_all.R`. It automatically calls `renv::restore()` to ensure environment consistency and renders notebooks `01` through `04` in sequential order within isolated sessions, halting immediately if any error occurs:

-   **From RStudio:** Open `run_all.R` and click Source.
-   **From Console:** Run `source("run_all.R")`.

#### Option B: step-by-step manual execution

If you prefer running or knitting notebooks individually from RStudio:

1.  Restore the environment manually once in the console:

``` r
renv::restore() 
```

2.  Knit the notebooks (`01_EDA_and_preprocessing.Rmd` → `04_clinical_synthesis.Rmd`) in strict numerical order. Each notebook depends directly on artifacts serialized into `data/split/`, `data/processed/`, or `data/models/` by preceding steps; skipping or reordering notebooks breaks the execution chain.

    ***Note:** Always knit to HTML. Attempting to knit to PDF or Word may fail due to missing system dependencies (such as LaTeX engines or Office converters) that are intentionally omitted from this minimal setup.*

In both workflows, compiled HTML reports will be generated directly within the `notebooks/` directory, allowing you to review the fully rendered analyses alongside their executed outputs.

## Methodological principles

A few design decisions run consistently across all four notebooks, and are worth stating up front:

-   **Zero data leakage:** the test set (20% of each dataset) is isolated in the very first block of Notebook 01 and remains untouched until the external evaluation in Notebook 04. Every statistic used for imputation, transformation, or standardization is calculated exclusively on the training partition.
-   **Recipe specialization over a single shared pipeline:** standardization and outlier winsorization are mathematically necessary for centroid-based and parametric models, but irrelevant, or even counterproductive, for tree-based ones. Building three separate recipes instead of forcing a single pipeline preserves both mathematical validity and clinical interpretability (e.g., reading a C5.0 split in real biological units instead of a dimensionless Z-score).
-   **Clinical cost asymmetry:** a false negative in an oncological or cardiovascular triage context is a fundamentally different error than a false positive. This shapes the entire evaluation strategy, from the metric hierarchy used during hyperparameter selection (Sensitivity \> PR-AUC \> F1-Score \> ROC-AUC \> Accuracy) to the explicit cost matrix built into the C5.0 models.
-   **Honest reporting over metric maximization:** every model is evaluated against the same hidden test set, and every result is reported as obtained, including the ones that reveal a model's failure to generalize. The cardiovascular scenario, in particular, is used throughout the project to illustrate the limits of standard algorithms against a genuinely overlapped, non-linear biological space.

## Summary of findings

For breast cancer, every modeling paradigm converges on the same conclusion: the underlying biological signal is largely linear and well separated, so a simple, interpretable model (Penalized Logistic Regression at an adjusted threshold) performs on par with more complex non-linear architectures. Even some unsupervised algorithms (K-Means and K-Medians) manage to geometrically recover the separation between malignant and benign cases, without ever seeing the target.

For cardiovascular risk, the opposite pattern holds at every stage of the pipeline, from the unsupervised exploration to the final external evaluation: the three risk tiers overlap heavily in feature space, and no algorithm evaluated here, including Random Forest, achieves a clinically deployable result across the full three-class problem. The high-risk class is consistently the best discriminated, while the low-risk class is consistently the hardest to capture, often disappearing entirely from a model's predictions.

The full methodological reasoning, the exact metrics behind these conclusions, and the clinical interpretation of each result are developed in detail across the four notebooks, closing with a dedicated synthesis section in `04_clinical_synthesis.Rmd`.
