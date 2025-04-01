# 2025-swpo-mls-meeting

Model development for 2025 Southwest Pacific Ocean striped marlin stock assessment conducted during the joint Pacific Community (SPC) and NOAA Fisheries stock assessment modeling meeting held in Honolulu, Hawai'i from 20-24 January 2025.

Any model runs contained within this repository represent preliminary model explorations and should not be used '*as is*' for the consideration of management advice. Only models that have been presented and accepted to the Western and Central Pacific Fisheries Commission (WCPFC) Scientific Committee should be used as the basis for management advice.

## Getting Started

This project uses `renv` for package management to ensure reproducibility across environments. The repository is also configured to work with GitHub Codespaces for cloud-based development.

### Running the Models Locally

First, [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) this repository to your local machine.

#### Using Base R

1. Open an R terminal (recommended [version 4.4.0](https://cloud.r-project.org/) with RTools 4.4 installed)
2. Set your working directory to the cloned repository:
   ```r
   setwd("path/to/repository/")
   ```
3. Source the `.Rprofile` file to bootstrap `renv`:
   ```r
   source(".Rprofile")
   ```
4. If `renv` doesn't bootstrap automatically, run:
   ```r
   renv::restore()
   ```
5. Follow the prompts to install all required packages (including `data.table`, `magrittr`, `ggplot2`, `r4ss`, `FLR4MFCL`, and other dependencies)

#### Using RStudio

1. Open the repository as an [RStudio Project](https://bookdown.org/ndphillips/YaRrr/projects-in-rstudio.html)
2. RStudio should automatically detect the `renv` configuration and prompt for package installation
3. If `renv` doesn't bootstrap automatically, run:
   ```r
   renv::restore()
   ```

#### Using Visual Studio Code

1. Configure VS Code to work with R using the [R extension](https://github.com/REditorSupport/vscode-R)
2. Follow the specific [configuration steps for renv projects](https://github.com/REditorSupport/vscode-R/wiki/Working-with-renv-enabled-projects)
3. Open the repository folder in VS Code
4. Open an R terminal, which should prompt `renv` to bootstrap
5. If `renv` doesn't bootstrap automatically, run:
   ```r
   renv::restore()
   ```

## Running the Models with GitHub Codespaces

GitHub Codespaces provides a cloud-based development Linux environment with all dependencies pre-configured:

1. [Open a Codespace](https://docs.github.com/en/codespaces/developing-in-a-codespace/creating-a-codespace-for-a-repository#creating-a-codespace-for-a-repository) for this repository using the default options
2. Initial Codespace creation may take 15-20 minutes to complete
3. Once ready, open an R terminal in the Codespace
4. The `renv` package should bootstrap automatically
5. If `renv` doesn't bootstrap automatically, run:
   ```r
   renv::restore()
   ```
::: {.callout-caution}
Please note that the Linux version (*3.30.23.1*) of the most recent Stock Synthesis executable distributed with the repository `executables\stock-synthesis\3.30.23.1\ss3_linux` does not work with the flavor of Linux used in GitHub Codespaces so the preceding version (*3.30.22.1*) should be used `executables\stock-synthesis\3.30.22.1\ss3_linux`.
:::  ```

## Project Structure

The project is organized with the following directory structure:

```
project/
├── assets/
│   └── static/         # Output plots and figures
├── code/               # Main workflow scripts
│   └── helper-fns/     # Helper functions sourced by main scripts
├── data/               # Input data files
├── executables/
│   └── stock-synthesis/ # Stock Synthesis executables
├── html-dashboard/     # HTML viewer output
├── models/
│   ├── mfcl/           # MFCL model files
│   └── stock-synthesis/ # Stock Synthesis model files
└── shiny/              # Shiny app for interactive model exploration
```

- **`assets/static/`**: Contains all generated plots for reports and presentations
- **`code/`**: Houses the main workflow scripts for model development and visualization
- **`code/helper-fns/`**: Contains utility functions that support the main workflow
- **`data/`**: Stores input data files including catch histories, size compositions, and age data
- **`executables/`**: Contains different versions of Stock Synthesis executable files
- **`html-dashboard/`**: Destination for generated HTML viewer files
- **`models/`**: Stores input/output files for both MFCL and Stock Synthesis models
- **`shiny/`**: Contains files for the interactive Shiny application

## Model Workflow and Scripts

The main workflow scripts are located in the `code/` directory. These scripts progress from initial model conversion through refinement stages to the final model implementation.

### Model Development Scripts

- **`01-make-baseline-stock-synthesis-model.r`**: Creates the baseline Stock Synthesis model from MFCL model version "1979_20p3".
  - Converts from MFCL to annual model with length-based selectivity
  - Maintains single sex configuration and fixes initial F
  - Preserves seasonal index for Q4

- **`02-mls-ss3-1979.r`**: Updates the initial MLS model with:
  - Resolves NZ recreational selectivity by sharing terminal bins
  - Adds tail compression
  - Adds extra observation standard error for CPUE index
  - Estimates end logit parameter for selectivity
  - Downweights size composition data

- **`03-chg-selex-1979.r`**: Adjusts selectivity parameters for specific fleets (LL 2, AU REC 9, Index 15)

- **`04-start-1952.r`**: Extends the model period back to 1952
  - Adds historical catch data
  - Includes NZ recreational weight composition data
  - Changes parameterization of extra index observation error

- **`06-exclude-more-comp.r`**: Excludes additional composition data:
  - Drops index (15) length composition
  - Drops fishery 2 length composition
  - Drops fishery 2 weight composition for 1979

- **`07-catch-uncertainty.r`**: Incorporates catch uncertainty for pre-1979 Fleet 2

- **`12-CAAL-old-growth-SD.r`**: Implements conditional age-at-length (CAAL) model
  - Incorporates age data
  - Adds ageing error matrix
  - Enables estimation of growth parameters L1, L2, K

### Visualization and Analysis

- **`00-html-viewer.r`**: Generates HTML files for viewing individual model results, creating interactive dashboards that display model fits, parameter estimates, and diagnostic plots.

- **`00-pres-plots.r`**: Creates plots for presentation, focusing on model comparison time series. Generates sequential comparisons of depletion estimates across different model configurations to illustrate how each modeling decision affects population trajectory estimates.

- **`00-report-plots.r`**: Produces comprehensive plots for report documents, including:
  - Time series comparisons (SSB, depletion, recruitment)
  - Selectivity curve visualization by fleet
  - Index fit diagnostics with variance adjustments
  - Length and weight composition fits across time periods and fleets
  - Biological parameter visualization (growth, maturity)
  - Catch history and model fits to observed catch

- **`00-shiny-gather.r`**: Collects summary files from all models for the Shiny application, standardizing outputs for cross-model comparison.

- **`00-shiny-viewer.r`**: Script to launch the Shiny app for interactive model exploration, allowing users to dynamically filter and compare results across all model configurations.

## Helper Functions

The project includes numerous utility functions in the `code/helper-fns/` directory that support the main workflow:

### Model Management

- **`clean_dir.r`**: Cleans unnecessary files from model directories to save disk space, keeping only essential files needed for analysis.
- **`is_linux_os.r`**: Detects the operating system to ensure compatibility with the appropriate Stock Synthesis executable.
- **`make_ss_output.r`**: Generates Stock Synthesis output files in a controlled environment, facilitating consistent results across different systems.
- **`ss_transfer.r`**: Transfers model files between directories with appropriate executable selection based on the operating system.

### Data Processing

- **`parse-frq.r`**: Parses MFCL frequency (`.frq`) files into a simplified S4 class object, extracting catch, effort, penalty, length frequency, and weight frequency data.
- **`make-selex-par.r`**: Contains helper functions to define different selectivity parameter configurations, including:
  - `make_size_selex_par_24`: For double-normal selectivity
  - `make_size_selex_par_5`: For mirrored selectivity
  - `make_size_selex_par_1`: For logistic selectivity
  - `make_age_selex_par_11`: For non-parametric age-based selectivity

### Reporting and Visualization

- **`ss_model_html_viewer.r`**: Generates an HTML dashboard for interactive model exploration, creating a comprehensive visual representation of model outputs.
- **`summarize_ss_model.r`**: Extracts key model parameters, diagnostics, and estimates into standardized CSV files for easier analysis and comparison.
- **`summarize_fleets.r`**: Summarizes fleet structure, selectivity, and data components for each fishery in the model.
- **`report-helpers.r`**: Contains specialized functions for data extraction and visualization:
  - `extract_mfcl_quantities`: Extracts population metrics from MFCL models
  - `get_population_numbers`: Extracts numbers-at-age from Stock Synthesis models
  - `plot_model_comparison_ts`: Creates time series comparison plots across models
  - `plot_model_comparison_selex`: Visualizes selectivity patterns across models
  - `plot_index_comparison`: Compares fits to abundance indices
  - `plot_composition_comparison`: Displays fits to length and weight composition data
  - `plot_catch_comparison`: Shows catch fits across models
  - `plot_model_comparison_bio`: Compares biological parameters like growth and maturity

## Usage Examples

Here are some common usage examples for working with the repository:

### Generating HTML Dashboard for a Model

```r
# Create HTML dashboard for model "12-CAAL-old-growth-SD"
model_dir <- file.path("models", "stock-synthesis", "12-CAAL-old-growth-SD")
source("code/00-html-viewer.r")
ss_model_html_viewer(model_dir)
```

### Comparing Models with Time Series Plots

```r
# Compare depletion time series across multiple models
source("code/00-report-plots.r")
plot_model_comparison_ts(
  model_ids = c("07-catch-uncertainty", "12-CAAL-old-growth-SD"),
  model_stem = file.path("models", "stock-synthesis"),
  categories = c("Depletion (D)"),
  show_se = TRUE,
  save_plot = TRUE,
  save_dir = file.path("assets", "static"),
  plot_name = "depletion_comparison"
)
```

### Using the Interactive Shiny App

```r
# First gather model summary data
# Note that this step is only needed if you have added new models to the repository
source("code/00-shiny-gather.r")

# Then launch the Shiny app
source("code/00-shiny-viewer.r")
```

### Creating a New Model from an Existing One

```r
# Transfer base model to a new directory and make modifications
from_dir <- file.path("models", "stock-synthesis", "12-CAAL-old-growth-SD")
to_dir <- file.path("models", "stock-synthesis", "my-new-model")
source("code/helper-fns/ss_transfer.r")
ss_transfer(from_dir, to_dir, file.path("executables", "stock-synthesis"))

# Make modifications to data or control files
# ...
# Write out the files ...

# Run the model
run(dir = to_dir, exe = "ss3.exe")

# Summarize the results
source("code/helper-fns/summarize_ss_model.r")
summarize_ss_model(to_dir)
```


## License

The code contained in this repository is licensed under the GNU GENERAL PUBLIC LICENSE version 3 ([GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)).

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
