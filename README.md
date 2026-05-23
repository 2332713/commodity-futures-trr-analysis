# R Empirical Analysis

This folder contains a self-contained R workflow for reproducing the empirical
figures and tables from a prepared analysis dataset.

## Contents

- `data/analysis_dataset.csv`: analysis dataset
- `scripts/run_empirical_analysis.R`: R script for generating figures and tables
- `output/`: generated results

## Usage

Run the empirical analysis from the project root:

```bash
Rscript scripts/run_empirical_analysis.R
```

On Windows, if R is not on your system path, use the full path to `Rscript.exe`:

```powershell
D:\Soft\R-4.5.2\bin\x64\Rscript.exe scripts\run_empirical_analysis.R
```

The script reads `data/analysis_dataset.csv` and writes the generated
figures and tables to `output/`.
