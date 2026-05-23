# Tail-Risk Resilience TVP-VAR Analysis

This R code implements a TVP-VAR-based tail-risk resilience model for Chinese
commodity futures. It estimates time-varying responses of downside tail risk to
pessimistic tone shocks, ranks commodities by resilience, and evaluates
downside-risk warning performance.

## Contents

- `data/analysis_dataset.csv`: commodity-date analysis panel
- `scripts/run_empirical_analysis.R`: TVP-VAR analysis script
- `output/`: generated results

The analysis dataset contains the following input variables:

- `date`: trading date
- `commodity`: commodity name
- `symbol`: futures symbol
- `sector`: commodity sector
- `close`: closing price
- `volume`: trading volume
- `open_interest`: open interest
- `return`: daily log return
- `tr`: downside tail-risk measure
- `tone`: pessimistic tone measure

## Usage

Run the empirical analysis from the project root:

```bash
Rscript scripts/run_empirical_analysis.R
```

On Windows, if R is not on your system path, use the full path to `Rscript.exe`:

```powershell
D:\Soft\R-4.5.2\bin\x64\Rscript.exe scripts\run_empirical_analysis.R
```

The script reads `data/analysis_dataset.csv`, estimates the TVP-VAR model, and
writes the generated figures and tables to `output/`.
