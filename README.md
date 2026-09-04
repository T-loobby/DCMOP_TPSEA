# DCMOP_TPSEA

Minimal MATLAB package for running **TPSEA** on the **DCF1--DCF10** dynamic constrained multi-objective benchmark suite.

## Run

Open MATLAB in the repository root and execute:

```matlab
Main
```

This runs TPSEA once on every problem from DCF1 through DCF10 and saves results under `Results/TPSEA/`.

For a short installation check:

```matlab
Main(true)
```

## Requirements

- MATLAB
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox (used by the reference-front implementations of DCF4 and DCF8)

The population size is fixed at 100 because TPSEA's history model assumes `N = 100`.
