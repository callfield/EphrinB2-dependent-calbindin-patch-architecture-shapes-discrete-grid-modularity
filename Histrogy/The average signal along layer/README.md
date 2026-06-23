# Average Signal Along Layer

Small demo data are included in `demo`.

Run the script from this folder. Update `dataDir` and the CSV file name in `Run_NormSection.m` before running the analysis on new data.

## Normalized Layer Signal

Run with MATLAB:

```matlab
Run_NormSection
```

### Expected outcome

Running `Run_NormSection.m` calculates a normalized intensity profile along the selected layer axis. The input intensity profile is resampled into `nSec` equally spaced sections, normalized to a 0-1 range, and smoothed with a moving average.

The main output is the MATLAB workspace variable:

- `Norm_intensity`

`Norm_intensity` contains the normalized signal profile across the layer. This variable can be used for plotting or for downstream comparison of average signal distribution along the layer.
