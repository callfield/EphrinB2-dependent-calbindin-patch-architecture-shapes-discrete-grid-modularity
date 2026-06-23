# Quantitative Analysis Of CalB Patches

Small demo data are included in `demo`.

Run the script from this folder. Update the input sample list, ROI files, and CSV paths before running the analysis on new data.

## CalB Patch Quantification

Run with MATLAB:

```matlab
Run_clusteranalysis
```

### Expected outcome

Running `Run_clusteranalysis.m` detects CalB patch peaks from ROI-based intensity profiles. The script first performs intensity and background normalization, then detects peaks along the analyzed axis and along the superficial-deep axis.

The main output files are written to the `results` folder:

- `results/pks.xlsx`
- `results/locs.xlsx`
- `results/width.xlsx`
- `results/peaks_from_varey.xlsx`
- `results/peak_center.xlsx`
- `results/sng_SD-axis_pks.xlsx`
- `results/sng_SD-axis_locs.xlsx`
- `results/sng_SD-axis_width.xlsx`
- `results/sng_SD-axis_peaks_from_varey.xlsx`

These files contain detected peak intensities, peak locations, estimated peak widths, peak prominence values, peak center positions, and superficial-deep axis peak measurements.

The `results` folder also contains `.jpg` quality-control images showing:

- Raw intensity profiles.
- Background-normalized signal profiles.
- Gaussian-filtered profiles.
- Detected peaks.
- ROI overlays on the original tissue images.

These images can be used to confirm that the detected CalB patches match the input ROI and fluorescence signal.
