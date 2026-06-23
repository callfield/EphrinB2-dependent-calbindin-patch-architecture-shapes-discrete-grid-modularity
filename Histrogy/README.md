# Histology Analysis

This folder contains histology analysis scripts for CalB patch quantification, layer-normalized signal profiles, and Iba1/GFAP immunofluorescence quantification.

Run each analysis from the folder indicated below. Update the input paths and file names in each script before running them on new data.

## Quantitative Analysis Of CalB Patches

Folder: `Quantitative analysis of CalB patches`

Run with MATLAB:

```matlab
Run_clusteranalysis
```

### Expected outcome

The script detects CalB patch peaks from ROI-based intensity profiles and exports peak intensity, location, width, prominence, and center-position measurements to the `results` folder. It also generates quality-control images for raw intensity profiles, background-normalized signal profiles, Gaussian-filtered profiles, detected peaks, and ROI overlays.

See `Quantitative analysis of CalB patches/README.md` for details.

## Average Signal Along Layer

Folder: `The average signal along layer`

Run with MATLAB:

```matlab
Run_NormSection
```

### Expected outcome

The script calculates a normalized intensity profile along the selected layer axis. The main output is the MATLAB workspace variable `Norm_intensity`, which can be used for plotting or downstream comparison of signal distribution along the layer.

See `The average signal along layer/README.md` for details.

## Iba1 And GFAP

Folder: `Iba1 and GFAP`

Run with MATLAB:

```matlab
Iba1_quantitate
GFAP_quantitate
```

### Expected outcome

The scripts quantify normalized Iba1 and GFAP immunofluorescence intensities for negative-control and Casp/Cre-positive samples. They generate slice-wise summary plots and report group-level statistical comparisons in the MATLAB command window.

See `Iba1 and GFAP/README.md` for details.
