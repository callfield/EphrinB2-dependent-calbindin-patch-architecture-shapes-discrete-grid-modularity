# Iba1 And GFAP Quantification

Small demo data are included in `demo`.

Run each script from this folder. Update the input CSV paths and group definitions before running the analysis on new data.

## Iba1 Signal Quantification

Run with MATLAB:

```matlab
Iba1_quantitate
```

### Expected outcome

Running `Iba1_quantitate.m` quantifies normalized Iba1 immunofluorescence intensity for negative-control and Casp/Cre-positive samples. For each slice, the target signal intensity is divided by the corresponding reference intensity, and the values are normalized to the mean of the negative-control group.

The main output file is:

- `Iba1_normalIntensity_Slice.emf`

The plot shows normalized Iba1 intensity for control and Casp/Cre-positive samples with group means, SEM, and individual slice-level data points. The MATLAB command window also reports the statistical comparison between groups, including the p value, t statistic, degrees of freedom, and group mean +/- SD.

## GFAP Signal Quantification

Run with MATLAB:

```matlab
GFAP_quantitate
```

### Expected outcome

Running `GFAP_quantitate.m` quantifies normalized GFAP immunofluorescence intensity for negative-control and Casp/Cre-positive samples. For each slice, the target signal intensity is divided by the corresponding reference intensity, and the values are normalized to the mean of the negative-control group.

The main output file is:

- `GFAP_normalIntensity_Slice.emf`

The plot shows normalized GFAP intensity for control and Casp/Cre-positive samples with group means, SEM, and individual slice-level data points. The MATLAB command window also reports the statistical comparison between groups, including the p value, t statistic, degrees of freedom, and group mean +/- SD.
