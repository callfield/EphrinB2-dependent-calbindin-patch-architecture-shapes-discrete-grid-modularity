Small damo data included ＠Iba1 and GFAPg\demo


For quantification IbaI signal

1, Run“Iba1_quantitate.m” with matlab



For quantification GFAP signal

1, Run“GFAP_quantitate.m” with matlab



Expected outcome

Running `Iba1_quantitate.m` or `GFAP_quantitate.m` quantifies the normalized immunofluorescence intensity for negative-control and Casp/Cre-positive samples. For each slice, the target signal intensity is divided by the corresponding reference intensity, and the values are normalized to the mean of the negative-control group.

The scripts generate slice-wise summary plots:
- `Iba1_normalIntensity_Slice.emf`
- `GFAP_normalIntensity_Slice.emf`

Each plot shows the normalized intensity for control and Casp/Cre-positive samples with group means, SEM, and individual slice-level data points. The MATLAB command window also reports the statistical comparison between groups, including the p value, t statistic, degrees of freedom, and group mean +/- SD.