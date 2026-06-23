# CalB MEC Calcium Imaging Analysis

Small demo data may be available in `MEC/demo`. The demo data include IDPS-processed CSV files and behavior tracking data generated with `OF_analysis_PCI_1.m`.

Run each script from the folder indicated in the corresponding section. Update placeholder paths in the scripts before running them on new data.

## Animal Tracking

Folder: `01_Each_animal`

Run with MATLAB:

```matlab
OF_analysis_PCI_1
```

### Expected outcome

Running `OF_analysis_PCI_1.m` extracts animal tracking data from the behavior video and detects calcium-recording periods from the lamp video. The script saves `RawTrack.mat`, which contains the raw tracking data and recording-period information.

The script also creates a `Track` folder with quality-control images:

- `Track/all rec track.jpg`
- `Track/all moving rec track.jpg`
- `Track/all rec velocity.jpg`

These images can be used to confirm the extracted trajectory, movement-filtered trajectory, and velocity trace.

## Cell Detection And Rate Map

Calcium imaging data should be processed with Inscopix Data Processing Software before running these scripts.

Folder: `01_Each_animal`

Run with MATLAB:

```matlab
OF_analysis_PCI_2
OF_analysis_PCI_3
```

### Expected outcome

Running `OF_analysis_PCI_2.m` aligns tracking data with the calcium imaging data exported from Inscopix Data Processing Software. It generates processed calcium and behavior files:

- `ST_PCI_Ca_behav_track.csv`
- `ST_PCI_dF.csv`
- `ST_PCI_noDup_dF.csv`
- `ST_PCI_noDup_CellPos.csv`
- `Original_Cell_ID_noDup.csv`
- `ST_dF_grid_aut_data.mat`

The script also generates quality-control images such as:

- `Track/all track.jpg`
- `Track/all moving track.jpg`
- `Removed Cells.jpg`

Running `OF_analysis_PCI_3.m` detects calcium events, generates spatial firing-rate maps, computes autocorrelograms, calculates grid scores, estimates grid scale and orientation, and identifies grid cells. The main output is an updated `ST_dF_grid_aut_data.mat` containing event timing, moving-event timing, rate maps, occupancy maps, autocorrelograms, grid scores, shuffled grid-score thresholds, z-scored grid scores, grid-cell identities, grid scale, grid orientation, and grid width.

The script also generates per-cell activity plots, firing maps, autocorrelogram summaries, grid-cell summary images, and distribution plots for grid scale, grid orientation, and grid width.

## Summarize Grid Analysis

Run with MATLAB from the CalB folder:

```matlab
MakeDirFile
Run_GridBregma
```

### Expected outcome

Running `MakeDirFile.m` creates or updates `Data.mat` with the directory paths for each animal and trial in group 1 and group 2. These paths are used by downstream group-level analysis scripts.

Running `Run_GridBregma.m` converts cell-position and grid-cell measurements into anatomical coordinates relative to bregma. The script creates the `fromBregma` folder and appends `GROUP1` and `GROUP2` variables to `Data.mat`, which contain cell-level anatomical position and grid-property summaries for each animal and trial.

## Spatial Modulation Of Calcium Activity

Folder: `02_Spatial_cell`

Run with MATLAB:

```matlab
Run_conditional_entropy_Zscore
```

Optional visualization:

```matlab
Run_Visualise_spatial_info
```

### Expected outcome

Running `Run_conditional_entropy_Zscore.m` computes spatial modulation of calcium activity for each cell. The script saves:

- `Spatial_Info.mat`
- `NormSpatial_Info.xlsx`

`Spatial_Info.mat` contains spatial information metrics for group 1 and group 2. `NormSpatial_Info.xlsx` contains normalized spatial information values for each group.

The script also creates the `conditional_entropy_zscore` folder with group-comparison CDF plots, including:

- `Average spatial info per spike cdf GROUP1 vs GROUP2.pdf`
- `Normalized spatial info cdf GROUP1 vs GROUP2.pdf`

Running `Run_Visualise_spatial_info.m` creates per-cell spatial-information visualizations in `allcell` and saves `All_Info.mat` for group-level visualization outputs.

## Border Selectivity

Folder: `03_Border`

Run with MATLAB:

```matlab
Run_BorderScore_v2
```

### Expected outcome

Running `Run_BorderScore_v2.m` computes border scores for each cell and summarizes the results by group. The main outputs are:

- `Data_Borderv2.mat`
- `BorderScore_v2.xlsx`
- `Border Score v2.jpg`

`Data_Borderv2.mat` stores the group-wise border-score variables. `BorderScore_v2.xlsx` contains border-score values for group 1 and group 2. The figure output shows group-level border-score distributions and summary comparisons.

Per-cell visualization files are generated in folders such as `RateMap`, `RateMap_vsV1`, and `Visualise`.

## Head-Direction Tuning

Folder: `04_Head_Direction`

Run with MATLAB:

```matlab
Run_a_direction_Analysis_nofiler_Raylign
Run_b_nofilterHD_pilot_RL
```

### Expected outcome

Running `Run_a_direction_Analysis_nofiler_Raylign.m` extracts head and body angle information from DLC-tracked body-part coordinates and saves `Angle.mat` for each trial.

Running `Run_b_nofilterHD_pilot_RL.m` computes head-direction tuning scores and saves:

- `Data_HD.mat`
- `HDScore.xlsx`

Per-cell direction tuning figures are generated in each trial folder under `Direction/Head_direction_nofilter` and in the group-level `nofilter_HDcell` folder. The output includes no-speed-filter and moving-only head-direction metrics.

## Synchronization Analysis

Folder: `05_Synchronization`

Before running these scripts, manually define CalB-patch cell groups in `Island.mat`.

Run with MATLAB:

```matlab
Run_a_Ex_DV_meandF_Island
Run_b_SubCorrelation_Island
Run_c_Allcellpair_summary
```

### Expected outcome

Running `Run_a_Ex_DV_meandF_Island.m` generates example mean dF/F traces for intra- and trans-CalB+ patch cells. The main output folder is:

- `dFex_230609`

Running `Run_b_SubCorrelation_Island.m` computes subpopulation correlation metrics for intra- and trans-patch cell pairs. The script creates:

- `Ratio`
- `MeanCor`
- `Data.mat`

Running `Run_c_Allcellpair_summary.m` summarizes all cell-pair correlation results across groups. The main outputs include:

- `Corr_cellpair.xlsx`
- `Corr_cellpair_allEB2.csv`
- `ViolinPlot_ALLmat_vs_ALLeb_mat.pdf`
- `Sum`
