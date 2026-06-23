# Reelin MEC Calcium Imaging Analysis

Small demo data are included in `MEC/demo`. The demo data include IDPS-processed CSV files and behavior tracking data generated with `OF_analysis_PCI_1.m`.

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

Running `OF_analysis_PCI_3.m` detects calcium events, generates spatial firing-rate maps, computes autocorrelograms, calculates grid scores, estimates grid scale and orientation, and identifies grid cells. The main output is an updated `ST_dF_grid_aut_data.mat` containing event timing, moving-event timing, rate maps, occupancy maps, autocorrelograms, grid scores, shuffled grid-score thresholds, z-scored grid scores, grid-cell identities, grid scale, and grid orientation.

The script also generates per-cell activity plots, firing maps, autocorrelogram summaries, grid-cell summary images, and distribution plots for grid scale, grid orientation, and grid width.

## Summarize Grid Analysis

Run with MATLAB from the Reelin folder:

```matlab
a_MakeDirFile
b_Run_GridBregma
c_MakeTable
```

### Expected outcome

Running `a_MakeDirFile.m` creates or updates `Data.mat` with the directory paths for each animal and trial in group 1 and group 2. These paths are used by downstream group-level analysis scripts.

Running `b_Run_GridBregma.m` converts cell-position and grid-cell measurements into anatomical coordinates relative to bregma. The script creates the `fromBregma` folder and appends `GROUP1` and `GROUP2` variables to `Data.mat`, which contain cell-level anatomical position and grid-property summaries for each animal and trial.

Running `c_MakeTable.m` combines per-trial grid-analysis outputs with anatomical position, head-direction score, spatial information score, border score, and grid-module identity. The main outputs are:

- `Group1_CellTables.mat`
- `Group2_CellTables.mat`

The script also creates linked per-trial files such as `GSrate_map.mat`, `dF_data.mat`, and `Trk_wTS.mat` for rate maps, calcium traces, and tracking data with timestamps.

## Spatial Modulation Of Calcium Activity

Folder: `02_Spatial_cell`

Run with MATLAB:

```matlab
Run_conditional_entropy_Zscore
```

### Expected outcome

Running `Run_conditional_entropy_Zscore.m` computes spatial modulation of calcium activity for each cell. The script saves:

- `Spatial_Info.mat`
- `NormSpatial_Info.xlsx`

`Spatial_Info.mat` contains spatial information metrics for group 1 and group 2. `NormSpatial_Info.xlsx` contains normalized spatial information values for each group.

The script also creates the `conditional_entropy_zscore` folder with group-comparison CDF plots, including:

- `Average spatial info per spike cdf GROUP1 vs GROUP2.pdf`
- `Normalized spatial info cdf GROUP1 vs GROUP2.pdf`

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

Per-cell direction tuning figures are generated in each trial folder under `Direction/Head_direction_nofilter`. The output includes no-speed-filter and moving-only head-direction metrics.

## Grid Scale

Folder: `05_G_Scale`

Run with MATLAB:

```matlab
Run_Scatter_scale
```

### Expected outcome

Running `Run_Scatter_scale.m` summarizes the relationship between grid scale and anatomical position. The script generates group-wise scatter plots with fitted slopes and a region-adjusted comparison:

- `GScale_ScatterwithSlope gp1 all trial.pdf`
- `GScale_ScatterwithSlope gp2 all trial.pdf`
- `G_scale_regionAdjusted_all_trial.emf`

## Grid Field

Folder: `06_G_Field`

Run with MATLAB:

```matlab
Run_Scatter_field
```

### Expected outcome

Running `Run_Scatter_field.m` summarizes the relationship between grid field width and anatomical position. The script generates group-wise scatter plots with fitted slopes and a region-adjusted comparison:

- `G_field_scatter_with_slope_group1_all_trial.emf`
- `G_field_scatter_with_slope_group2_all_trial.emf`
- `G_field_regionAdjusted_all_trial.emf`

## Grid Orientation

Folder: `07_G_Orient`

Run with MATLAB:

```matlab
Run_Scatter_ori
```

### Expected outcome

Running `Run_Scatter_ori.m` summarizes grid orientation across anatomical position and group. The script generates:

- `G_ori_scatter_group1_all_trial.pdf`
- `G_ori_scatter_group2_all_trial.pdf`
- `G_ori_KSD_all_trial.pdf`

These figures show grid-orientation distributions and group-wise spatial trends.

## Identification Of Grid Modules

Folder: `08_GridModule`

Run with MATLAB:

```matlab
Run_a_Eachtrial_GridScaleRatio_2D
Run_b_Define_2DModule
Run_c_ksdNorm_scale_colorscatter
Run_d_Min_Gscale
```

### Expected outcome

Running these scripts identifies grid modules based on trial-wise grid-scale relationships. The main output is:

- `GridMod.mat`

`GridMod.mat` contains grid-scale ratios, module assignments, module-specific grid scales, and related group-wise summary variables.

Additional outputs include:

- `NormGS_Gmodtrial/NormGS.pdf`
- `NormGS_Gmodtrial/NormGS ScatterwithSlope group1.pdf`
- `NormGS_Gmodtrial/NormGS ScatterwithSlope group2.pdf`
- `0-150_mod1_GridS.xlsx`

These outputs summarize normalized grid scale, module-related spatial gradients, and minimum grid-scale values in the 0-150 um anatomical range.

## Speed Modulation Analysis Between Grid Modules

Folder: `09_Grid_module_speed`

Run with MATLAB:

```matlab
Run_a_Speed_vs_meandFCorr2_NewP
Run_b_Speed_zScore_Hz_each1
Run_c_Speed_zScore_Hz_each5
Run_d_zSpeedeach1_eachModule
Run_e_zSpeedeach1_eachModule_Linfit
```

### Expected outcome

Running these scripts quantifies how grid-cell activity and grid-module activity vary with running speed. The main outputs are:

- `SpeedCorr2_NewP.mat`
- `ZSpeedHz_each1.mat`
- `ZSpeedHz.mat`
- `ZScored_each1_gp1_all.xlsx`
- `ZScored_each1_gp2_all.xlsx`
- `ZScored_Speed_DIF.xlsx`
- `SpeedLine_fit.xlsx`

The scripts also generate speed-response figures, including z-scored activity versus dorsal-ventral speed bins and module-wise speed-modulation plots:

- `ZScored_all_trial_each1/ZScored_DVspeed_*.pdf`
- `ZScored_all_trial/ZScored_DVspeed_*.jpg`
- `GridMod_ZScored_DVspeed_gp1.pdf`
- `GridMod_ZScored_DVspeed_gp2.pdf`

## Spatial Specificity Of Population Activities Of Grid Cells

Folder: `10_SpatialUniqness`

Run with MATLAB:

```matlab
Run_a_SpatialUniquness_FullVec_Grid
Run_Fig_SpatialUniquness_Grid
```

### Expected outcome

Running `Run_a_SpatialUniquness_FullVec_Grid.m` computes spatial uniqueness metrics for population activity of grid cells and saves:

- `SpatialUniquness_Data.mat`

Running `Run_Fig_SpatialUniquness_Grid.m` summarizes the decoding-error distributions and exports Excel tables and figures, including:

- `GridMaxFR_MaxErroDist_*.xlsx`
- `MaxErrorDist_log_cdf_GROUP1vsGROUP2_*.pdf`
- `MaxErrorDist_log_hist_GROUP1vsGROUP2_*.jpg`
- `MaxErrorDist_nolog_cdf_GROUP1vsGROUP2_*.jpg`
- `MaxErrorDist_nolog_hist_GROUP1vsGROUP2_*.jpg`
- `MaxErrorDist_nolog_hist_GROUP1vsGROUP2_*.pdf`
- `MaxErrorDist_nolog_hist_GROUP1vsGROUP2_*.emf`

These outputs quantify how uniquely grid-cell population activity represents spatial position.

## Animal Position Decoding From Non-Grid Spatial Cells

Folder: `11_Decoding`

Run with MATLAB and Python:

```matlab
a_RunBatchExportVariable
```

```bash
python "b_decoding batch bestSpatialInfo v5.1 Spatial_nogrid.py"
python "c_decoding dist Summary.py"
python "d_decoding dist Summary 2.py"
```

### Expected outcome

Running `a_RunBatchExportVariable.m` exports the variables required for Python-based decoding analysis.

Running the Python scripts trains and evaluates LSTM-based animal-position decoding from non-grid spatial cells. The decoding script generates per-session decoding results and diagnostic figures such as:

- `*_lstm_results.npz`
- `*_shuffle_convergence.jpg`
- `*_error_distribution.jpg`

The summary scripts generate Excel files containing decoding-error distribution statistics:

- `lstm_noGrid_bestSPi50cell_errorDist_with_stats.xlsx`
- `lstm_noGrid_bestSPi50cell_errorDist_with_stats2.xlsx`

## Grid Phase Difference

Folder: `12_G_PhaseDiff`

Run with MATLAB:

```matlab
Run_GridPhaseDiff
```

### Expected outcome

Running `Run_GridPhaseDiff.m` computes pairwise grid-phase differences and their relationship to physical distance between cells. The main outputs are:

- `Grid_phase.mat`
- `GridPhase_GROUP1.mat`
- `GridPhase_GROUP2.mat`

The script also generates phase-distance plots and all-pair grid-phase summary figures in the `Phase_vs_distance` and `AllP_Phase_Diff` folders.

## Grid Scale-Field Width Relationship And Module Distance

Folder: `13_GridFigures`

Run with MATLAB:

```matlab
a_GridScore
b_FiringRate
c_GridScale
d_DVpositionVsRawGridScale
e_GridScaleVsWidth
f_Fig_ScalePerWidth
g_SpeedModRate
h_GridPhaseDistribution
i_GridOri
```

### Expected outcome

These scripts generate publication-ready summary figures for grid-cell properties and module organization. The main outputs include:

- `GridScore.emf`
- `AnimalWise_Frate *.emf`
- `GridScale_2.emf`
- `DVposVsScale.emf`
- `EB2_Scatter_GridScaleVsWidth_All_260511.emf`
- `FigS GridScalePerWidth_Group2.emf`
- `SpeedModRatio_Group2.emf`
- `FigS GridPhaseDist_Group2.emf`
- `FigS_GridOriDiff_Group2.emf`

Together, these figures summarize grid score, firing rate, grid scale, dorsal-ventral position, grid field width, speed modulation, grid-phase distribution, and grid-orientation differences.
