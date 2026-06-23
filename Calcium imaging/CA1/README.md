# CA1 Calcium Imaging Analysis

This folder contains the CA1 calcium-imaging analysis pipeline for animal tracking, cell validation, place-cell analysis, and Bayesian decoding.

Before running the scripts, update the hard-coded data and output paths in each script as needed. Run each script from the folder indicated below.

## Animal Tracking

Folder: `1. AnimalTrack`

Run the ROI-definition scripts from `1. AnimalTrack/DefineROI`:

```matlab
DefineArenaROI_260202
DefineLightROI_260202
LightDuration_260202
```

Then run the tracking scripts from `1. AnimalTrack`:

```matlab
Run_VideoTrack_260202
Run_TrackRefine_260205
```

### Expected outcome

The arena ROI, LED ROI, LED-on frames, raw mouse trajectory, and refined trajectory are generated for each behavioral video.

Main outputs include:

- Arena ROI files in each `Track` folder, such as `cropRect_*.mat` and `Arena_GUI_*.jpg`.
- LED ROI and LED-duration files in each `Light` folder, such as `Light_cropRect_*.mat`, `Light_ROI_*.jpg`, `RecFrames_2_*.mat`, and cropped LED videos.
- Raw tracking files in each `Track` folder, such as `Track_*.mat`, background images, and tracking summary images.
- Refined trajectory files in `TrackResults_260204`, including `Track *.jpg`, `Velocity *.jpg`, `Velocity_Epoch *.jpg`, and `TrackPos_*.mat`.

## Cell Validation

Folder: `2. CellValidation`

Run:

```matlab
CellValidation_seq_260204
```

### Expected outcome

Cells are validated using spatial overlap, correlation, skewness, and trace-quality criteria. Accepted cells and their aligned calcium traces are saved for downstream place-cell analysis.

Main outputs include:

- Validated cell files in `results_Validation`, such as `CellTraces_Validated_*.mat`.
- Accepted-cell variables including `PeakIdx_accepted`, `AcceptedMetric`, `full_ts`, `interp_traces_accepted`, `image_data`, `CellImage_acceptedCells`, and `Boundaries_cleaned`.
- Validation figures exported with prefixes such as `CellTraces_Validated_*_fig_originalImage`, `*_fig_originalDetection`, `*_fig_Detection_skew`, and `*_fig_Detection_skew_corr`.

## Place Cell Analysis

Folder: `3. PlaceCellAnalysis`

Run:

```matlab
ReDefineArenaROI_260202
Analysis_260204_seq_6MAD_AdaptiveBin
CleaningSavedData_260206
StatAnalysis_PlaceCell_260524
```

### Expected outcome

Arena boundaries are refined, calcium events are detected from validated traces, position-aligned rate maps are generated, place-cell metrics are calculated, and summary figures are created.

Main outputs include:

- Refined arena ROI files, such as `Arena_Redefined_*.jpg` and `Arena_Redefined_*.mat`.
- Spatial-analysis files, such as `SpaceAnalysis_*.mat`, containing rate maps, spike-detection results, spatial information, coherence, field metrics, and reliability metrics.
- Field-analysis and bootstrap/shuffle outputs used to evaluate place-cell reliability and spatial tuning.
- Compact cleaned analysis files, such as `CleanedAnalysisData_*.mat`, for downstream statistics.
- Summary data and figures from `StatAnalysis_PlaceCell_260524`, including `PlaceCellInfo_260218.mat`, place-cell percentage plots, cell-count plots, running-behavior plots, field metrics, stability metrics, reliability metrics, and speed-modulation figures.

## Bayesian Decoding

Folder: `4. BayesianDecoding`

Run:

```matlab
BayesianDecoding_260317
Fig_Bayesian_260525
```

### Expected outcome

Bayesian position decoding is performed using CA1 calcium activity, and decoding performance is summarized across cell populations and spatial regions.

Main outputs include:

- Decoding result folders grouped by cell count, such as `BayesResults_260317_BorderCenter/*cells`.
- Per-population decoding figures, such as `BayesDecoding_*.png` and `BayesDecoding_*.fig`.
- Decoding result files, such as `bayesian_decoding_results.mat` and `bayesian_decoding_summary.csv`.
- Scatter plots and final summary figures from `Fig_Bayesian_260525`, including decoding-error, z-scored decoding-error, and accuracy PDF files.
