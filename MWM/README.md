# Morris Water Maze Analysis

Small demo data are included in `demo`. The demo data contain behavior tracking data used by the Morris water maze analysis scripts.

Before running the scripts, load or define the required variables for each analysis, such as `Track`, `AnimalID`, `numAnimal`, `Group1`, and `Group2`. Run each script from the folder indicated below.

## Turn Behavior

Folder: `01_Turn_bahavior`

Run with MATLAB:

```matlab
Run_Quantitate_shortTurn
```

### Expected outcome

Running `Run_Quantitate_shortTurn.m` detects short-turn events from animal tracking data and summarizes navigation behavior around the target platform during the probe trial.

The script detects trajectory self-crossing points, defines turn segments, excludes very long turns and very short turns, and visualizes the first detected turns for each animal.

Main outputs include:

- `OmitLongTurn crossTurn4 Probe Track animal*.jpg`
- `MinDistance_ToPlatform.xlsx`
- `Ratio_Target.xlsx`
- `CumsSumRatio_Target.xlsx`

`MinDistance_ToPlatform.xlsx` contains the minimum distance to the target platform for detected turns. `Ratio_Target.xlsx` contains target-quadrant occupancy ratios during each turn. `CumsSumRatio_Target.xlsx` contains cumulative target-quadrant occupancy ratios across turns. Each Excel file is separated into `Group1` and `Group2` sheets.

## Swimming Direction

Folder: `02_Swimming_direcition`

Run with MATLAB:

```matlab
Run_Probe_Track_angle
```

### Expected outcome

Running `Run_Probe_Track_angle.m` estimates the initial swimming direction during the first 5 seconds of the probe trial by fitting a line to the early trajectory of each animal.

The main output is:

- `Angel_first5.xlsx`

`Angel_first5.xlsx` contains the estimated initial swimming angle for each animal. The filename is kept as written in the script.

Together, these outputs quantify early search direction, turn structure, distance to the former platform location, and target-zone occupancy during probe-trial behavior.
