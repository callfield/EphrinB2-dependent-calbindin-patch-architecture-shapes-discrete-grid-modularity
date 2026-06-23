Small damo data (bahavior track date)included ＠MWM\demo

Quantiate Turn
@ Turn_bahavior
1, Run“Run_Quantitate_shortTurn.m” with matlab

Swimming direction
@ Swimming_direcition
1, Run“Run_Probe_Track_angle.m” with matlab

Expected outcome

This folder contains two Morris water maze analyses: turn behavior during the probe trial and initial swimming direction.

Running `Run_Quantitate_shortTurn.m` detects short-turn events from animal tracking data and summarizes navigation behavior around the target platform. The script generates trajectory images for each animal showing the probe-track path and the first detected turns after excluding long turns.

Main output files:
- `OmitLongTurn crossTurn4 Probe Track animal*.jpg`: probe-track images with detected turn segments.
- `MinDistance_ToPlatform.xlsx`: minimum distance to the target platform for each detected turn, separated by group.
- `Ratio_Target.xlsx`: target-quadrant occupancy ratio during each turn, separated by group.
- `CumsSumRatio_Target.xlsx`: cumulative target-quadrant occupancy ratio across turns, separated by group.

Running `Run_Probe_Track_angle.m` estimates the initial swimming direction during the first 5 seconds of the probe trial by fitting a line to the early trajectory of each animal. The main output is:
- `Angel_first5.xlsx`: estimated initial swimming angle for each animal.

Together, these outputs quantify early search direction, turn structure, distance to the former platform location, and target-zone occupancy during probe-trial behavior.