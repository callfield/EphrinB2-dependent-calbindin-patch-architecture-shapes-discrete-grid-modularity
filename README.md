# README

This repository contains analysis code and small demo datasets associated with the manuscript:

**EphrinB2-dependent calbindin patch architecture shapes discrete grid modularity in the medial entorhinal cortex**

This manuscript is currently under peer review. A final publication citation and DOI will be added after publication.

## Authors

Naoki Yamamoto<sup>1,4</sup>, Hisayuki Osanai<sup>1,4</sup>, Kritika Ramesh<sup>1</sup>, Indrajith R. Nair<sup>1</sup>, Mark Henkemeyer<sup>2,3</sup>, Sachie K. Ogawa<sup>1,3</sup> & Takashi Kitamura<sup>1,2,3,*</sup>

<sup>1</sup>Department of Psychiatry, University of Texas Southwestern Medical Center, Dallas, TX 75390, USA  
<sup>2</sup>Department of Neuroscience, University of Texas Southwestern Medical Center, Dallas, TX 75390, USA  
<sup>3</sup>Peter O'Donnell Brain Institute, University of Texas Southwestern Medical Center, Dallas, TX 75390, USA  
<sup>4</sup>Equally contributed  
<sup>*</sup>Corresponding author. Email: Takashi.Kitamura@UTSouthwestern.edu

## Overview

The code in this repository was used to analyze histology, Morris water maze behavior, and calcium-imaging data for the manuscript. The folders are organized by experimental modality. Each major analysis folder contains its own `README.md` with the expected inputs, run order, and expected outputs.

Small demo datasets are included in selected `demo` folders. Large raw datasets are not included in this repository.

## Repository Structure

```text
EphrinB2_GridModule
|-- Histrogy
|   |-- Quantitative analysis of CalB patches
|   |-- The average signal along layer
|   `-- Iba1 and GFAP
|-- MWM
|   |-- 01_Turn_bahavior
|   `-- 02_Swimming_direcition
`-- Calcium imaging
    |-- MEC
    |   |-- Reelin
    |   `-- CalB
    `-- CA1
```

## Main Analyses

### Histology

Folder: [`Histrogy`](Histrogy/README.md)

This folder contains MATLAB scripts for:

- Quantitative analysis of CalB patch structure.
- Normalized signal profiles along cortical layers.
- Iba1 and GFAP immunofluorescence quantification.

Main outputs include peak-measurement spreadsheets, normalized intensity profiles, summary figures, and group-comparison statistics.

### Morris Water Maze

Folder: [`MWM`](MWM/README.md)

This folder contains MATLAB scripts for probe-trial behavior analysis:

- Short-turn detection and turn-by-turn target-zone occupancy.
- Initial swimming direction estimation during the first 5 seconds of the probe trial.

Main outputs include trajectory images and Excel summary files for turn behavior, target occupancy, distance to platform, and initial swimming angle.

### MEC Calcium Imaging

Folders:

- [`Calcium imaging/MEC/Reelin`](Calcium%20imaging/MEC/Reelin/README.md)
- [`Calcium imaging/MEC/CalB`](Calcium%20imaging/MEC/CalB/README.md)

These folders contain MATLAB and Python scripts for MEC calcium-imaging analyses, including animal tracking, cell detection, grid-cell identification, grid-module analysis, spatial modulation, border selectivity, head-direction tuning, speed modulation, phase-difference analysis, synchronization analysis, and decoding.

Main outputs include processed tracking files, calcium trace files, rate maps, grid-score summaries, module assignments, spatial-information metrics, head-direction metrics, synchronization metrics, decoding summaries, and manuscript-related figure outputs.

### CA1 Calcium Imaging

Folder: [`Calcium imaging/CA1`](Calcium%20imaging/CA1/README.md)

This folder contains MATLAB scripts for CA1 calcium-imaging analyses, including behavioral tracking, cell validation, place-cell analysis, and Bayesian decoding.

Main outputs include refined animal trajectories, validated calcium traces, place-cell metrics, rate maps, cleaned analysis data, Bayesian decoding results, and summary figures.

## Requirements

Most scripts are written in MATLAB. Some decoding analyses are written in Python.

Recommended software:

- MATLAB with commonly used analysis, statistics, image-processing, and signal-processing functions.
- Python 3 for decoding scripts in `Calcium imaging/MEC/Reelin/11_Decoding`.
- Inscopix Data Processing Software for preprocessing calcium-imaging data before running the MATLAB analysis scripts.
- ImageJ/Fiji for ROI preparation in histology analyses when ROI files are required.

Several scripts contain placeholder or hard-coded input paths. Update these paths before running the analysis on new data.

## Basic Usage

1. Open the target analysis folder in MATLAB or Python.
2. Read the `README.md` in that folder.
3. Update input paths, sample names, group definitions, and animal IDs as needed.
4. Run the scripts in the order listed in the corresponding folder-level README.
5. Check the generated `.mat`, `.xlsx`, `.csv`, `.jpg`, `.pdf`, `.emf`, or `.fig` outputs described in the Expected outcome section.

## Demo Data

Small demo data are included in selected folders:

- `Histrogy/Quantitative analysis of CalB patches/demo`
- `Histrogy/The average signal along layer/demo`
- `Histrogy/Iba1 and GFAP/demo`
- `MWM/demo`
- `Calcium imaging/MEC/Reelin/01_Each_animal/demo`

These demo files are intended to help reviewers inspect input formats and test representative parts of the analysis workflow.

## Notes For Reviewers

This repository is provided for peer review and reproducibility assessment of the associated unpublished manuscript. Some scripts were written for the manuscript-specific directory structure and may require path edits before use on another computer. Folder-level README files describe the expected output of each analysis step.

## Citation

Please cite the final published article when it becomes available. Until publication, refer to the manuscript title:

Yamamoto N, Osanai H, Ramesh K, Nair IR, Henkemeyer M, Ogawa SK, Kitamura T. **EphrinB2-dependent calbindin patch architecture shapes discrete grid modularity in the medial entorhinal cortex**. Manuscript under peer review.
