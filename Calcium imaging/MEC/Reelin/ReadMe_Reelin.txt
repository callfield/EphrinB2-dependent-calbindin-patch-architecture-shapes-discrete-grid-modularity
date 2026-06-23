Small damo data (IDPS processed csv and bahavior track date aquired with OF_analysis_PCI_1.m)included ＠MEC\demo

Animal traking 
@01_Each_animal
1, Run tracking code “OF_analysis_PCI_1.m” with matlab


Cell detection and rate map
Calcium imaging (process with Inscopix Data Processing Software )


Identification of grid modules
@01_Each_animal
1, Run “OF_analysis_PCI_2.m” with matlab 
2, Run “OF_analysis_PCI_3.m” with matlab 


Summarize Grid analysis 
1, Run “a_MakeDirFile.m” with matlab 
2, Run “b_Run_GridBregma.m” with matlab 
3, Run “c_MakeTable.m” with matlab 

Spatial modulation of calcium activity
@ 02_Spatial_cell
1, Run “Run_conditional_entropy_Zscore.m” with matlab 


Border selectivity 
@ 03_Border
1, Run “Run_BorderScore_v2.m” with matlab 



Head-direction tuning
@ 04_Head_Direction
1, Run “Run_a_direction_Analysis_nofiler_Raylign.m” with matlab 
2, Run “Run_b_nofilterHD_pilot_RL.m” with matlab 



Grid_Scale
@ 05_G_Scale
1, Run “Run_Scatter_scale.m” with matlab 

Grid_Field
@ 06_G_Field
1, Run “Run_Scatter_field.m” with matlab 

Grid_Orientation
@ 07_G_Orient
1, Run “Run_Scatter_ori.m” with matlab 



Identification of grid modules
@ 08_GridModule
1, Run “Run_a_Eachtrial_GridScaleRatio_2D.m” with matlab 
2, Run “Run_b_Define_2DModule.m” with matlab 
3, Run “Run_c_ksdNorm_scale_colorscatter.m” with matlab 
4, Run “Run_d_Min_Gscale.m” with matlab 


Speed modulation analysis between Grid module
@ 09_Grid_module_speed
1, Run “Run_a_Speed_vs_meandFCorr2_NewP.m” with matlab 
2, Run “Run_b_Speed_zScore_Hz_each1.m” with matlab 
3, Run “Run_c_Speed_zScore_Hz_each5.m” with matlab 
4, Run “Run_d_zSpeedeach1_eachModule.m” with matlab 
5, Run “Run_e_zSpeedeach1_eachModule_Linfit.m” with matlab 


Spatial specificity of population activities of grid cells
@ 10_SpatialUniqness
1, Run “Run_a_SpatialUniquness_FullVec_Grid.m” with matlab 
2, Run “Run_b_Fig_SpatialUniquness_Grid.m” with matlab 



Animal position decoding from non-grid spatial cells
@ 11_Decoding
1, Run “Run_GridPhaseDiff.m” with matlab 
2, Run “b_decoding batch bestSpatialInfo v5.1 Spatial_nogrid.py" with python 
3, Run “c_decoding dist Summary.py" with python 
4, Run “d_decoding dist Summary 2.py" with python 


Grid_PhaseDifference
@ 12_G_Orient
1, Run “Run_GridPhaseDiff.m” with matlab 


Grid scale-field width relationship and Module distance
@ 13_otherGridFigures
1, Run "a_GridScore.m” with matlab for Ex Fig1 and 4
2, Run "b_FiringRate.m” with matlab for Ex Fig1 and 4
3, Run "c_GridScale.m” with matlab for Ex Fig1 and 4
4, Run "d_DVpositionVsRawGridScale.m” with matlab for Ex Fig1 and 4
5, Run "e_GridScaleVsWidth.m” with matlab for Ex Fig4
6, Run "f_Fig_ScalePerWidth.m” with matlab for Ex Fig4
7, Run "g_SpeedModRate.m” with matlab for Ex Fig4
8, Run "h_GridPhaseDistribution.m” with matlab for Ex Fig5
9, Run "i_GridOri.m” with matlab for Ex Fig5


