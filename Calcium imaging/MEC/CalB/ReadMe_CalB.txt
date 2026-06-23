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
1, Run “MakeDirFile.m” with matlab 
2, Run “Run_GridBregma.m” with matlab 


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



Synchronization analysis of intra- and trans-CalB+ patch cells
@ 05_Synchronization
1, Run “Run_a_Ex_DV_meandF_Island.m” with matlab 
2, Run “Run_b_SubCorrelation_Island.m” with matlab 
3, Run “Run_c_Allcellpair_summary.m” with matlab 









