close all
clear; 
path('function')
nSec=50;
CDir=pwd;
dataDir = 'data/for/csv'; %intensity files

Norm_intensity=Norm_Section(nSec,dataDir, "datafile.csv");

