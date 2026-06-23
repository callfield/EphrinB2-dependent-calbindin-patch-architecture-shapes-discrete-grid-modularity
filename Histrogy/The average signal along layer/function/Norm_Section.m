function Norm_signal = Norm_Section(nSec,dataDir,Name)



   
    rowData=csvread( [dataDir, char(Name)] ,1,1);


    mecDepth = length(rowData);
    ttt = mecDepth/nSec:mecDepth/nSec:mecDepth;
    Sectioned_rowData= interp1(1:1:mecDepth, rowData,ttt , 'pchip');
    tmp = Sectioned_rowData - min(Sectioned_rowData);
    Norm_signal =  tmp/max(tmp);
    Norm_signal=movmean(Norm_signal,5);
