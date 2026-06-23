close all
clear; 

CDir=pwd;
cd(CDir)
addpath("function")
mkdir results
%% Reading csv files
[file, path] = uigetfile('*.xlsx','select List file (1)'); % sample name file
[num, txt] = xlsread(file);
dataDir = 'PATH/FOR/CSV/AND/IMAGE/AND/ROI';
Snum = length(txt); % number of samples

% reading signal intensity files
rStr = 1 ; % Reading start (pixel)
rEnd = 2566; % Reading end (pixel)
figxlim = [200 1700]; 
figylim = [0 80];
filenameS = cell(1,Snum);
filenameBG = cell(1,Snum);
dataS = cell(1,Snum);
dataBG = cell(1,Snum);

filenameS = txt(:,1);
filenameBG = txt(:,2);
for fi=1:1:Snum
    dataS{fi} = csvread( [dataDir, char(filenameS{fi})] ,1,0);
    dataS2{fi} =  dataS{fi}(:,:); % use for ROI drowing
    dataS{fi} =  dataS{fi}(rStr:rEnd,:);
    dataBG{fi} = csvread( [dataDir, char(filenameBG{fi})] ,1,0);
    dataBG{fi} =  dataBG{fi}(rStr:rEnd,:);
end



%% background 
 
I_SBratio = intergroup_SBratio(dataS, dataBG, Snum, figxlim, [0 5], filenameS); %signa l-to-background ratio


%% smooth by gaussian filter

winsize = 115; % gaussian filter window(pixcel)
I_gaussianfilt = gaussian_filter_ny(I_SBratio, winsize, Snum, figxlim, [0 5], dataS, filenameS);


cutlength = 100; % cut first pixcels to ignore the effect of gaussian window.
                               
minpeakprominence = 0.5; % set minimum height of intensity prominence.

[pks, locs, w, p_from_varey, wxPk] = prominence_ny(I_gaussianfilt, Snum, figxlim, [0 5], dataS, cutlength, minpeakprominence,filenameS); % prominence analysis (peak intensity etc.)


P = cell(20,Snum);L = cell(20,Snum);W = cell(20,Snum);PV = cell(20,Snum);C = cell(20,Snum);
P(1,:) = filenameS; L(1,:) = filenameS; W(1,:) = filenameS; PV(1,:) = filenameS ; C(1,:) = filenameS ;
for k=1:1:Snum
    P(2:1+length(pks{k}), k) = num2cell(pks{k});
    L(2:1+length(pks{k}), k) = num2cell(locs{k});
    W(2:1+length(pks{k}), k) = num2cell(w{k});
    PV(2:1+length(pks{k}), k) = num2cell(p_from_varey{k});
    C(2:1+length(pks{k}), k) = num2cell( (  wxPk{k}(1: length(pks{k})) + wxPk{k}(length(pks{k})+1:length(pks{k})*2)  )/2);
end
xlswrite("results/pks.xlsx", P);
xlswrite("results/locs.xlsx", L);
xlswrite("results/width.xlsx", W);
xlswrite("results/peaks_from_varey.xlsx", PV);
xlswrite("results/peak_center.xlsx", C);



%% Roi drowing
filenameTif = cell(1,Snum);
filenameROI = cell(1,Snum);
filenameGfil = cell(1,Snum);
filenamePeak = cell(1,Snum);
filenameRaw = cell(1,Snum);
filenameDraw =cell(1,Snum);
[num, txt] = xlsread("results/peak_center.xlsx");
pixelsize = 0.65; % um/pix


for fi=1:1:Snum
    filenameTif{fi} = strcat(erase(filenameS{fi}, ".csv"), ".tif");
    dataTif_raw{fi} = imread([dataDir, char(filenameTif{fi})]);
    dataTif{fi} = imadjustn(dataTif_raw{fi},stretchlim(dataTif_raw{fi}(:),[0.001 0.999]) ) ;% Optimize contrast
    filenameROI{fi} = strcat(erase(filenameS{fi}, ".csv"), ".roi");
     dataROI{fi} = ReadImageJROI([dataDir,  char(filenameROI{fi})]);
    ROIcoordinates{fi} = dataROI{fi}.mnCoordinates;
    linetemp{fi} = dataROI{fi}.mnCoordinates;
end





for j=1:3:Snum
    [ha, pos] = tight_subplot(1, 3, [.0 .0], [.01 .1], [.01 .01]) %tight_subplot(Nh, Nw, gap, marg_h, marg_w)
    for ii=1:1:3
        fi = ii+j-1;
        %subplot(1,3,ii)
        axes(ha(ii))
        image(dataTif{fi});
        axis image
        axis off
        ylim([0 size(dataTif{fi},1)]);
        xlim([0 size(dataTif{fi},2)]);

        % interpolation
        hold on
        x = linetemp{fi} (:,1);
        y = linetemp{fi} (:,2);

        dist1 = sum(sqrt(diff(x).^2+diff(y).^2));
        %N = 1000;
        N= length(dataS2{fi});
        gradient = [0 0; 0 0];
        [yy, xx] = cubicSpline2d(N,[x, y],gradient); % cubic spline interporation

        % draw ROI
        ml = 100;
        start_points = [xx; yy];
        start_points = start_points(:, 1:end-1);
        goal_points = [xx; yy]; 
        goal_points = goal_points(:, 2:end);
        v = goal_points - start_points;
        v = ml*normc(v);
        xo1 = xx(1:end-1)+v(2,:); xo1 = [xo1, xx(end)+v(2,end)];
        xo2 = xx(1:end-1)-v(2,:); xo2 = [xo2, xx(end)-v(2,end)];
        yo1 = yy(1:end-1)-v(1,:); yo1 = [yo1, yy(end)-v(2,end)];
        yo2 = yy(1:end-1)+v(1,:); yo2 = [yo2, yy(end)+v(2,end)];

       

         
         plot(xo1(100:rEnd-1), yo1(100:rEnd-1), 'w--');
         plot(xo2(100:rEnd-1), yo2(100:rEnd-1), 'w--');
         plot([xo1(100), xo2(100)],[yo1(100), yo2(100)], 'w--');
         plot([xo1(rEnd-1), xo2(rEnd-1)],[yo1(rEnd-1), yo2(rEnd-1)], 'w--'); %avoid line end; line ends tend to have a irregular slope

        
        % um to pix
        locs2 = locs{fi}/pixelsize;
        w2 = w{fi}/pixelsize;
        c2 = rmmissing(num(:,fi))/pixelsize;

        %Peak Mark
        B = cumsum(sqrt(diff(xx).^2+diff(yy).^2));
        Idx = knnsearch(B',locs2);
        % plot(xx(Idx),yy(Idx),'pw', 'MarkerSize',7)
        plot(xo1(Idx),yo1(Idx),'pw', 'MarkerFaceColor','w',...
            'MarkerSize',3)
        title(strcat(erase(filenameS{fi}, ".csv")));
 
        %Plot peak width
        wl =  c2 - w2/2;
        wh =  c2 + w2/2;
        Idx_wl = knnsearch(B', wl);
        Idx_wh = knnsearch(B', wh);
              
        for k = 1:1:length(Idx_wh)
           plot([xx(Idx_wl(k)), xx(Idx_wh(k))], [yy(Idx_wl(k)), yy(Idx_wh(k))], '-w') ;
        end

    end
    set(gcf,'color',[1 1 1]);
    set(gcf,'InvertHardCopy','off');
    
    print(strcat("results/Roi_", erase(filenameS{j}, ".csv"), ".jpg"), '-djpeg', '-r500')
    close all    
end


for fi=1:3:Snum
    filenameGfil{fi} = strcat("results/Gfiltered_SBratio ",erase(filenameS{fi}, ".csv"), ".jpg");
    dataGfil{fi} = imread(filenameGfil{fi});
    filenamePeak{fi} = strcat("results/peak ",erase(filenameS{fi}, ".csv"), ".jpg");
    dataPeak{fi} = imread( filenamePeak{fi});
    filenameRaw{fi} = strcat("results/raw ",erase(filenameS{fi}, ".csv"), ".jpg");
    dataRaw{fi} = imread( filenameRaw{fi});
    filenameDraw{fi} = strcat("results/Roi_",erase(filenameS{fi}, ".csv"), ".jpg");
    dataDraw{fi} = imread( filenameDraw{fi});
end



  
for fi=1:3:Snum
     [ha, pos] = tight_subplot(2, 2, [0 0], [.01 .01], [.01 .01]) %tight_subplot(Nh, Nw, gap, marg_h, marg_w)
    axes(ha(1))
    image(dataDraw{fi})
    axis off
    axes(ha(2))
    image(dataRaw{fi})
    axis off
    axes(ha(3))
    image(dataPeak{fi})
    axis off
     axes(ha(4))
    image(dataGfil{fi})
    axis off
    % print(strcat("Sammary_", erase(filenameS{fi}, ".csv"), ".jpg"), '-djpeg', '-r800')
    close all   
end

