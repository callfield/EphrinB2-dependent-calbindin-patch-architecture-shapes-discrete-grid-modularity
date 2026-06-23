close all
clear; 
CDir=pwd;
cd(CDir)
addpath("function")
%%  

[num, txt] = xlsread('results/peak_center.xlsx');
dataDir = '../intensity/';
Snum = length(txt); % number of samples

filenameS = cell(1,Snum);
filenameTif = cell(1,Snum);
filenameROI = cell(1,Snum);

for fi=1:1:Snum
    filenameS{fi} = txt{fi};
    filenameTif{fi} = strcat(erase(txt{fi}, ".csv"), ".tif");
    
    dataTif_raw{fi} = imread([dataDir, char(filenameTif{fi})]);
    dataTif_adj{fi} = imadjustn(dataTif_raw{fi},stretchlim(dataTif_raw{fi}(:),[0.001 0.999]) ) % Optimize contrast
    filenameROI{fi} = strcat(erase(txt{fi}, ".csv"), ".roi");
    dataROI{fi} = ReadImageJROI([dataDir,  char(filenameROI{fi})]);
    ROIcoordinates{fi} = dataROI{fi}.mnCoordinates;
    linetemp{fi} = dataROI{fi}.mnCoordinates;
    
    filenamebgROI{fi}= strcat(erase(txt{fi}, ".csv"), ".bg.roi");
    bgROI{fi} = ReadImageJROI([dataDir,  char(filenamebgROI{fi})]);
    bgROIcoordinates{fi} = bgROI{fi}.mnCoordinates;
    bglinetemp{fi} = bgROI{fi}.mnCoordinates;
end



%% Roi drowing
dv = 0; % analysis range of dorsoventral axis
sd = 250; % analysis range of superficial-deep axis
pixelsize = 0.65; % um/pix
winsize = 115; % gaussian filter window(pixcel)
gw = gausswin(winsize);
sgw = sum(gw)


P = cell(20,Snum);L = cell(20,Snum);W = cell(20,Snum);PV = cell(20,Snum);
P(1,:) = filenameS; L(1,:) = filenameS; W(1,:) = filenameS; PV(1,:) = filenameS ;
meanP = cell(1,Snum)    

for fi=1:1:Snum
 %% Backgroung    
    Plocs = rmmissing(num(:,fi))   

         % interpolation       
    x = bglinetemp{fi} (:,1);
    y = bglinetemp{fi} (:,2);
    dist1 = sum(sqrt(diff(x).^2+diff(y).^2));
    N = 1000;
    gradient = [0 0; 0 0];
    [bgyy, bgxx] = cubicSpline2d(N,[x, y],gradient); % cubic spline interporation

    % draw ROI
    start_points = [bgxx; bgyy];
    start_points = start_points(:, 1:end-1);
    goal_points = [bgxx; bgyy]; 
    goal_points = goal_points(:, 2:end);
    v = goal_points - start_points;
    B = cumsum(sqrt(diff(bgxx).^2+diff(bgyy).^2));
    bgst = knnsearch(B',500/0.65);
    bged = knnsearch(B',1500/0.65);
    meanBG = improfile(dataTif_raw{fi},[bgxx(bgst) bgxx(bged)],[bgyy(bgst) bgyy(bged)]);


    % um to pix
    locs2 = Plocs/pixelsize;

    % interpolation       
    x = linetemp{fi} (:,1);
    y = linetemp{fi} (:,2);
    dist1 = sum(sqrt(diff(x).^2+diff(y).^2));
    N = 1000;
    gradient = [0 0; 0 0];
    [yy, xx] = cubicSpline2d(N,[x, y],gradient); % cubic spline interporation

    % draw ROI
    start_points = [xx; yy];
    start_points = start_points(:, 1:end-1);
    goal_points = [xx; yy]; 
    goal_points = goal_points(:, 2:end);
    v = goal_points - start_points;
    v = sd*normc(v);
    deep = (sd+150)*normc(v);
    xo1 = xx(1:end-1)+v(2,:); xo1 = [xo1, xx(end)+v(2,end)];
    xo2 = xx(1:end-1)-deep(2,:); xo2 = [xo2, xx(end)-deep(2,end)];
    yo1 = yy(1:end-1)-v(1,:); yo1 = [yo1, yy(end)-v(2,end)];
    yo2 = yy(1:end-1)+deep(1,:); yo2 = [yo2, yy(end)+deep(2,end)];

    %Peak Mark &  Visualise analized region
       subplot('Position', [0.02 0.02 0.35 0.9])
    image(dataTif_adj{fi});
    axis image
    axis off
    ylim([0 size( dataTif_adj{fi},1)])
    xlim([0 size( dataTif_adj{fi},2)])
    hold on
    B = cumsum(sqrt(diff(xx).^2+diff(yy).^2));
    Idx = knnsearch(B',locs2);
    plot(xx(Idx),yy(Idx),'pw', 'MarkerFaceColor','w',...
    'MarkerSize',2)
    title(strcat(erase(filenameS{fi}, ".csv")));

    st = knnsearch(B',locs2 - dv)
   
        for i = 1:1:length(locs2)
        plot([xo1(st(i)), xo2(st(i))],[yo1(st(i)), yo2(st(i))], 'w:');
          end
    
     plot([bgxx(bgst) bgxx(bged)],[bgyy(bgst) bgyy(bged)], 'b--o');    
        
    if isempty(Plocs)  
    else
        meanP = cell(1,length(locs2))

        for i = 1:1:length(locs2)
            sz = size(improfile(dataTif_raw{fi},[xo2(Idx(i)) xo1(Idx(i))],[yo2(Idx(i)) yo1(Idx(i))]));
            c2 = zeros(sz)  ;
            n = length(improfile(dataTif_raw{fi},[xo2(Idx(i)) xo1(Idx(i))],[yo2(Idx(i)) yo1(Idx(i))]));
            a = locs2(i)-dv : 1 :locs2(i)+dv ;
            idx = knnsearch(B',a(:) );

                for k = 1
                c = improfile(dataTif_raw{fi},[xo2(idx(k)) xo1(idx(k))],[yo2(idx(k)) yo1(idx(k))],n) ;
                c2 = c2 + c;
                end


            meanP{i} = c2
          

                bsuv =  imdivide(meanP{i}, mean2(meanBG));
             ga = filter(gw, sgw, bsuv);
   
            % Display mean singal intensity around peaks  
            subplot(length(locs2),3, 2 + (i-1)*3)
            plot(meanP{i}(:,:,1))
              [pks,locs,w,p_from_varey] =  findpeaks_ho(meanP{i}(:,:,1), 'Annotate','extents','SortStr', 'descend', 'NPeaks', 1)% save variables for prominence analysis
            findpeaks(meanP{i}(:,:,1),'Annotate','extents', 'SortStr','descend', 'NPeaks', 1) % plot prominence analysis
          
            legend('hide');



            [pks,locs,w,p_from_varey] =  findpeaks_ho(ga(:,:,1), 'Annotate','extents','SortStr', 'descend', 'NPeaks', 1)% save variables for prominence analysis

            subplot(length(locs2),3,   3 + (i-1)*3);
            plot(ga(:,:,1))
            findpeaks(ga(:,:,1),'Annotate','extents', 'SortStr','descend', 'NPeaks', 1) % plot prominence analysis
            text(locs, pks+0.5 ,num2str(round(p_from_varey,3,'significant')),'FontSize',5)
            text(locs, pks-1,num2str(round(w,3,'significant')),'FontSize',5)
            legend('hide');
            ylim([3 15]);
            if isempty(pks)
                P(1+i, fi) = num2cell(0);
                L(1+i, fi) = num2cell(0);
                W(1+i, fi) = num2cell(0);
                PV(1+i, fi) = num2cell(0);
            else
                P(1+i, fi) = num2cell(pks);
                L(1+i, fi) = num2cell(locs);
                W(1+i, fi) = num2cell(w);
                PV(1+i, fi) = num2cell(p_from_varey);
            end
       end

        set(gcf,'color',[1 1 1]);
        set(gcf,'InvertHardCopy','off');
        subplot(length(locs2),3, 2);
        title("Intnsity plot (SD-axis)");
        subplot(length(locs2),3, 3);
        title("G filtered Signal Ratio (Int/mean Int)");
       end  

    close all

    end
     
xlswrite("results/sng_SD-axis_pks.xlsx", P);
xlswrite("results/sng_SD-axis_locs.xlsx", L);
xlswrite("results/sng_SD-axis_width.xlsx", W);
xlswrite("results/sng_SD-axis_peaks_from_varey.xlsx", PV);


