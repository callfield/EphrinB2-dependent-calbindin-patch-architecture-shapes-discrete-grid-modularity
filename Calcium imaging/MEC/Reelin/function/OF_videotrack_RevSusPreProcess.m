function trk = OF_videotrack_RevSusPreProcess(vidDir,vidFileName,  vw, vh,Thresould)

mkdir Track
% vidDir = '' ;
% vidFileName{1} = 'Trial_1.avi';
% fr = 10 % Frame rate number for video analyse;
% vw = 100 % acrutal arena width(cm)
% vh = 100 % acrutal arena high(cm)
% Thresould=160;

vNum = length(vidFileName);
trk = cell(1,vNum); 

    for i=1:1:vNum;

        vidObj{i} = VideoReader( [vidDir, char(vidFileName{i})] );
        vidObj{i}
        
        numFrames=vidObj{i}.NumFrames
        fm = vidObj{i}.FrameRate;
        vidHeight = vidObj{i}.Height;
        vidWidth = vidObj{i}.Width;
        s = struct('cdata',zeros(vidHeight,vidWidth));
        


         vidObj{i}.CurrentTime = 0;
         dt = 1/fm; k = 1;
         tmp_video=zeros(vidHeight,vidWidth);
         while hasFrame(vidObj{i})
            s(k).cdata = rgb2gray(readFrame(vidObj{i}));

                img=im2double(s(k).cdata);% comvert to grayscale to matrix
                img=abs(img-1);%imvert
                tmp_video=tmp_video+img; % sum for later averaging

            k = k+1;
         end 
         Ave_Img =tmp_video/numFrames;


     % create thebinalised  video  with defined fps
     writerObj = VideoWriter( strcat(vidDir, "Bined_",char(vidFileName{i})));

      writerObj.FrameRate = fm;
       % open the video writer
       open(writerObj);
        M=zeros(numFrames,1); area_array=cell(numFrames,1);biarea=cell(numFrames,1);
        susp0=[];
        center = zeros(numFrames,2);
        for k=1:1:numFrames
            
            img=abs(im2double(s(k).cdata)-1); %comvert to grayscale to matrix and imvert
            img2=(img-Ave_Img); % subtract from image
            I2bi=mat2gray(img2)*255; % covert to matrix to grayscale image
            I2bi=medfilt2(I2bi);

            I2bi(I2bi<Thresould)=0;%made binary image
            I2bi(I2bi>=Thresould)=1; %made binary image


            %tracking maximum area's particle
            bwl = bwlabel(I2bi); %label binary object
            biarea{k} = regionprops(bwl, 'Area', 'BoundingBox', 'Centroid'); %#ok<MRPBW>
            if  isempty(vertcat(biarea{k}.Area)) ==0     
                area_array{k} = vertcat(biarea{k}.Area);    
                M(k,1) = max(area_array{k}); %find maximum area
           % Mark position of centroid
                c=vertcat(biarea{k}.Centroid);
                id=find(area_array{k}==max(area_array{k}));
                x=round(c(id(1),1));% sometimes there are same maxima
                y=round(c(id(1),2));
                xx=x-5:x+5;yy=y-5:y+5;
                xx=xx(xx>0&xx<vidWidth);
                yy=yy(yy>0&yy<vidHeight);
                I2bi(yy,xx)=0.7;
                center(k,:) =c(id(1),:); % centroid of closest area  
           
  %  Deifinition of suspecious frame #1 two object
                if length(area_array{k})>1& min(maxk(area_array{k},2))>100
                 susp0 = [susp0;k]; 
                end
                        
               end
            %imshow(I2bi)
            writeVideo(writerObj, I2bi);
        end
        close(writerObj);% close the video writer
        


%  Deifinition of suspecious frame #2 too big
        ss=find(M>1500);
        susp0=[susp0;ss];
        susp0=unique(susp0);
        
        
%  Deifinition of suspecious frame #3 rescue some frame (not too big&not too fast)
        psw = vw / vidObj{i}.Width; %[cm/pix]
        psh = vh / vidObj{i}.Height; %[cm/pix]
        
        x = center(:,1);
        y = center(:,2);
        spx = diff(x)/dt; spx = [0; spx]; %pix/s
        spy = diff(y)/dt; spy = [0; spy]; %pix/s
        spy = spy*psh; %cm/s
        spx = spx*psw; %cm/s
        v = sqrt(spx.^2 + spy.^2); 
       
        susp2=[];tmp_M=M;tmp_M(susp0)=[];
        for k=susp0.'    
            if k>1&k<numFrames
                if  sum(maxk(area_array{k},2))<prctile(tmp_M,95)&v(k)<40&v(k-1)<40&v(k+1)<40
                    susp2=[susp2;k];
                end
            else
                if  sum(maxk(area_array{k},2))<prctile(tmp_M,95)&v(k)<40
                    susp2=[susp2;k];
                end
            end
        end
       susp=setdiff(susp0,susp2);
        

 

%% Remove detection errer frame         
        
       susp=[susp;susp+1];% add frame affected (to cal velocity) by suspecious frame
       susp=unique(susp);
       susp=susp(susp<=numFrames);
       
        v(susp) = zeros(length(susp),1);
        trk{i} = [x*psw, y*psh, v];
        trk{i}(susp,:)=zeros(length(susp),3);
        
        tq = linspace(0, vidObj{i}.Duration-dt, numFrames);
        tq = tq.';
        trk{i} = [tq, trk{i}];

   %     csvwrite(strcat("Track\Track ", erase(vidFileName{i}, [".avi"]), ".csv"), trk{i});

        
%% visualise

        h = plot(trk{i}(:,1), trk{i}(:,4), 'k');
        xlabel('time [sec]')
        ylabel('velocity [cm/sec]')
        print(strcat("Track\Velocity ", erase(vidFileName{i}, ["converted_", ".avi"]), ".jpg"), '-djpeg', '-r0')
        close
        
        colormap(jet)
        xx=x*psw; xx(susp)=[];yy=y*psh;yy(susp)=[];vv=v;vv(susp)=[];
        cplot(xx,yy,vv,'LineWidth',3);
        ylim([0 vh]);
        xlim([0 vw]);
        colorbar;
        print(strcat("Track\Track ", erase(vidFileName{i}, ["converted_", ".avi"]), ".jpg"), '-djpeg', '-r0')
        close
   
        % check
        tmp_M=M;
        tmp_M(susp)=[];
        histogram(tmp_M, 150)
        xlabel('Area (pic)')
        ylabel('number of frame')
         print(strcat("Track\Nosusp_Area_", erase(vidFileName{i}, ["converted_", ".avi"]), ".jpg"), '-djpeg', '-r0')
        close

    
    end
    
    
    
end

