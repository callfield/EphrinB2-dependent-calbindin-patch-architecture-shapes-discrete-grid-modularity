function [Rec,RecStart] = OF_videotrack_lampPreProcess(vidDir2,LvidFileName,  vw, vh)
%%
mkdir Track

vNum = length(LvidFileName);
Rec = cell(1,vNum); 
RecStart = cell(1,vNum);

    for i=1:1:vNum;

        vidObj{i} = VideoReader( [vidDir2, char(LvidFileName{i})] );
        vidObj{i}

        fm = vidObj{i}.FrameRate;
        vidHeight = vidObj{i}.Height;
        vidWidth = vidObj{i}.Width;
        s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
            'colormap',[]);

         vidObj{i}.CurrentTime = 0;
         dt = 1/fm; k = 1;
        while hasFrame(vidObj{i})
                s(k).cdata = rgb2gray(readFrame(vidObj{i}));
            k = k+1;
        end

        Min_Img=abs(im2double(s(1).cdata)-1);%imvert
        mm=80
        k = 1;Rec{i} = [];

        while k<=length(s)

             I2bi=s(k).cdata ; % covert to matrix to grayscale image
             I2bi(s(k).cdata<180)=0;%made binary image
             I2bi(s(k).cdata>=180)=1; %made binary image
           
           
            
            %tracking maximum area's particle
            bwl = bwlabel(I2bi); %label binary object
            biarea = regionprops(bwl, 'Area', 'BoundingBox', 'Centroid'); %#ok<MRPBW>
            area_array = vertcat(biarea.Area);
           M = max(area_array(area_array>mm)); %find maximum area     
            if  isempty(M) ==0       
               Rec{i} = [Rec{i};k];
             end
   
            k = k+1;
        end

        % find file showing rec start end
         RecStart{i}=[]; RecStart{i}(1,1) = Rec{i}(1,1); t =1;k=2;
        for k = 2:1:length(Rec{i})
            if Rec{i}(k)-1 ~= Rec{i}(k-1)
                 RecStart{i}(t+1,1) = Rec{i}(k);
                 RecStart{i}(t,2) = Rec{i}(k-1);
                 t=t+1;
            end
        end
        [row col]=size(RecStart{i});
       RecStart{i}(row,2) = Rec{i}(end);
       RecStart{i}(:,3) = RecStart{i}(:,2)-RecStart{i}(:,1)+1;
        
        
        
    end