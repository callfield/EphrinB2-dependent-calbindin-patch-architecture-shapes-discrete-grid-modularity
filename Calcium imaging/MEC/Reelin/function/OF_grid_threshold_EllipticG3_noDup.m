function [Gth, Rand_Grid_Score]= OF_grid_threshold_EllipticG3_noDup(RandNum, lk,vlim, Trk, GSoccup_m, caFr, numFrames, numCells,  vw, vh, bin, win,w,h)

w_binNum=ceil(vw/2);
h_binNum=ceil(vh/2);
c_row = w_binNum-1;c_col = h_binNum-1;
%Transform Cartesian coordinates to Polar coordinates
[oy, ox] = ndgrid((w_binNum-1):-1:-(w_binNum-1),-(h_binNum-1):1:(h_binNum-1)); % make Cartesian coordinates, be careful
[theta,rho] = cart2pol(ox,oy);

ran_lk=cell(numCells,1);
activ_m=cell(numCells,1);GSactiv_m=cell(numCells,1);GSrate_map=cell(numCells,1);
Rand_Grid_Score=zeros(numCells,RandNum);

Fpoint=zeros(numFrames,numCells);
for k=1:1:numCells 
        Fpoint(round(lk{k}*caFr),k)=1;
end


%%
 for rrtt = 1:1:RandNum 
rrtt
   ran_Fpoint=circshift(Fpoint,randsample(numFrames,1));
    for k=1:1:numCells 
       find(ran_Fpoint(:,k)==1);
       x=Trk(find(ran_Fpoint(:,k)==1),2);y=Trk(find(ran_Fpoint(:,k)==1),3);v=Trk(find(ran_Fpoint(:,k)==1),4);
    %    ran_lk{k} = sort(randsample(numFrames,length(lk{k})));
    %    x=Trk(ran_lk{k},2);y=Trk(ran_lk{k},3);v=Trk(ran_lk{k},4);
        x=x(find(v>vlim));y=y(find(v>vlim));
        ran_lk{k}=length(x);
        activ_m{k}=zeros(ceil(vw/bin),ceil(vh/bin));
        for i=bin:bin:vw
            for l=bin:bin:vh
                dd=sqrt(((x-i).^2)+((y-l).^2));
                activ_m{k}(ceil((vh-l)/bin)+1,ceil(i/bin))=length(find(dd<=win));
            end
        end
        GSactiv_m{k} = imgaussfilt(activ_m{k},2,'FilterSize', 5);
        GSrate_map{k} = GSactiv_m{k}./(GSoccup_m/caFr);
    end
    
    Autoc = cell(numCells,1);  
    for k=1:1:numCells ;
    %        if ran_lk{k} > 80 ; ;
                for dx=-w_binNum:1:w_binNum
                    for dy=-h_binNum:1:h_binNum
                        aut1=zeros();aut2=zeros();aut3=zeros();aut4=zeros();aut5=zeros();aut6=zeros();aut7=zeros();
                        n=0;
                           for x= 1:1:w_binNum
                               if 0<(x + dx) && (x + dx)*bin <=vw
                                   for y=1:1:h_binNum
                                      if 0<(y+dy) && (y+dy)*bin <= vh
                                          %Becareful!! (y,x) and y is start from upper left
                                          r1 = GSrate_map{k}(y, x);
                                          r2 = GSrate_map{k}(y+dy, x+dx);
                                        aut1 = aut1 + r1*r2 ;
                                        aut2 = aut2 + r1 ;
                                        aut3 = aut3 + r2 ;
                                        aut4 = aut4 + r1^2 ;
                                        aut5 = aut5 + r1 ;
                                        aut6 = aut6 + r2^2 ;
                                        aut7 = aut7 + r2;    
                                        n=n+1 ;
                                      end    
                                   end
                               end
                           end
                       if n>20
                           Autoc{k}(dy+h_binNum, dx+w_binNum) = (n*aut1-aut2*aut3)/(sqrt(n*aut4 - aut5^2) * sqrt(n*aut6-aut7^2)) ;                    
                       end

               end
            end

     %   end
    end

    Dis_Center_Autoc=cell(numCells,1);%gaus_Autoc{k}=cell(numCells,1);
    in_r=cell(numCells,1);
    T_Autoc=cell(numCells,1);
    rot_T_Autoc=cell(numCells,6);
    Grid_Score=NaN(numCells,1);
    trim_Autoc=cell(numCells,1);
    rGS=zeros(numCells,round(vh));

    for k=1:1:numCells 
    %     if ran_lk{k} > 80 ;   
         B = imgaussfilt(Autoc{k},3); % smoothed
        [cent, varargout]=FastPeakFind(B*255,0); %find peaks in autocorrecogram
      %  [cent, varargout]=FastPeakFind(-1*Autoc{k}*255,0); ; %find peaks in autocorrecogram
        Angle=reshape(theta.*varargout,[],1);
        Dist=reshape(rho.*varargout,[],1);
        tmp=[Angle, Dist];
        G_all_peaks=tmp(tmp(:,1)~=0&tmp(:,2)~=0,:);
    if length(G_all_peaks)>=6
        %G_all_peaks=tmp(find(tmp(:,2)>=10),:); % Remove peak around (<10cm) to center
        G_6peaks=G_all_peaks(knnsearch(G_all_peaks(:,2),0,'K',6),:);
        in_r{k}=bin*min(G_6peaks(:,2))-10;

    %Correct Elliptical distortions (Brandon et al., Science, 2011)
        G_3peaks=G_6peaks(G_6peaks(:,1)>0,:);
        r1=min(G_3peaks(:,2));
        R1=max(G_3peaks(:,2));
        min_Peak=G_3peaks(G_3peaks(:,2)==r1,:);
        max_Peak=G_3peaks(G_3peaks(:,2)==R1,:);
      % case1: minor axis=closet peak
        c = max_Peak(1,1)-min_Peak(1,1);%the shortest angle between the vector pointing to this field and the vector pointing to the farthest field
    %            Speculate_Major_axis=sqrt( (-cos(c)^2)/ ((sin(c)/r1)^2-R1^-2));
        Speculate_Major_axis=sqrt( (cos(c+pi/2)^2)/ (R1^-2-(sin(c+pi/2)/r1)^2));            
    % compress along the major axis
        if Speculate_Major_axis/r1<2 &Speculate_Major_axis^2>0&length(G_6peaks)>=6 % do not correct large distortion and % When peaks did not fit with ellipse
            phi_1=min_Peak(1,1)+pi/2; % 
            tilt_Autoc1=imrotate(Autoc{k},-phi_1*180/pi,'crop');% correct tilt
            T = [r1/Speculate_Major_axis    0  0; 0    1 0; 0    0  1];
            t_aff = affine2d(T);
            Cor_Autoc1 = imwarp(tilt_Autoc1,t_aff);  % compress along the major axis
        else
            Cor_Autoc1 = Autoc{k};
        end   



        [cor_y cor_x]=size(Cor_Autoc1); 
        [t_oy, t_ox] = ndgrid(-(cor_y-1)/2:1:(cor_y-1)/2,(cor_x-1)/2:-1:-(cor_x-1)/2); % make Cartesian coordinates, be careful
        [t_theta,t_rho] = cart2pol(t_ox,t_oy);
         B = imgaussfilt(Cor_Autoc1,3); % smoothed
        [cent, varargout]=FastPeakFind(B*255,0); %find peaks in autocorrecogram
         Angle=reshape(t_theta.*varargout,[],1);
        Dist=reshape(t_rho.*varargout,[],1);
        tmp=[Angle, Dist];
        G_all_peaks1=tmp(tmp(:,1)~=0&tmp(:,2)~=0,:);
        %G_all_peaks=tmp(find(tmp(:,2)>=10),:); % Remove peak around (<10cm) to center
        G_6peaks1=G_all_peaks1(knnsearch(G_all_peaks1(:,2),0,'K',6),:);
        in_r{k}=bin*min(G_6peaks1(:,2))-10;

    % trim inner and outer mask 
        if isempty(in_r{k})==1 % when only one specific peak
            in_r{k}=30
        end
        in = t_rho*bin > round(in_r{k});
        out = t_rho*bin < round(in_r{k})+20;
        inout = in.*out ;
        trim_Autoc{k} = Cor_Autoc1.*inout;
        trim_Autoc{k} = trim_Autoc{k}./inout;
        r60= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},60,'crop'));
        r120= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},120,'crop'));
        r30= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},30,'crop'));
        r90= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},90,'crop'));
        r150= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},150,'crop'));

        tilt_GS1 = min([r60, r120]) - max([r30, r90, r150]);

    % case2: major axis=farthest peak
        Speculate_Minor_axis=sqrt( (sin(c)^2)/ (1/r1^2-(cos(c)/R1)^2));


        % compress along the major axis
        if R1/Speculate_Minor_axis<2 &Speculate_Minor_axis^2>0 &length(G_6peaks)>=6% do not correct large distortion
            phi_2=max_Peak(1,1); % 
            tilt_Autoc2=imrotate(Autoc{k},-phi_2*180/pi,'crop');% correct tilt
            T = [Speculate_Minor_axis/R1    0  0; 0    1 0; 0    0  1];
            t_aff = affine2d(T);
            Cor_Autoc2 = imwarp(tilt_Autoc2,t_aff);  % compress along the major axis
        else
            Cor_Autoc2 = Autoc{k};
        end   

        [cor_y cor_x]=size(Cor_Autoc2); 
        [t_oy, t_ox] = ndgrid(-(cor_y-1)/2:1:(cor_y-1)/2,(cor_x-1)/2:-1:-(cor_x-1)/2); % make Cartesian coordinates, be careful
        [t_theta,t_rho] = cart2pol(t_ox,t_oy);
         B = imgaussfilt(Cor_Autoc2,3); % smoothed
        [cent, varargout]=FastPeakFind(B*255,0); %find peaks in autocorrecogram
         Angle=reshape(t_theta.*varargout,[],1);
        Dist=reshape(t_rho.*varargout,[],1);
        tmp=[Angle, Dist];
        G_all_peaks2=tmp(tmp(:,1)~=0&tmp(:,2)~=0,:);
        G_6peaks2=G_all_peaks2(knnsearch(G_all_peaks2(:,2),0,'K',6),:);
        in_r{k}=bin*min(G_6peaks2(:,2))-10;

    % trim inner and outer mask 
        if isempty(in_r{k})==1 % when only one specific peak
            in_r{k}=30
        end
        in = t_rho*bin > round(in_r{k});
        out = t_rho*bin < round(in_r{k})+20;
        inout = in.*out ;
        trim_Autoc{k} = Cor_Autoc2.*inout;
        trim_Autoc{k} = trim_Autoc{k}./inout;
        r60= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},60,'crop'));
        r120= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},120,'crop'));
        r30= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},30,'crop'));
        r90= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},90,'crop'));
        r150= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},150,'crop'));

        tilt_GS2 = min([r60, r120]) - max([r30, r90, r150]);
    % chose bigger GS as Grid score
        rGS(k,round(in_r{k})+20) = max(tilt_GS1,tilt_GS2);

        if  isempty(max(rGS(k,rGS(k,:)~=0)))==0
            Grid_Score(k,1) = max(rGS(k,rGS(k,:)~=0));
        end
end
 %       Grid_Score(k,1) = max(rGS(k,rGS(k,:)~=0));
        %    end
    end

    Rand_Grid_Score(:,rrtt)= Grid_Score;
 end
 
a= Rand_Grid_Score;
a(isnan(a))=[];
Gth = prctile(reshape(a,[],1),95);
% {
a= Rand_Grid_Score;
b=zeros(RandNum,2);
for i=1:1:RandNum
    aa=a(:,1:i);
    aa(aa==0)=[];
    b(i,1)=prctile(reshape(aa,[],1),95);
    b(i,2)=mean(a(:,1:i),'all');
end
    
plot(b(:,1))
hold on
plot(b(:,2))
legend('p95','mean')
 print(strcat("Rand_Gth.jpg"), '-djpeg', '-r0')
 close
 
%imagesc(Autoc{k})


%}
    


