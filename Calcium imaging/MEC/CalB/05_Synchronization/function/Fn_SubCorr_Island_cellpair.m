function  [data1 data2 data3 data4] =Fn_SubCorr_Island_cellpair(DATA,ISLAND);
% data1: ALL-ALL island
% data2: Intra Island(I1-I1, I2-I2, I3-I3)
% data3: non Island
% data4:  trans Island(I1-I1, I2-I2, I3-I3)


% SigfCorr_ratio: ratio of sig. corr on each cell pair
% 1, initial cell;2-3, respond cell pair; 4, each corr
% NETWORK: cell pair matrix for each network
% 1:IntraIsland, 2:nonIsland, 3:InterIsland,
% clf
NumCell=size(DATA,1);


% All island
mIsland = find(ISLAND>0);
All_Island = nan(NumCell, NumCell, NumCell);
for q = 1:length(mIsland)
    for qq = 1:length(mIsland)
        i = mIsland(q);
        j = mIsland(qq);
        All_Island(:, i, j) = 1;
    end
end




% IntraIsland 
Label=unique(ISLAND);
IntraIsland=cell(length(Label)-1,1);
II=cell(length(Label)-1,1);
for l=1:1:length(Label)-1
    IntraIsland{l}=nan(NumCell,NumCell,NumCell);
    tmp=find(ISLAND==Label(l+1));
    for q=1:1:length(tmp)
        for qq=q+1:1:length(tmp)
            i=tmp(q);ii=tmp(qq);
            IntraIsland{l}(:,i,ii)=1;
        end
    end
    II{l}=find(ISLAND==l);
end

% nonIsland pair
nonIsland=nan(NumCell,NumCell,NumCell);
tmp=find(ISLAND==0);
for q=1:1:length(tmp)
    for qq=q+1:1:length(tmp)
        i=tmp(q);ii=tmp(qq);
        nonIsland(:,i,ii)=1;
    end
end

% trans Island pair
Label=unique(ISLAND);
tII=[];
if max(Label)>1
    xxx=1;
    for l=2:1:length(Label)
        for ll=l+1:1:length(Label)
           xxx=xxx+1;
        end
    end

    TransIsland=cell(xxx-1,1);
    tII=cell(xxx-1,1);
    xxx=1;
    for l=2:1:length(Label)
        for ll=l+1:1:length(Label)
        TransIsland{xxx}=nan(NumCell,NumCell,NumCell);
        tmp1=find(ISLAND==Label(l));
        tmp2=find(ISLAND==Label(ll));
    
        for q=1:1:length(tmp1)
            for qq=1:1:length(tmp2)
                i=tmp1(q);ii=tmp2(qq);
                TransIsland{xxx}(:,i,ii)=1;
            end
        end
        tII{xxx}=find(ISLAND==l|ISLAND==ll);
        xxx=xxx+1;

        end
    end
end



% 
% 

% data1: ALL-ALL island
% data2: Intra Island(I1-I1, I2-I2, I3-I3)
% data3: non Island
% data4:  trans Island(I1-I1, I2-I2, I3-I3)


data1 = reshape(DATA(mIsland, :, :) .* All_Island(mIsland, :, :) * 100,[],1);
data1 = data1(~isnan(data1));


data2=[];
for l=1:1:length(Label)-1
    tmp = reshape(DATA(II{l},:,:).*IntraIsland{l}(II{l},:,:)*100,[],1) ;
    data2 = [data2;tmp];
end
data2 = data2(~isnan(data2));

data3=[];
tmp=find(ISLAND==0);     
data3 = reshape(DATA(tmp,:,:).*nonIsland(tmp,:,:)*100,[],1);
data3 = data3(~isnan(data3));


data4=[];
for l=1:1:length(tII)
    tmp=reshape(DATA(tII{l},:,:).*TransIsland{l}(tII{l},:,:)*100,[],1);
    data4 = [data4;tmp];
end
data4 = data4(~isnan(data4));


 
   




end
