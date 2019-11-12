%% Compare the functional connectivity and the structural connectivity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unzip the structural connection data
struct_conn_path = '.\'; 
pa1='weights.csv.gz';
filenames = gunzip(pa1,struct_conn_path);
% weights= csvread([struct_conn_path,'\weights.csv']); 
pa1 = 'nodes.csv.gz';
filenames = gunzip(pa1,struct_conn_path);
% nodes = csvread([struct_conn_path,'\nodes.csv']);

% Input your save path
savepath = 'XXX\'; 

% Input your cluster number
Cluster_num = 12;

% Note: A lot of memory needs to be consumed for loading the connection data. We read the data in batches.
% If your computer memory less than 16GB, set batch_size less than 500
batch_size = 500;

% Named format: 1)left hemisphere:'L_XXX'; 2)right hemisphere:'R_XXX';
struct_name = cell(Cluster_num,1);
struct_name{1}= 'MS'; struct_name{2}= 'L_BF'; struct_name{3}= 'L_FL'; struct_name{4}= 'L_PtA'; struct_name{5}= 'L_V'; struct_name{6}= 'RSP';
struct_name{7}= 'M2'; struct_name{8}= 'HL'; struct_name{9}= 'R_FL'; struct_name{10}= 'R_PtA'; struct_name{11}= 'R_V'; struct_name{12}= 'R_BF';

% Registration ratio
% The size of the mouse's brain differs from that of Allen mouse brain, so you need to scale Allen mouse brain.
h_scale = 1.4386;  % From left to right
v_scale = 1.5606;  % From anterior to posterior

% Input the translation between the Allen mouse template and your mouse brain diagram.
h_t = 18;
v_t = 35;
%% Run
load('mean_cluster_total.mat')
load('cluster_group.mat')
li=zeros(size(mb));
li(mb~=0) = cluster_index;
%figure;imagesc(li); colormap jet 
se = strel('disk', 1);
li3 = zeros(size(li));
for i=1:max(cluster_index)
    li2 = li;  li2(li2~=i)=0;
    li2 = imerode(li2,se);
    li3 = li3+li2;
end
%figure; imagesc(li3); colormap jet
targetmask_ind=csvread('targetmask.csv'); 
targetmask_ind = reshape(targetmask_ind,448962,3);
targetmask_ind = targetmask_ind(:,[2,1,3]);
[Isocotexmask, ~] = nrrdread('Isocotexmb_100um.nrrd'); % load the Isocortex Mask
% The cortex obtained by optical imaging in our paper does not contain PL, ILA, and ORB,
% and these areas overlapped with the MO area in the top view. So we need to get rid of the PL, ILA, and ORB regions.
[PL_ILA_ORBmask, ~] = nrrdread('PL_ILA_ORB_100um.nrrd');
Isocotexmask =Isocotexmask-PL_ILA_ORBmask; 

pro_result_bilateral = cell(Cluster_num,1);
weights= csvread([struct_conn_path,'\weights.csv']); 
nodes = csvread([struct_conn_path,'\nodes.csv']);
nodes = reshape(nodes,448962,size(weights,2))'; 
sourcemask_ind=csvread('sourcemask.csv');
% (132,80,114)
sourcemask_ind=sourcemask_ind'; 
sourcemask_ind = sourcemask_ind(:,[2,1,3]);

for clu_ind=1:length(struct_name)
    str_split = split(struct_name{clu_ind},'_');    
    pro_result_total_R = [];
    pro_result_total_L = [];   
    m_com = 0;
    m_flip = 0;
    m_right = 1;
    count_num = 1;
    if length(str_split)==1
        %  Merge left and right hemispheres
        m_com = 1;
        count_num = [count_num,2];
    else if str_split{1}=='L'
            % flip to right hemisphere
            m_flip = 1;
        end
    end
    for i=count_num
        if m_com && (i==1)
            m_right = 0;
            m_flip = 1;
        end
        if m_com && (i==2)
            m_right = 1;
            m_flip = 0;
        end
        li4 = li3;
        if m_com && m_right
            li4(:,1:round(1+size(mb,2)/2)) = 0;
        end
        if m_com && ~m_right
            li4(:,round(1+size(mb,2)/2):end) = 0;
        end
        if m_flip
            li4 = li4(:,end:-1:1);
        end
        pic = zeros(round(v_scale*132),round(h_scale*114));
        pic(v_t:v_t+size(mb,1)-1, h_t:h_t+size(mb,2)-1) = li4;
        li4 = single(imresize(pic,[132,114],'nearest'));
        %figure; imagesc(li4); colormap jet
        pro_mean = [];
        pro_result_total = [];
        [clu_ind, i]
        [A_P,L_R] = find(li4==clu_ind);
        coord_reg = [];
        for j=1:length(A_P)
            S_I = find(Isocotexmask(:,A_P(j),L_R(j))==1);
            coord_temp = [S_I, A_P(j)*ones(size(S_I)), L_R(j)*ones(size(S_I))];
            coord_reg = [coord_reg; coord_temp];
        end
        sour_ind = [];
        for j=1:size(coord_reg,1)
            sour_ind = [sour_ind; find(sum((sourcemask_ind-coord_reg(j,:)).^2,2)==0)];
        end        
        startp = 1; endp=min([batch_size,length(sour_ind)]);
        projection = zeros(1,size(nodes,2));
        for j=1:ceil(length(sour_ind)/batch_size)  
            (j-1)*batch_size+1
            projection = projection + sum(weights(sour_ind(startp:endp),:)*nodes,1);
            startp = endp+1;
            endp = min([endp+batch_size,length(sour_ind)]);
        end        
        pro_result = zeros(size(Isocotexmask));
        for j=1:size(targetmask_ind,1)
            indtt = targetmask_ind(j,:);
            pro_result(indtt(1),indtt(2),indtt(3)) = projection(j);
        end
        pro_result = pro_result.*Isocotexmask;
        if m_flip
            pro_result = pro_result(:,:,end:-1:1);
        end
        pro_result_total = pro_result;
        h_sum = squeeze(sum(pro_result));
        h_sum = (h_sum-min(h_sum(:)))/(max(h_sum(:))-min(h_sum(:)));
        h_sum2 = imresize(h_sum,[round(v_scale*132),round(h_scale*114)],'nearest');
        h_sum2 = h_sum2(v_t:v_t+size(mb,1)-1, h_t:h_t+size(mb,2)-1);
        h_sum2 = h_sum2.*double(mb);
        %figure; imagesc(h_sum2); colormap hot
        if m_right
            pro_result_total_R = pro_result_total;
        else
            pro_result_total_L = pro_result_total;
        end
    end
    
    if m_com % Merge left and right hemispheres
        pro_result_bilateral{clu_ind} = pro_result_total_L/max(pro_result_total_L(:))...
            + pro_result_total_R/max(pro_result_total_R(:));
    else
        pro_result_bilateral{clu_ind} = pro_result_total/max(pro_result_total(:));
    end
    pro_result_bilateral{clu_ind} = pro_result_bilateral{clu_ind}/max(pro_result_bilateral{clu_ind}(:));
end

%% Show the projection density maps
for clu_ind= 1:length(pro_result_bilateral)
    close all
    pro_result = pro_result_bilateral{clu_ind};
    h_sum = squeeze(sum(pro_result));
    h_sum = (h_sum-min(h_sum(:)))/(max(h_sum(:))-min(h_sum(:)));
    h_sum2 = imresize(h_sum,[round(v_scale*132),round(h_scale*114)],'nearest');
    h_sum2 = h_sum2(v_t:v_t+size(mb,1)-1, h_t:h_t+size(mb,2)-1);
    h_sum2 = h_sum2.*double(mb);
    figure; imagesc(h_sum2); colormap hot
    
    set(gcf, 'InvertHardCopy', 'off'); 
    set(gcf,'position',[50,50,350,350])
    set(gca,'xtick',[],'ytick',[]);
    RemovePlotWhiteArea(gca);
    pause(1)
    fig = gcf;
    fig.PaperPositionMode = 'auto'; 
    saveas(gcf,[savepath,struct_name{clu_ind},'-proden.tif'])
    saveas(gcf,[savepath,struct_name{clu_ind},'-proden.fig'])
    pause(2)
end

%% Calculate the spatial correlation coefficient
rr = [];
for clu_ind= 1:Cluster_num
    pro_result = pro_result_bilateral{clu_ind};
    h_sum = squeeze(sum(pro_result));
    h_sum = (h_sum-min(h_sum(:)))/(max(h_sum(:))-min(h_sum(:)));
    h_sum2 = imresize(h_sum,[round(v_scale*132),round(h_scale*114)],'nearest');
    h_sum2 = h_sum2(v_t:v_t+size(mb,1)-1, h_t:h_t+size(mb,2)-1);
    h_sum2 = h_sum2.*double(mb);
    %figure; imagesc(h_sum2); colormap hot
    
    li=zeros(size(mb));
    li(mb~=0)=mean_cluster_total_temp(clu_ind,:);
    li = imresize(li,[size(h_sum2,1),size(h_sum2,2)],'nearest');
    li(li<0) = 0;
    % figure; imagesc(li); colormap jet
    
    rr = [rr; corr(single(li(:)), h_sum2(:))];
end
rr=rr';

save([savepath,'PD_VS_FC_reg.mat'],'rr','pro_result_bilateral','mb')

