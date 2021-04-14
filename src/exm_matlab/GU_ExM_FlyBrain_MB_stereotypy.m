%% GU_ExM_FlyBrain_MB_stereotypy

%% calculate min otsu threhold
rt = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/subsamples';
cd(rt);
% data = GU_loadConditionData
%%
% T = zeros(1,numel(data));
parfor i = 1:numel(data)
    im = readtiff(data(i).framePaths{1});
    %     T(i) = thresholdOtsu(im(im>0 & im<prctile(im(:),99.5)));
    GU_ExM_calcPunctaDensityVol_1chVol(data(i).framePaths{1}, im, ...
        'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [192],'ExpansionFactor', 3.99, 'Threshold', T(i),...
        'PixelSize', 0.102,'PixelSizeZ', 0.18,'OTSUMaxPer', 99)
end


%%

% location of n5
rt{1} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F1_L/Stitch_Igor/Stitch_decon_blend/export.n5';
rt{2} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F1_R/Stitch_Igor/Stitch_decon_blend/export.n5';
rt{3} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F2_L/Stitch_Igor/Stitch_decon_blend/export.n5';
rt{4} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F3_L/Stitch_Igor/Stitch_decon_blend/export.n5';
rt{5} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F3_R/Stitch_Igor/Stitch_decon_blend/export.n5';
rt{6} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/RBPF1_R/Stitch_Igor/Stitch_decon_blend/export.n5_stage';

srt{1} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex1_F1L_z0p18/Tiles/';
srt{2} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex2_F1R_z0p18/Tiles/';
srt{3} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex3_F2L_z0p18/Tiles/';
srt{4} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex4_F3L_z0p18/Tiles/';
srt{5} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex5_F3R_z0p18/Tiles/';
srt{6} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex6_RBPF1R_z0p18/Tiles/';


ax = [4019, 4325, 4316, 4631, 4637, 5305];
ay = [3221, 3224, 4261, 3199, 3187, 3184];
az = [2754, 2808, 2787, 3216, 3286, 3772];
channel = [0, 0, 0, 0,0,1];
% T = [518.3020  705.7216  905.8000  803.0235  759.2078  979.1373]; %%99.9 prct
% T = [355.2471  471.9843  552.2941  496.3412  492.1647  667.8235]; % 99.5
T = [311.1529  398.2549  453.0000  405.6824  402.0471  550.1137]; % 99 prct
jid = 1;
for i = 1:numel(rt)
    ipath = rt{i};
    opath = srt{i};
    
    if ~exist(opath, 'dir')
        mkdir(opath)
    end
    
    mx = ax(i); %1
    my = ay(i); %7
    mz = az(i); % 41
    OL = 50;% overlap size
    VS = 500;% volume size
    MT = T(i);
    ch = channel(i);
    EF = '3.99';
    MVV = '192';
    PXY = '0.102';
    PZ = '0.18';
    
    s = repmat(VS, [1,3]);
    nx = ceil(mx/(s(1)-OL));
    ny = ceil(my/(s(2)-OL));
    nz = ceil(mz/(s(3)-OL));
    
    % generate x y z cropping coordinates
    clear xmin xmax ymin zmin ymax zmax nn ex
    
    xmin = zeros(nx*ny*nz,1);
    xmax = zeros(nx*ny*nz,1);
    ymin = zeros(nx*ny*nz,1);
    ymax = zeros(nx*ny*nz,1);
    zmin = zeros(nx*ny*nz,1);
    zmax = zeros(nx*ny*nz,1);
    
    nn = 1;
    for k = 1:nz
        for j = 1:ny
            for l = 1:nx
                xmin(nn) = (1+((l-1)*(s(1)-OL)));
                xmax(nn) = (s(1)*l-(OL*(l-1)));
                ymin(nn) = (1+((j-1)*(s(2)-OL)));
                ymax(nn) = (s(2)*j-(OL*(j-1)));
                zmin(nn) = (1+((k-1)*(s(3)-OL)));
                zmax(nn) = (s(3)*k-(OL*(k-1)));
                
                nn = nn + 1;
            end
        end
    end
    mx = num2str(mx);
    my = num2str(my);
    mz = num2str(mz);
    OL = num2str(OL);
    VS = num2str(VS);
    MT = num2str(MT);
    ch = num2str(ch);
    
    for ex = 1:nn-1
        %         GU_ExM_JaneliaCluster_calcPunctaDensityVol_1ch_generic(num2str(ex), mx, my, mz, OL, VS, MT, ch, EF, MVV, PXY, PZ, ipath, opath)
        
        submit = ['bsub -n 1 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "nc82' num2str(jid) '" -o ~/logs/nc82' ...
            num2str(jid) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  ...
            num2str(jid) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_calcPunctaDensityVol_1ch_generic ' ...
            num2str(ex) ' ' mx ' ' my ' ' mz ' ' OL ' ' VS ' ' MT ' ' ch ' ' EF ' ' MVV ' ' PXY ' ' PZ ' ' ipath ' ' opath ''''];
        [~, ~] = system(submit,'-echo');
        jid = jid+1;
    end
end


%% combine all the data - from the recalculated set

OL = 50;% overlap size
VS = 500;% volume size

drt{1} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex1_F1L_z0p18/Tiles/Analysis_1ch/192/DataStructures/';
drt{2} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex2_F1R_z0p18/Tiles/Analysis_1ch/192/DataStructures/';
drt{3} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex3_F2L_z0p18/Tiles/Analysis_1ch/192/DataStructures/';
drt{4} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex4_F3L_z0p18/Tiles/Analysis_1ch/192/DataStructures/';
drt{5} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex5_F3R_z0p18/Tiles/Analysis_1ch/192/DataStructures/';
drt{6} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex6_RBPF1R_z0p18/Tiles/Analysis_1ch/192/DataStructures/';

srt{1} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex1_F1L_z0p18/Tiles/';
srt{2} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex2_F1R_z0p18/Tiles/';
srt{3} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex3_F2L_z0p18/Tiles/';
srt{4} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex4_F3L_z0p18/Tiles/';
srt{5} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex5_F3R_z0p18/Tiles/';
srt{6} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex6_RBPF1R_z0p18/Tiles/';

ax = [4019, 4325, 4316, 4631, 4637, 5305];
ay = [3221, 3224, 4261, 3199, 3187, 3184];
az = [2754, 2808, 2787, 3216, 3286, 3772];
nds = zeros(1,numel(drt));
for k = 1:numel(drt)
    allData = [];
    
    mx = ax(k); %1
    my = ay(k); %7
    mz = az(k); % 41
    OL = 50;% overlap size
    VS = 500;% volume size
    MT = T(k);
    ch = channel(k);
    EF = '3.99';
    MVV = '192';
    PXY = '0.102';
    PZ = '0.18';
    
    s = repmat(VS, [1,3]);
    nx = ceil(mx/(s(1)-OL));
    ny = ceil(my/(s(2)-OL));
    nz = ceil(mz/(s(3)-OL));
    
    % generate x y z cropping coordinates
    clear xmin xmax ymin zmin ymax zmax nn ex
    
    xmin = zeros(nx*ny*nz,1);
    xmax = zeros(nx*ny*nz,1);
    ymin = zeros(nx*ny*nz,1);
    ymax = zeros(nx*ny*nz,1);
    zmin = zeros(nx*ny*nz,1);
    zmax = zeros(nx*ny*nz,1);
    
    nn = 1;
    for dd = 1:nz
        for j = 1:ny
            for l = 1:nx
                xmin(nn) = (1+((l-1)*(s(1)-OL)));
                xmax(nn) = (s(1)*l-(OL*(l-1)));
                ymin(nn) = (1+((j-1)*(s(2)-OL)));
                ymax(nn) = (s(2)*j-(OL*(j-1)));
                zmin(nn) = (1+((dd-1)*(s(3)-OL)));
                zmax(nn) = (s(3)*dd-(OL*(dd-1)));
                
                nn = nn + 1;
            end
        end
    end
    nds(k) = nn-1;
    
    % loop through and merge all data tiles
    parfor_progress(nds(k));
    for i = 1:nds(k)
        fn = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_192_densityData.mat'];
        
        variableInfo = who('-file', [drt{k} filesep fn]);
        PF = ismember('filteredData', variableInfo);
        
        if PF
            t = load([drt{k} filesep fn], 'filteredData');
            idx = t.filteredData(:,1)>OL/2 &  t.filteredData(:,2)>OL/2 & t.filteredData(:,3)>OL/2 & ...
                t.filteredData(:,1)<=(VS-OL) &  t.filteredData(:,2)<=(VS-OL) & t.filteredData(:,3)<=(VS-OL);
            
            t.filteredData = t.filteredData(idx,:);
            t.filteredData(:,1) = t.filteredData(:,1) + xmin(i)-1;
            t.filteredData(:,2) = t.filteredData(:,2) + ymin(i)-1;
            t.filteredData(:,3) = t.filteredData(:,3) + zmin(i)-1;
            allData = [allData; t.filteredData];
        end
        parfor_progress;
    end
    parfor_progress(0);
    save([srt{k} 'mergednc82FilteredDataNotZcorrected.mat'], 'allData');
end

%% read out masked positions

mrt{1} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F1_L/Stitch_Igor/Stitch_decon_blend/slice-tiff-binned/ch0_a3_mask.tif';
mrt{2} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F1_R/Stitch_Igor/Stitch_decon_blend/slice-tiff-binned/a3_mask.tif';
mrt{3} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F2_L/Stitch_Igor/Stitch_decon_blend/slice-tiff-binned/a3_mask.tif';
mrt{4} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F3_L/Stitch_Igor/Stitch_decon_blend/slice-tiff-binned/a3_mask.tif';
mrt{5} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/F3_R/Stitch_Igor/Stitch_decon_blend/slice-tiff-binned/a3_mask.tif';
mrt{6} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/RBPF1_R/Stitch_Igor/Stitch_decon_blend/slice-tiff-binned_stage/a3_mask.tif';

srt{1} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex1_F1L_z0p18/Tiles/';
srt{2} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex2_F1R_z0p18/Tiles/';
srt{3} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex3_F2L_z0p18/Tiles/';
srt{4} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex4_F3L_z0p18/Tiles/';
srt{5} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex5_F3R_z0p18/Tiles/';
srt{6} = '/groups/betzig/betziglab/4Stephan/181009_MBalpha3/GokulWorkspace/Ex6_RBPF1R_z0p18/Tiles/';
bf = [8,8,5]; % mask binning factors in x,y,z
MB_nc82 = zeros(1,numel(mrt));
MB_nc82_density = zeros(1,numel(mrt));
maskVol = zeros(1,numel(mrt));
ef = 3.99;
pxy = 0.102/ef;
pz = 0.180/ef;
zAniso = pz/pxy;

for k = 1:numel(mrt)
    clear t AD X Y Z Xq Yq Zq idx mask
    tic
    t = load([srt{k} 'mergednc82FilteredDataNotZcorrected.mat'], 'allData');
    mask = (readtiff(mrt{k}));
    maskVol(k) = nnz(mask(:)) * pxy^2 * pz * bf(1) * bf(2) * bf(3);
    % adjust nc82 localization locations to match binned mask volume
    AD(:,1) = round(t.allData(:,1)/bf(1));
    AD(:,2) = round(t.allData(:,2)/bf(2));
    AD(:,3) = round(t.allData(:,3)/bf(3));
    
    [sy, sx, sz] = size(mask);
    [X, Y, Z] = meshgrid(1:sx, 1:sy, 1:sz);
    Xq = AD(:,1);
    Yq = AD(:,2);
    Zq = AD(:,3);
    idx = interp3(X,Y,Z,mask,Xq,Yq,Zq,'nearest');
    MB_nc82(k) = nnz(idx);
    MB_nc82_density(k) = nnz(idx)/maskVol(k);
    
    toc
end

figure
bar(1:6, MB_nc82_density);
xticks(1:6);
xticklabels({'F1L', 'F1R','F2L' , 'F3L', 'F3R', 'RBPF1'})


figure
bar(1:6, maskVol);
xticks(1:6);
xticklabels({'F1L', 'F1R','F2L' , 'F3L', 'F3R', 'RBPF1'})