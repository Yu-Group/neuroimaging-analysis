% /groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing
nn = 10936;
for i = 1:nn-1
    % GU_ExM_JaneliaCluster_SampleScanProcessing171129YFPSample(ex)
    
    submit = ['bsub -n 16 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "SS' num2str(i) '" -o ~/logs/SS' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_SampleScanProcessing171129YFPSample ' num2str(i)  ''''];
    [stat, res] = system(submit,'-echo');
end

%% use this after deskewing
for i = 1:nn-1
    % GU_ExM_JaneliaCluster_SampleScanProcessing171129YFPSample(ex)
    
    submit = ['bsub -n 16 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "SS' num2str(i) '" -o ~/logs/SS' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_SampleScanDeconRot171129YFPSample_v2 ' num2str(i)  ''''];
    [stat, res] = system(submit,'-echo');
end

%%

% for DS stitched tiles
mx = 26629;
my = 10663;
mz = 5842;
OL = 50;
s = [750,750,750];

% % for raw tiles
% mx = 10440;
% my = 10663;
% mz = 5842;
% OL = 100;
% s = [500,500,500];

nx = ceil(mx/(s(1)-OL));
ny = ceil(my/(s(2)-OL));
nz = ceil(mz/(s(3)-OL));

theta = 32.2;
xyPixelSize = 0.104;
dz = 0.340;
dx = cosd(theta)*dz/xyPixelSize;


% generate x y z cropping coordinates
clear xmin xmax ymin zmin ymax zmax

xmin = zeros(nx*ny*nz,1);
xmax = zeros(nx*ny*nz,1);
ymin = zeros(nx*ny*nz,1);
ymax = zeros(nx*ny*nz,1);
zmin = zeros(nx*ny*nz,1);
zmax = zeros(nx*ny*nz,1);

xminDS = zeros(nx*ny*nz,1);
xmaxDS = zeros(nx*ny*nz,1);

nn = 1;
nxOut = ((s(3)-1)*dx);
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            if j <ny
                ymax(nn) = (s(2)*j-(OL*(j-1)));
            else
                ymax(nn) = my;
            end
            
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            if k < nz
                zmax(nn) = (s(3)*k-(OL*(k-1)));
            else
                zmax(nn) = mz;
            end
            
            %              nxOut = ((zmax(nn)-zmin(nn))*dx);
            xmin(nn) = (1+((i-1)*(s(1)-OL)));
            %
            %             xminDS(nn) = xmin(nn)-(k-1)*(nxOut - (nxOut*OL/s(3)));
            xminDS(nn) = xmin(nn)-(k-1)*(nxOut - (nxOut*OL/s(3))) +  nxOut * (1-(zmax(nn)-zmin(nn)+1)/s(3));
            if i <nx
                xmax(nn) = (s(1)*i-(OL*(i-1)));
            else
                xmax(nn) = mx;
            end
            %             xmaxDS(nn) = xmax(nn) + nxOut - (k-1)*(nxOut - (nxOut*OL/s(3)));
            xmaxDS(nn) = xmax(nn) + nxOut - (k-1)*(nxOut - (nxOut*OL/s(3))) +  nxOut * (1-(zmax(nn)-zmin(nn)+1)/s(3));
            nn = nn + 1;
        end
    end
end

%% figures

figure, scatter3(xminDS, ymin, zmin, '.')
hold on
scatter3(xmaxDS, ymax, zmax, '.')
hold on
scatter3(mean([xmaxDS, xminDS]'), mean([ymax,ymin]'), mean([zmax,zmin]'), '.')
%%
figure,
plot(nxOut * (1-(zmax-zmin+1)/s(3)))
%% find missing files
rt = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/matlab_decon_16bit/';
% rt = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/DSR_dx0.1229_dz0.16515/';
% rt = 'U:\igor\Rui\Wholebrainsynapse\ch1\Analysis\DataStructures\';
fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name};
CompletedJobsIdx = zeros(1, nn-1);
for ex = 1:nn-1
    fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
    CompletedJobsIdx(ex) = ismember(fn,fnlist);%exist(fn(2:end), 'file');
end

job2resubmit = find(CompletedJobsIdx==0);

for i = job2resubmit
    submit = ['bsub -n 16 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "SS' num2str(i) '" -o ~/logs/SS' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_SampleScanProcessing171129YFPSample ' num2str(i)  ''''];
    [stat, res] = system(submit,'-echo');
end

%% get histogram

rt = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/matlab_decon/MIPs/';
fn = dir([rt '*.tif']);
fn = {fn.name};
uBound = zeros(1, numel(fn));
parfor_progress(numel(fn));
parfor i = 1:numel(fn)
    im = readtiff([rt fn{i}]);
    uBound(i) = prctile(im(:), 99.9);
    parfor_progress;
end
parfor_progress(0);
%% rescale the data to 16 bit based ont he 99.9th percentile of the data
maxR{1} = 512658;
maxR{2} = 386767;%round(prctile(uBound(uBound>0),99.9)); % 512658 for ch0; 386767 for ch1
rt{1} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/DSR_dx0.1229_dz0.16515/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/DSR_dx0.1229_dz0.16515/';

srt{1} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/DSR_dx0.1229_dz0.16515_16bit/';
srt{2} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/DSR_dx0.1229_dz0.16515_16bit/';

rt2{1} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/matlab_decon/';
rt2{2} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/matlab_decon/';

srt2{1} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/matlab_decon_16bit/';
srt2{2} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/matlab_decon_16bit/';

a = 1;
for j = 1:2
    for ex = 1:nn-1 %job2resubmit%
        fn = [rt{j} num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
        fn2 = [rt2{j} num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
        
        sfn = [srt{j} num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
        sfn2 = [srt2{j} num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
        
        %     GU_scaleContrast_16bit(volpath, rangeInmin, rangeInmax)

        if exist(fn, 'file') %&& ~exist(sfn, 'file')
        submit = ['bsub -n 1 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "RS' num2str(a) '" -o ~/logs/RS' num2str(a) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(a) ' ; /groups/betzig/home/upadhyayulas/GU_scaleContrast_16bit ' fn ' ' num2str(0) ' ' num2str(maxR{j})   ''''];
        [~, ~] = system(submit,'-echo');
        a = a+1;
        end
        
        if exist(fn2, 'file') %&& ~exist(sfn2, 'file')
        submit = ['bsub -n 1 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "RS' num2str(a) '" -o ~/logs/RS' num2str(a) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(a) ' ; /groups/betzig/home/upadhyayulas/GU_scaleContrast_16bit ' fn2 ' ' num2str(0) ' ' num2str(maxR{j})   ''''];
        [~, ~] = system(submit,'-echo');
        a = a+1;
        end
    end
end
%% re-rotate data without cropping

%GU_ExM_JaneliaCluster_SampleScanRotation171129YFPSample

for i = 1:nn-1
    submit = ['bsub -n 16 -R"affinity[core(1)]" -R"select[broadwell]" -J' ' "SS' num2str(i) '" -o ~/logs/SS' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_JaneliaCluster_SampleScanRotation171129YFPSample ' num2str(i)  ''''];
    [stat, res] = system(submit,'-echo');
end

%% rename DSR 16bit files
clear rt
rt{1} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/DSR_dx0.1229_dz0.16515/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/DSR_dx0.1229_dz0.16515/';
% rt{1} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/matlab_decon_16bit/';
% rt{2} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/matlab_decon_16bit/';
mkdir([rt{1} 'allRenamedTiles' filesep])
mkdir([rt{2} 'allRenamedTiles' filesep])
theta = 32.2;
xyPixelSize = 0.104;
dz = 0.340;
zxRatio = sind(abs(theta))*dz/xyPixelSize;
dx = cosd(abs(theta))*dz/xyPixelSize;
mx = 26629;
my = 10663;
mz = 5842;
fileSufix = '';
for i = 1:numel(rt)
    GU_ExM_WriteRotatedCoordinates(rt{i}, mx, my, mz, zxRatio, theta, fileSufix)
end

%% 
clear rt
rt = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/OTFtest/';
fn = 'Scan_Iter_0002_CamB_ch1_CAM1_stack0000_560nm_0000000msec_0004819447msecAbs_004x_003y_000z_0002t';
tic
% im = readtiff([rt 'DS_crop.tif']);
im = readtiff([rt fn '.tif']);
toc
im = im(:,:, 1:500);
% deskew
theta = 32.2;
dz = 0.34;
xyPixelSize = 0.104;
tic;
volout = deskewFrame3D(im, theta, dz, xyPixelSize,'Crop', true,'reverse', true); toc
writetiff((volout), [rt fn '_ds32bit.tif']);

%% decon
psfrt = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/PSFs_objective/';
PSFpath = [psfrt '488PSF_OS_100nmbd_03.tif'];

xyPixelSize = 0.104;
dzPSF = 0.1;
nIter = 15;
zAniso = sind(32.2)*0.34/xyPixelSize;
Background = 100;

tic
volout = RLdecon_for_ExMPipeline(im, [rt fn '_ds.tif'], PSFpath, Background, nIter, dzPSF, zAniso, 0, [], 0, ...
            xyPixelSize, 0, 0, 0, 0, [0,0,1],0, []) ;toc
%         rsvo = uint16(scaleContrast(volout,[0, prctile(volout(:),99.9)], [0 2^16-1]));
%         writetiff(rsvo, [rt 'decon_crop.tif']);
        writetiff(volout, [rt 'decon_crop_32bit.tif']);
        mask = logical(im);
        se = strel('disk', 50);
            mask = imerode(mask,se);
            
            volout(~mask) = 0;
            clear mask
             writetiff(volout, [rt 'decon_crop_32bit_erroded.tif']);