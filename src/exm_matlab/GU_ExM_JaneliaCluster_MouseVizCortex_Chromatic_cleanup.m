function GU_ExM_JaneliaCluster_MouseVizCortex_Chromatic_cleanup(ex)

ex = str2double(ex);

sigVox(1) = GU_calcGaussianIntegral(1, [3.2 3.2 3.2]);

mx = 10727;
my = 39519;
mz = 5413;
OL = 50;
s = repmat(750, [1,3]);

nx = ceil(mx/(s(1)-OL));
ny = ceil(my/(s(2)-OL));
nz = ceil(mz/(s(3)-OL));

% generate x y z cropping coordinates
clear xmin xmax ymin zmin ymax zmax

xmin = zeros(nx*ny*nz,1);
xmax = zeros(nx*ny*nz,1);
ymin = zeros(nx*ny*nz,1);
ymax = zeros(nx*ny*nz,1);
zmin = zeros(nx*ny*nz,1);
zmax = zeros(nx*ny*nz,1);

nn = 1;
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            xmin(nn) = (1+((i-1)*(s(1)-OL)));
            if i <nx
                xmax(nn) = (s(1)*i-(OL*(i-1)));
            else
                xmax(nn) = mx;
            end
            
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
            
            nn = nn + 1;
        end
    end
end

% rt = '/groups/betzig/betziglab/Ved/Data/20170613_ExM/VisualYFPmbpcaspr/Stitch_Igor/restitching-incremental/decon-export/';
cmd0 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/groups/betzig/betziglab/Ved/Data/20170613_ExM/VisualYFPmbpcaspr/Stitch_Igor/restitching-incremental/decon-export/export.n5" --output "/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch0" --channel 0 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
cmd1 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/groups/betzig/betziglab/Ved/Data/20170613_ExM/VisualYFPmbpcaspr/Stitch_Igor/restitching-incremental/decon-export/export.n5" --output "/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch1" --channel 1 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
cmd2 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/groups/betzig/betziglab/Ved/Data/20170613_ExM/VisualYFPmbpcaspr/Stitch_Igor/restitching-incremental/decon-export/export.n5" --output "/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch2" --channel 2 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
%
%
% run script
vol1rt = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch2/';
vol2rt = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch0/';
vol3rt = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/ch1/';

fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
fnv1 = [vol1rt fn]; % casper
fnv2 = [vol2rt fn]; % yfp
fnv3 = [vol3rt fn]; % mbp

% casper
if ~exist(fnv1,'file')
    tic
    [~, commandOut1] = system(cmd2);
    toc
    GU_ExM_zOffset3Dframe(fnv1, 1,6);
    %     GU_ExM_zOffset3Dframe(fnv1, 0,5);
end

% yfp
if ~exist(fnv2, 'file')
    tic
    [~, commandOut2] = system(cmd0);
    toc
    %     GU_ExM_zOffset3Dframe(fnv2, 1,6);
end

% mpb
if ~exist(fnv3, 'file')
    tic
    [~, commandOut3] = system(cmd1);
    toc
    GU_ExM_zOffset3Dframe(fnv3, 1,1);
end

% depuncate and calculate the local maxima
for ni = 1:numel(sigVox)
    GU_ExM_calcPunctaDensityVol(fnv1, fnv2,...
        'FindLocalMaxima', false, 'MinThreshold', [1400,1100], 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [round(sigVox(ni)), 1008],'ExpansionFactor', 3.62); % sigma 3.9
    % T = 1400 for casper; T = 1100 for yfp; T = 200 for myelin
    
    GU_ExM_calcPunctaDensityVol_1ch(fnv3,...
        'FindLocalMaxima', false, 'MinThreshold',200, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', 1008,'ExpansionFactor', 3.62)
end
% delete extracted tiff
delete(char(fnv1))
delete(char(fnv2))
delete(char(fnv3))