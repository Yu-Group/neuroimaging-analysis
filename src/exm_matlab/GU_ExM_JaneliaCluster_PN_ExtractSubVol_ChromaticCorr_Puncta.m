function GU_ExM_JaneliaCluster_PN_ExtractSubVol_ChromaticCorr_Puncta(ex)

ex = str2double(ex);

% addpath(genpath('/groups/betzig/home/upadhyayulas/Software/GU_PrivateRepository'));
% addpath(genpath('/groups/betzig/home/upadhyayulas/Software/GU_Repository'));
% rt_vol1 = '/nrs/saalfeld/igor/170217_SomatoYFP_HB/stitching/flip-x/confidence-interval/restitching/restitching-3rd-phase/n5-export-decon/export.n5/c2';
% rt_vol2 = '/nrs/saalfeld/igor/170217_SomatoYFP_HB/stitching/flip-x/confidence-interval/restitching/restitching-3rd-phase/n5-export-decon/export.n5/c0';

sigVox(1) = GU_calcGaussianIntegral(1, [2.5,2.5,2.5]);
% sigVox(2) = GU_calcGaussianIntegral(1, [2.9,2.9,2.9]);
% sigVox(3) = GU_calcGaussianIntegral(1, [3.2 3.2 3.2]);
% sigVox(4) = GU_calcGaussianIntegral(1, [3.6 3.6 3.6]);


mx = 9664;
my = 6721;
mz = 5186;
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

cmd0 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/nrs/saalfeld/igor/PN3new2-test/stitch-slabs/combined-manually-adjusted/restitching-test/iter0/decon-export/export.n5" --output "/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/ch0" --channel 0 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
cmd1 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/nrs/saalfeld/igor/PN3new2-test/stitch-slabs/combined-manually-adjusted/restitching-test/iter0/decon-export/export.n5" --output "/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/ch1" --channel 1 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
cmd2 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/nrs/saalfeld/igor/PN3new2-test/stitch-slabs/combined-manually-adjusted/restitching-test/iter0/decon-export/export.n5" --output "/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/ch2" --channel 2 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];


% run script
vol1rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/ch1/';
vol2rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/ch2/';
vol3rt = '/groups/betzig/betziglab/4Stephan/171120_PNnc82Syd1/ch0/';

fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
fnv1 = [vol1rt fn]; % nc82
fnv2 = [vol2rt fn]; % syd1
fnv3 = [vol3rt fn]; % PN

% nc82
if ~exist(fnv1,'file')
    tic
    [~, commandOut1] = system(cmd1);
    toc
    GU_ExM_zOffset3Dframe(fnv1, 1,1);
%     GU_ExM_zOffset3Dframe(fnv1, 0,5);
end

if ~exist(fnv2, 'file')
    tic
    [~, commandOut2] = system(cmd2);
    toc
    GU_ExM_zOffset3Dframe(fnv2, 1,6);
end

if ~exist(fnv3, 'file')
    tic
    [~, commandOut3] = system(cmd0);
    toc
%     GU_ExM_zOffset3Dframe(fnv3, 0,6);
end


% depuncate and calculate the local maxima
for ni = 1:numel(sigVox)
    % nc82 + PN
    GU_ExM_calcPunctaDensityVol(fnv1, fnv3,...
        'FindLocalMaxima', true, 'MinThreshold', [576,756],'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [round(sigVox(ni)), 1008],'ExpansionFactor', 3.99);
    % Syd1 + PN
    GU_ExM_calcPunctaDensityVol(fnv2, fnv3,...
        'FindLocalMaxima', true, 'MinThreshold', [119,756],'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [round(sigVox(ni)), 1008],'ExpansionFactor', 3.99);
end
% delete extracted tiff
% delete(char(fnv1))
% delete(char(fnv2))
% delete(char(fnv3))