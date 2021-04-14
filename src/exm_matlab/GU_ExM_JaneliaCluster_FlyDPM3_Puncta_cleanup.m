function GU_ExM_JaneliaCluster_FlyDPM3_Puncta_cleanup(ex)

ex = str2double(ex);

sigVox(1) = GU_calcGaussianIntegral(1, [2.5 2.5 2.5]);

mx = 10617;
my = 10712;
mz = 6638;
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
cmd0 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/groups/betzig/betziglab/4Stephan/160919_DPM3/Stitch_Igor/decon-export/export.n5" --output "/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch0" --channel 0 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
cmd1 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/groups/betzig/betziglab/4Stephan/160919_DPM3/Stitch_Igor/decon-export/export.n5" --output "/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch1" --channel 1 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];

%
% run script
vol1rt = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch1/';
vol2rt = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch0/';

fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
fnv1 = [vol1rt fn]; % nc82
fnv2 = [vol2rt fn]; % yfp

% yfp
if ~exist(fnv2, 'file')
    tic
    [~, commandOut2] = system(cmd0);
    toc
end

% mpb
if ~exist(fnv1, 'file')
    tic
    [~, commandOut3] = system(cmd1);
    toc
end

% depuncate and calculate the local maxima
for ni = 1:numel(sigVox)
    GU_ExM_calcPunctaDensityVol(fnv1, fnv2,...
        'FindLocalMaxima', true, 'MinThreshold', [1100,2500], 'Verbose', true,'OTSUGaussKernel', 0.5,'AssSphereSize', 7,...
        'MinVoxelVolume', [round(sigVox(ni)), 1008],'ExpansionFactor', 3.99); % sigma 3.9
end
% delete extracted tiff
delete(char(fnv1))
delete(char(fnv2))
