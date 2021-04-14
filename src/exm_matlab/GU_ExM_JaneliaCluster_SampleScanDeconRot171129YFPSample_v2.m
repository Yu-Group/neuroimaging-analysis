function GU_ExM_JaneliaCluster_SampleScanDeconRot171129YFPSample_v2(ex)

% Dec 27, 2017; Gokul U
if ischar(ex)
    ex = str2double(ex);
end
% GU_ExM_SampleStackDeskewRotateWorkspace.m
%
rt = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/DSw/allRenamedTiles/export.n5/';
rt2 = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/DSw/allRenamedTiles/export.n5/';
psfrt = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/PSFs_objective/';
savertch0 = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/';
savertch1 = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/';

mx = 26629;
my = 10663;
mz = 5842;
OL = 50;
s = [750,750,750];

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
%
psffn{1} = [psfrt '488PSF_OS_100nmbd_03.tif'];
psffn{2} = [psfrt '560PSF_OS_100nmbd_02.tif'];
dz = 0.340; % in nm
opts = {'Crop', 1,...
    'xyPixelSize', 0.104,...
    'Deskew', 0,...
    'SkewAngle', 32.2,...
    'Rotate', 1,...
    'Reversed', 1,...
    'Deconvovle', 1,...
    'dzPSF', 0.1,...
    'Background', 100,...
    'nIter', 15,...
    'Save16bit', 0,...
    'EdgeArtifacts', 0};


cmd{1} = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "' rt '" --output "' savertch0 '" --channel 0  --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
cmd{2} = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "' rt2 '" --output "' savertch1 '" --channel 0 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];

fnv{1} = [savertch0 fn];
fnv{2} = [savertch1 fn];

for j = 1:numel(cmd)
    if ~exist(fnv{j},'file')
        tic
        [~, ~] = system(cmd{j});
        toc
    end
    vol = readtiff(fnv{j});
    SV(j) = sum(vol(:));
end

Proceed = sum(SV(:));

for j = 1:numel(cmd)
    if Proceed
        GU_ExM_SampleStackDeskewRotate(fnv{j}, dz, opts{:},'PSFpath', psffn{j});
    end
    delete(char(fnv{j}))
end
%%