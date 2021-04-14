function GU_ExM_JaneliaCluster_calcPunctaDensityVol_1ch_generic(ex, mx, my, mz, OL, VS, MT, ch, EF, MVV, PXY, PZ, ipath, opath)

ex = str2double(ex);
MT = str2double(MT);% minimum threshold
EF =  str2double(EF);% expansion factor
MVV =  str2double(MVV);% min voxel volume
PXY =  str2double(PXY);% 
PZ =  str2double(PZ);% 

% addpath(genpath('/groups/betzig/home/upadhyayulas/Software/GU_PrivateRepository'));
% addpath(genpath('/groups/betzig/home/upadhyayulas/Software/GU_Repository'));

mx = str2double(mx);
my = str2double(my);
mz = str2double(mz);
OL = str2double(OL); % overlap size
VS = str2double(VS); % volume size
s = repmat(VS, [1,3]);

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
            xmax(nn) = (s(1)*i-(OL*(i-1)));
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            ymax(nn) = (s(2)*j-(OL*(j-1)));
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            zmax(nn) = (s(3)*k-(OL*(k-1)));
            
            nn = nn + 1;
        end
    end
end

% cmd1 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/nrs/saalfeld/igor/illumination-correction/Sample1_C1/stitching/decon-warped-export/export.n5" --output "/nrs/saalfeld/igor/Rui/Wholebrainsynapse/ch1" --channel 1  --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];

% cmd0 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/nrs/saalfeld/igor/illumination-correction/Sample1_C1/stitching/decon-warped-export/export.n5" --output "/nrs/saalfeld/igor/Rui/Wholebrainsynapse/ch0"  --channel 0 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
cmd0 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "' ipath '" --output "' opath '" --channel ' ch ' --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
    

% run script
% ch0rt = '/nrs/saalfeld/igor/Rui/Wholebrainsynapse/ch0/';
fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
fn1 = [opath fn];


if ~exist(fn1, 'file')
    tic
    [~, commandOut2] = system(cmd0);
    toc
end

GU_ExM_calcPunctaDensityVol_1ch(fn1, ...
        'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [MVV],'ExpansionFactor', EF, 'MinThreshold', MT,...
        'PixelSize', PXY,'PixelSizeZ', PZ, 'OTSUMaxPer', 99);
    
% GU_ExM_calcSynapseDensityVol(fn2, fn1,'MinThreshold', [550,720],'Verbose', true,'OTSUGaussKernel', 0.5);

% delete extracted tiff
delete(char(fn1))

fprintf('Completed %d of %d \n', ex,numel(1:nn-1));