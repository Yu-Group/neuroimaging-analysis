function GU_ExM_JaneliaCluster_Mouse_ReCalcYFPChannel(ex)

ex = str2double(ex);

% addpath(genpath('/groups/betzig/home/upadhyayulas/Software/GU_PrivateRepository'));
% addpath(genpath('/groups/betzig/home/upadhyayulas/Software/GU_Repository'));
% rt_vol1 = '/nrs/saalfeld/igor/170217_SomatoYFP_HB/stitching/flip-x/confidence-interval/restitching/restitching-3rd-phase/n5-export-decon/export.n5/c2';
% rt_vol2 = '/nrs/saalfeld/igor/170217_SomatoYFP_HB/stitching/flip-x/confidence-interval/restitching/restitching-3rd-phase/n5-export-decon/export.n5/c0';
% 

sigVox(1) = 1008;
% sigVox(1) = 1008;
% sigVox(2) = 2500;
% sigVox(3) = 5000;
% sigVox(4) = 7500;

mx = 10580;
my = 66496;
mz = 4488;
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

% cmd1 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/nrs/saalfeld/igor/170217_SomatoYFP_HB/stitching/flip-x/confidence-interval/restitching/restitching-3rd-phase/n5-export-decon/export.n5" --output "/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2" --channel 2  --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
cmd0 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/nrs/saalfeld/igor/170217_SomatoYFP_HB/stitching/flip-x/confidence-interval/restitching/restitching-3rd-phase/n5-export-decon/export.n5" --output "/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0"  --channel 0 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];


% run script
vol2rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/';
fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
% vol1rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/';
fnv2 = [vol2rt fn];
% fnv1 = [vol1rt fn];

% if ~exist(fnv1,'file')
%     tic
%     [~, commandOut1] = system(cmd1);
%     toc
% end

if ~exist(fnv2, 'file')
    tic
    [~, commandOut2] = system(cmd0);
    toc
end

for ni = 1:numel(sigVox)
    GU_ExM_calcPunctaDensityVol_1ch(fnv2,...
        'FindLocalMaxima', false, 'MinThreshold', [1900],'Threshold', [1900],'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [sigVox(ni)],'ExpansionFactor', 3.62);
end
% delete extracted tiff
delete(char(fnv2))
% delete (char(fnv1))
