function GU_ExM_JaneliaCluster_WholeBrainCalcSurfaceAreaVolume(ex)

ex = str2double(ex);
% addpath(genpath('/groups/betzig/home/upadhyayulas/Software/GU_PrivateRepository'));
% addpath(genpath('/groups/betzig/home/upadhyayulas/Software/GU_Repository'));

mx = 15055;
my = 27964;
mz = 6647;
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
            xmax(nn) = (s(1)*i-(OL*(i-1)));
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            ymax(nn) = (s(2)*j-(OL*(j-1)));
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            zmax(nn) = (s(3)*k-(OL*(k-1)));
            
            nn = nn + 1;
        end
    end
end

cmd0 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/export.n5" --output "/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/SA_Vol"  --channel 0 --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];


% run script
ch0rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/SA_Vol/';
if ~exist(ch0rt, 'dir')
    mkdir(ch0rt);
end

sch0rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/SA_Vol/dataStructures/';
if ~exist(sch0rt, 'dir')
    mkdir(sch0rt);
end

fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
fn1 = [ch0rt fn];

if ~exist(fn1, 'file')
    tic
    [~, commandOut2] = system(cmd0);
    toc
end

mask = logical(readtiff(fn1));

[vol, sarea, VolcoordXYZ, SAcoordXYZ, iVolcoordXYZ] = GU_calc_surfaceAreaVolumeFromCellMask_frame(mask);
save([sch0rt fn(1:end-3) 'mat'], 'vol', 'sarea', 'VolcoordXYZ', 'SAcoordXYZ', 'iVolcoordXYZ');

% delete extracted tiff
delete(char(fn1))
