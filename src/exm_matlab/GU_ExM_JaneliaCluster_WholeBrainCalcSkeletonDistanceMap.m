function GU_ExM_JaneliaCluster_WholeBrainCalcSkeletonDistanceMap(ex)
if isstr(ex)
    ex = str2double(ex);
end
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
% run script
ch0rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/SA_Vol/dataStructures/';
fn = [ num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.mat'];
load([ch0rt fn], 'VolcoordXYZ');



sch0rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/SA_Vol/dataStructures/skel/';
if ~exist(sch0rt, 'dir')
    mkdir(sch0rt);
end

if ~isempty(VolcoordXYZ)
    volidx = sub2ind([750,750,750], VolcoordXYZ(:,2),VolcoordXYZ(:,1),VolcoordXYZ(:,3));
    vol = zeros([750,750,750], 'logical');
    vol(volidx) = 1;
    [cleanVol,~, ~, SkelcoordXYZ, SkelRadius] = GU_calcSkelDistMap(vol, 'Fill', true, 'minConnVox', 5000,'dierr', 7,'zAniso', 1.8557,'RemoveEdgeVoxels', OL);
    save([sch0rt fn(1:end-3) 'mat'], 'cleanVol', 'SkelcoordXYZ', 'SkelRadius','s','OL');
else
    cleanVol = [];
    SkelcoordXYZ = [];
    SkelRadius = [];
    s = [];
    OL = [];
    save([sch0rt fn(1:end-3) 'mat'], 'cleanVol', 'SkelcoordXYZ', 'SkelRadius','s','OL');
end
