function GU_ExM_JaneliaCluster_CorrectDANDensityCalc(ex)

ex = str2double(ex);


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


ch0rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0';
fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
ch1rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1';


matfile = [ch1rt filesep 'Analysis' filesep 'DataStructures' filesep fn '_densityData.mat'];
cleanDanpath = [ch0rt filesep 'Analysis' filesep 'clean' filesep fn '_clean.tif'];

GU_ExM_correctDANSynapseDensityCalc(matfile, cleanDanpath);