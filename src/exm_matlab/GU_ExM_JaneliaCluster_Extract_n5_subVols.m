function GU_ExM_JaneliaCluster_Extract_n5_subVols(ipath, opath, mx,my,mz,csx,csy,csz, OL,ch, subPanelIdx)

% function generates sub volumes from n5 dataset
% Gokul Upadhyayula, Jan 2018

if ischar(mx)
    mx = str2double(mx);
end

if ischar(my)
    my = str2double(my);
end

if ischar(mz)
    mz = str2double(mz);
end

if ischar(csx)
    csx = str2double(csx);
end

if ischar(csy)
    csy = str2double(csy);
end

if ischar(csz)
    csz = str2double(csz);
end

if ischar(OL)
    OL = str2double(OL);
end

if ischar(ch)
    ch = str2double(ch);
end

if ischar(subPanelIdx)
    subPanelIdx = str2double(subPanelIdx);
end
s = [csx,csy,csz];

nx = ceil(mx/(s(1)-OL));
ny = ceil(my/(s(2)-OL));
nz = ceil(mz/(s(3)-OL));

% generate x y z cropping coordinates
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

if subPanelIdx > 0
    ri = subPanelIdx;
else
    ri = 1:nn-1;
end

for ex = ri
    cmd1 = ['python /groups/betzig/betziglab/Tool_Codes/stitching/crop-n5.py --input "' ipath '" --output "' opath '" --channel ' num2str(ch) ' --min ' num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) ' --max ' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex))];
    
    % run script
    fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];
    ch1rt = opath;
    fn1 = [ch1rt fn];
    
    if ~exist(fn1,'file')
        tic
        [~, commandOut1] = system(cmd1);
        toc
    end
end