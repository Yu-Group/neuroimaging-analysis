function GU_ExM_CombineMouseYFPHomerCHS(yfpRT, yfpSuffix, HyRT, HySuffix, allHRT, allHRTSuffix, saveFolder, ex)

% this function combines volumes in different regions of the 16bit volume
% this assumes that the volumes are spatially non-overlapping and max projected
% Gokul Upadhyayula, Oct 2017
ex = str2double(ex);
RangeIn = {[0 15000], [0 2^16-1], [0 2^16-1]};
RangeOut = {[0 32767], [32768 49151], [49152 2^16-1]};

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

%yfp FN
fn{1} = [yfpRT filesep num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) yfpSuffix '.tif'];
%yfp associated Homer FN
fn{2} = [HyRT filesep num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) HySuffix '.tif'];
%all Homer FN
fn{3} = [allHRT filesep num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) allHRTSuffix '.tif'];

SavePath = [saveFolder filesep num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '_SpectralMixed.tif'];
tic
GU_CombineVolumesIntoSingleBitDepth(fn, RangeIn, RangeOut, SavePath);
toc

