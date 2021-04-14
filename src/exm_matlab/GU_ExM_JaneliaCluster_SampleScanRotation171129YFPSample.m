function GU_ExM_JaneliaCluster_SampleScanRotation171129YFPSample(ex)

if ischar(ex)
    ex = str2double(ex);
end
% GU_ExM_SampleStackDeskewRotateWorkspace.m
%

rt{1} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch0/';
rt{2} = '/groups/betzig/betziglab/4Stephan/171129_YFPsample4/Processing/ch1/';
mx = 26629;
my = 10663;
mz = 5842;
OL = 50;
s = [750,750,750];
%
%
% mx = 10440;
% my = 10663;
% mz = 5842;
% OL = 100;
% s = repmat(500,500,500);

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

dz = 0.340; % in nm
xyPixelSize = 0.104;
theta = 32.2;
zAniso = dz*sind(theta)/xyPixelSize;

fn = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '.tif'];


zf = cotd(abs(theta));
yf = 1;
xf = cosd(abs(theta)) + tand(abs(theta))*sind(abs(theta));
% rangeIn{1} = single([0 512658]);
% rangeIn{2} = single([0 386767]);

for j = 1:numel(rt)
    fnv = [rt{j} 'matlab_decon_16bit' filesep fn];
    dsrpath = [rt{j} filesep 'DSR_dx' num2str(xyPixelSize*xf) '_dz' num2str(xyPixelSize*zf)];
    if ~exist(dsrpath, 'dir')
        mkdir(dsrpath);
    end
    if exist(fnv, 'file')
        volout = readtiff(fnv);
        fprintf('Rotating Data...')
        tic
        volout = rotateFrame3D(volout, theta, zAniso, 'reverse', true,...
            'Crop', false, 'ObjectiveScan', false); toc
        
        fprintf('resampling Rotated Data...')
        tic
        volout = GU_resampleStack3D(volout, xf, yf, zf,'Interp', 'linear');toc
        
%         fprintf('rescaling to 16 bit...')
        %     tic
        %     volout = uint16(scaleContrast(volout, rangeIn{j}, [0 65535])); toc
        
        fprintf('writing rotated volume...')
        tic
        writetiff(uint16(volout), [dsrpath filesep fn]); toc
    end
end
%%