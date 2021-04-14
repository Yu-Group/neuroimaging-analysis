%% rotated size
mx = 10580;
my = 66496;
mz = 4488;
theta = 32;
zxRatio = .18/.097;
outSize = round([my mx*cos(theta)+mz*zxRatio*sin(abs(theta)) mz*zxRatio*cos(theta)+mx*sin(abs(theta))]);
hOS = outSize/2;
%
% %%
% mx = 750;
% my = 750;
% mz = 750;
% theta = 32;
% zxRatio = .18/.097;
% outSize = round([my mx*cos(theta)+mz*zxRatio*sin(abs(theta)) mz*zxRatio*cos(theta)+mx*sin(abs(theta))]);
% %%
clear cc
mx = 10580;
my = 66496;
mz = 4488*zxRatio;
OL = 50;
s = [750, 750, 750*zxRatio];

nx = ceil(mx/(s(1)-OL));
ny = ceil(my/(s(2)-OL));
nz = ceil(mz/(s(3)-OL*zxRatio));

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
            zmin(nn) = (1+((k-1)*(s(3)-OL*zxRatio)));
            if k < nz
                zmax(nn) = (s(3)*k-(OL*zxRatio*(k-1)));
            else
                zmax(nn) = mz;
            end
            nn = nn + 1;
        end
    end
end
theta = 32;

rt = 'X:\4Stephan\171104_Mousebrainsynapse\ch0\Analysis\1008\Rotated_indTiles\';
for ex = 1:nn-1
    cc(ex,1:3) = [(xmax(ex)-(xmax(ex)-xmin(ex)+1)/2), (ymax(ex)-(ymax(ex)-ymin(ex)+1)/2), (zmax(ex)-(zmax(ex)-zmin(ex)+1)/2)]; % current center x,y,z
end

occ = cc;

mx = 10580;
my = 66496;
mz = 4488;
OL = 50;
s = [750, 750, 750];

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
%%
% clear cc nc
% zf = cotd(32);
% yf = 1;
% xf = cosd(32) + tand(32)*sind(32);
% 
% 
% mx = 10580/xf;
% my = 66496;
% mz = 4488*zxRatio/zf;
% OL = 50;
% s = [750/xf, 750, 750*zxRatio/zf];
% 
% cx = mx/2;
% cy = my/2;
% cz = mz/2;
% 
% nx = ceil(mx/(s(1)-OL/xf));
% ny = ceil(my/(s(2)-OL));
% nz = ceil(mz/(s(3)-OL*zxRatio/zf));
% 
% % generate x y z cropping coordinates
% % clear xmin xmax ymin zmin ymax zmax
% 
% xminR = zeros(nx*ny*nz,1);
% xmaxR = zeros(nx*ny*nz,1);
% yminR = zeros(nx*ny*nz,1);
% ymaxR = zeros(nx*ny*nz,1);
% zminR = zeros(nx*ny*nz,1);
% zmaxR = zeros(nx*ny*nz,1);
% 
% nn = 1;
% for k = 1:nz
%     for j = 1:ny
%         for i = 1:nx
%             xminR(nn) = (1+((i-1)*(s(1)-OL/xf)));
%             if i <nx
%                 xmaxR(nn) = (s(1)*i-(OL/xf*(i-1)));
%             else
%                 xmaxR(nn) = mx;
%             end
%             yminR(nn) = (1+((j-1)*(s(2)-OL)));
%             if j <ny
%                 ymaxR(nn) = (s(2)*j-(OL*(j-1)));
%             else
%                 ymaxR(nn) = my;
%             end
%             zminR(nn) = (1+((k-1)*(s(3)-OL*zxRatio/zf)));
%             if k < nz
%                 zmaxR(nn) = (s(3)*k-(OL*zxRatio/zf*(k-1)));
%             else
%                 zmaxR(nn) = mz;
%             end
%             nn = nn + 1;
%         end
%     end
% end
% theta = 32;
% 
% for ex = 1:nn-1
%     cc(ex,1:3) = [(xmaxR(ex)-(xmaxR(ex)-xminR(ex)+1)/2), (ymaxR(ex)-(ymaxR(ex)-yminR(ex)+1)/2), (zmaxR(ex)-(zmaxR(ex)-zminR(ex)+1)/2)]; % current center x,y,z
% end
% % close all
% % figure, scatter3(cc(:,1),cc(:,2),cc(:,3));
% 
% [nc, R, t] = AxelRot(cc', theta*-1, [0,1,0], [cx+1,cy+1,cz+1]);
% nc = nc';
% close all; figure, scatter3(nc(:,1),nc(:,2),nc(:,3)); hold on; scatter3(occ(:,1),occ(:,2),occ(:,3));

%%
clear nc
zf = cotd(32);
yf = 1;
xf = cosd(32) + tand(32)*sind(32);


mx = 10580;
my = 66496;
mz = 4488*zxRatio;
OL = 50;

outSize = round([my mx*cos(theta)+mz*zxRatio*sin(abs(theta)) mz*zxRatio*cos(theta)+mx*sin(abs(theta))]);
% cx = outSize(2)/2;
% cy = outSize(1)/2;
% cz = outSize(3)/2;
cx = mx/2;
cy = my/2;
cz = mz/2;

theta = 32;

% close all
% figure, scatter3(cc(:,1),cc(:,2),cc(:,3));

[nc, R, t] = AxelRot(cc', theta*-1, [0,1,0], [cx+1,cy+1,cz+1]);
nc = nc';

% nc(:,1) = nc(:,1)/xf;
% nc(:,3) = nc(:,3)/zf;
close all; figure, scatter3(nc(:,1),nc(:,2),nc(:,3)); hold on; scatter3(occ(:,1),occ(:,2),occ(:,3)); hold on
nc(:,1) = nc(:,1)/xf;
nc(:,3) = nc(:,3)/zf;
scatter3(nc(:,1),nc(:,2),nc(:,3));
%% rename files based on the new center point

rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/RotatedStacks/';
rtcopy = [rt 'nonEmpty/'];
mkdir(rtcopy);
rtmat = [rt 'params/'];
hs = [750 750 750]/2; % half size
parfor_progress(nn-1);
FE = zeros(1, nn-1);

parfor i = 1:nn-1
    %      fnO = [num2str(occ(i,1)-hs(1)+1) ',' num2str(occ(i,2)-hs(2)+1) ',' num2str(occ(i,3)-hs(3)+1) '_' num2str(occ(i,1)+hs(1)) ',' num2str(occ(i,2)+hs(2)) ',' num2str(occ(i,3)+hs(3)) '_clean.tif'];
    fnO = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_SpectralMixed.tif'];
    if exist([rt fnO], 'file') && exist([rtmat fnO(1:end-3) 'mat'], 'file')
        Proceed = load([rtmat fnO(1:end-3) 'mat'], 'nVoxOccupied');
        if Proceed.nVoxOccupied > 0
            IMinfo = imfinfo([rt fnO]);
            hrsx = IMinfo(1).Width/2;
            hrsy = IMinfo(1).Height/2;
            hrsz = numel(IMinfo)/2;
            
            fnN = [num2str(nc(i,1)-hrsx+1) ',' num2str(nc(i,2)-hrsy+1) ',' num2str(nc(i,3)-hrsz+1) '_' num2str(nc(i,1)+hrsx) ',' num2str(nc(i,2)+hrsy) ',' num2str(nc(i,3)+hrsz) '_SpectralMixed.tif'];
            
            copyfile([rt fnO], [rtcopy fnN]);
            FE(i) = 1;
        end
    end
    parfor_progress;
end
parfor_progress(0);

%% rename rounded values
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/RotatedStacks/';
rtcopy = [rt 'nonEmptyRounded/'];
mkdir(rtcopy);
rtmat = [rt 'params/'];
hs = [750 750 750]/2; % half size
parfor_progress(nn-1);
FE = zeros(1, nn-1);

parfor i = 1:nn-1
    %      fnO = [num2str(occ(i,1)-hs(1)+1) ',' num2str(occ(i,2)-hs(2)+1) ',' num2str(occ(i,3)-hs(3)+1) '_' num2str(occ(i,1)+hs(1)) ',' num2str(occ(i,2)+hs(2)) ',' num2str(occ(i,3)+hs(3)) '_clean.tif'];
    fnO = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) '_SpectralMixed.tif'];
    if exist([rt fnO], 'file') && exist([rtmat fnO(1:end-3) 'mat'], 'file')
        Proceed = load([rtmat fnO(1:end-3) 'mat'], 'nVoxOccupied');
        if Proceed.nVoxOccupied > 0
            IMinfo = imfinfo([rt fnO]);
            hrsx = IMinfo(1).Width/2;
            hrsy = IMinfo(1).Height/2;
            hrsz = numel(IMinfo)/2;
            
            fnN = [num2str(round(nc(i,1)-hrsx+1)) ',' num2str(round(nc(i,2)-hrsy+1)) ',' num2str(round(nc(i,3)-hrsz+1)) '_' num2str(round(nc(i,1)+hrsx)) ',' num2str(round(nc(i,2)+hrsy)) ',' num2str(round(nc(i,3)+hrsz)) '_SpectralMixed.tif'];
            
            copyfile([rt fnO], [rtcopy fnN]);
            FE(i) = 1;
        end
    end
    parfor_progress;
end
parfor_progress(0);