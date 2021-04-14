function GU_ExM_WriteRotatedCoordinates(rt, mx, my, omz, zxRatio, theta, fileSufix)
%% rotated size

mz = omz*zxRatio;
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


for ex = 1:nn-1
    cc(ex,1:3) = [(xmax(ex)-(xmax(ex)-xmin(ex)+1)/2), (ymax(ex)-(ymax(ex)-ymin(ex)+1)/2), (zmax(ex)-(zmax(ex)-zmin(ex)+1)/2)]; % current center x,y,z
end

occ = cc;

mz = omz;
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
clear nc
zf = cotd(abs(theta));
yf = 1;
xf = cosd(abs(theta)) + tand(abs(theta))*sind(abs(theta));

mz = omz*zxRatio;
cx = mx/2;
cy = my/2;
cz = mz/2;

[nc, R, t] = AxelRot(cc', theta*-1, [0,1,0], [cx+1,cy+1,cz+1]);
nc = nc';

close all; figure, scatter3(nc(:,1),nc(:,2),nc(:,3)); hold on; scatter3(occ(:,1),occ(:,2),occ(:,3)); hold on
nc(:,1) = nc(:,1)/xf;
nc(:,3) = nc(:,3)/zf;
scatter3(nc(:,1),nc(:,2),nc(:,3));
%% rename files based on the new center point

rtcopy = [rt 'nonEmpty/'];
rtcopy2 = [rt 'allRenamedTiles/'];
mkdir(rtcopy);
mkdir(rtcopy2);
rtmat = [rt 'params/'];
hs = [750 750 750]/2; % half size
parfor_progress(nn-1);
FE = zeros(1, nn-1);

parfor i = 1:nn-1
    fnO = [num2str(xmin(i)) ',' num2str(ymin(i)) ',' num2str(zmin(i)) '_' num2str(xmax(i)) ',' num2str(ymax(i)) ',' num2str(zmax(i)) fileSufix '.tif'];
    if exist([rt fnO], 'file') 
        if exist([rtmat fnO(1:end-3) 'mat'], 'file')
            Proceed = load([rtmat fnO(1:end-3) 'mat'], 'nVoxOccupied');
        else 
            Proceed.nVoxOccupied = 0;
        end
        IMinfo = imfinfo([rt fnO]);
        hrsx = IMinfo(1).Width/2;
        hrsy = IMinfo(1).Height/2;
        hrsz = numel(IMinfo)/2;
        fnN = [num2str(nc(i,1)-hrsx+1) ',' num2str(nc(i,2)-hrsy+1) ',' num2str(nc(i,3)-hrsz+1) '_' num2str(nc(i,1)+hrsx) ',' num2str(nc(i,2)+hrsy) ',' num2str(nc(i,3)+hrsz) '_Rotated.tif'];
        if Proceed.nVoxOccupied > 0
            copyfile([rt fnO], [rtcopy fnN]);
            FE(i) = 1;
        end
        copyfile([rt fnO], [rtcopy2 fnN]);
    end
    parfor_progress;
end
parfor_progress(0);
sum(FE)