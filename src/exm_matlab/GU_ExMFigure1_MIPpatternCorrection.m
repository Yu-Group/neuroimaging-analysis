rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/Figures_ExLLSM/Figure1_Schematics+comparison/MIPs/';
cd(rt)

% fn = 'LLSM(SS)_interpX_MAX_ch0_Z2600-3150-1.tif';
fn = 'LLSM(SS)_MAX_ch0_Z2600-3150.tif';
im = readtiff(fn);

%%

close all

ker = 50;
% xz = mean(filterGauss2D(im(:,end-400:end),ker), 2);
xz = mean(filterGauss2D(im(:,end-2000:end-1000),ker), 2);
figure, plot(xz)

%%

%figure, imagesc(im)
m = 700;
nXZ = xz;
nXZ(nXZ > m) = m;
% nXZ = filterGauss1D(nXZ,ker);
nXZ = mat2gray(nXZ);
% nXZ = nXZ/max(nXZ);
figure, plot(nXZ)
cim = repmat(nXZ, [1, size(im,2)]);
%%
s = 700;
cim2 = cim*s;
cim2 = filterGauss2D(cim2,ker);
%
% cim2(cim<0.75) = 0.75;
nim = double(im)- cim2;
nim(nim<0) = 0;
nim = uint16(nim);

figure
imagesc(nim)
%%
writetiff(nim, [date '_' fn]);