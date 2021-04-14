function GU_ExM_JaneliaCluster_GenerateMaskedVolumes(fn, sliceN)

% Gokul Upadhyayula, 2017

% clean_DAN_rt = 'X:\4Stephan\171016_Flybrainsynapse\171016_100nmcluster\ch0\Analysis\clean';
% clean_nc82_rt = 'X:\4Stephan\171016_Flybrainsynapse\171016_100nmcluster\ch1\Analysis\clean';
% nc82_Data_rt = 'X:\4Stephan\171016_Flybrainsynapse\171016_100nmcluster\ch1\Analysis\DataStructures';
% matfile = [ch1rt filesep 'Analysis' filesep 'DataStructures' filesep fn '_densityData.mat'];
% cleanDanpath = [ch0rt filesep 'Analysis' filesep 'clean' filesep fn '_clean.tif'];

% min value = 0.001%; max value is 99%
minRange_DAN = 673;
maxRange_DAN = 3687;

minRange_nc82 = 524;
maxRange_nc82 = 2930;

sliceN = str2double(sliceN);

% stitched data dimensions

mx = 15055;
my = 27964;
mz = 6647;

% mask dimensions
sy = 1747;
sx = 940;
sz = 1661;

% mask downsampled
fy = 16.0155;%my/sy;
fx = 16.032;%mx/sx;
fz = mz/sz;


% Paths
DANrt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_DANAssociated/slice-tiff/ch0/'; %NOT corrected - need to use the corrected values
nonDANrt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_nonDANAssociated/slice-tiff/ch0/';

DANassMASKrt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_DANAssociated_masked/'; % these paths have data that are corrected for the DAN puncta
% nonDANassMASKrt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_nonDANAssociated_masked/';

DANmemrt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/clean/slice-tiff/ch0/';
nc82rt = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/slice-tiff/ch0/';

nonDANmaskedsaveRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/nc82_nonDANassociated_masked_rescaled_v3/';
DANmaskedsaveRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/nc82DANassociated_masked_rescaled_v3/';
nch82DAN_maskedsaveRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/nc82ALL_masked_rescaled_v3/';
DANmemMaskedSaveRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/DANclean_masked_rescaled_v3/';

mkdir(nonDANmaskedsaveRT);
mkdir(DANmaskedsaveRT);
mkdir(nch82DAN_maskedsaveRT);
mkdir(DANmemMaskedSaveRT);

% read mask
fprintf('reading mask...'); tic
mask = readtiff('/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/mask/CombinedMask.tif');
toc
maskslice = mask(:,:,ceil(sliceN/fz));

% read data
fprintf('reading data...'); tic
DANmask = logical(readtiff([DANassMASKrt fn]));

nonDAN = readtiff([nonDANrt fn]);

nonDAN(nonDAN>maxRange_nc82) = maxRange_nc82;
nonDAN(DANmask) = 0;

DAN = readtiff([DANrt fn]);

DAN(DAN>maxRange_nc82) = maxRange_nc82;
DAN(~DANmask) = 0;

DANmem = readtiff([DANmemrt fn]);
DANmem(DANmem>maxRange_DAN) = maxRange_DAN;

nc82 = readtiff([nc82rt fn]);
nc82(nc82>maxRange_nc82) = maxRange_nc82;
toc

fprintf('interpolating mask...'); tic
% interpolate mask
[y,x] = ndgrid(1:sy,1:sx);
[Y,X] = ndgrid(1:1/fy:sy,1:1/fx:sx);
imaskslice = interp2(x,y,double(maskslice),X,Y,'nearest');
imaskslice = uint8(imaskslice);
toc


% initialize new masked contrasted regions
nonDANmasked = zeros(size(imaskslice), 'uint16');
DANmasked = zeros(size(imaskslice), 'uint16');

nc82masked = zeros(size(imaskslice), 'uint16');
DANMEMmasked = zeros(size(imaskslice), 'uint16');

% loop through the 15 regions and rescale intensities
fprintf('rescaling intensity values...'); tic
u = unique(maskslice(:))';
for j = u
    nonDANmasked(imaskslice==j & nonDAN > 0) =  scaleContrast(nonDAN(imaskslice==j  & nonDAN > 0), [minRange_nc82, maxRange_nc82], [(j-1)*2^12, j*2^12]);
    DANmasked(imaskslice==j & DAN > 0) =  scaleContrast(DAN(imaskslice==j & DAN > 0), [minRange_nc82, maxRange_nc82], [(j-1)*2^12, j*2^12]);
    
    nc82masked(imaskslice==j & nc82 > 0) =  scaleContrast(nc82(imaskslice==j & nc82 > 0), [minRange_nc82, maxRange_nc82], [(j-1)*2^12, j*2^12]);
    DANMEMmasked(imaskslice==j & DANmem > 0) =  scaleContrast(DANmem(imaskslice==j & DANmem > 0), [minRange_DAN, maxRange_DAN], [(j-1)*2^12, j*2^12]);
end
toc

fprintf('writing masked data...'); tic
writetiff(uint16(nonDANmasked), [nonDANmaskedsaveRT fn]);
writetiff(uint16(DANmasked), [DANmaskedsaveRT fn]);

writetiff(uint16(nc82masked), [nch82DAN_maskedsaveRT fn]);
writetiff(uint16(DANMEMmasked), [DANmemMaskedSaveRT fn]);
toc