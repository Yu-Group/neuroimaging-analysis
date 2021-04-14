function GU_ExM_Janelia_MouseVisualCortex_Mask_SpectralMix(ex, df)

% Gokul Upadhyayula, Nov 2017

if ischar(df)
    df = str2double(df);
end

if ischar(ex)
    ex = str2double(ex);
end

se = strel('sphere', df);
fprintf('calculating file names...'); tic

rt0 = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/slice-tiff_chroma_cleaned_rotated/ch0/';
rt1 = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/slice-tiff_chroma_cleaned_rotated/ch1/';
rt2 = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/slice-tiff_chroma_cleaned_rotated/ch2/';


% read tiff
mx = 12228;
my = 39519;

% read tiff
fprintf('reading volume...'); tic
im0 = zeros(my,mx,11, 'logical');
im1 = zeros(my,mx,11, 'uint16');
im2 = zeros(my,mx,11, 'uint16');
nf = [(ex-5):(ex+5)]-1;

for i = 1:numel(nf)
    im0(:,:,i) = logical(readtiff([rt0 [num2str(nf(i)) '.tif']]));
    im1(:,:,i) = readtiff([rt1 [num2str(nf(i)) '.tif']]);
    im2(:,:,i) = readtiff([rt2 [num2str(nf(i)) '.tif']]);
end
toc

fprintf('dilating YFP volume...'); tic
im0 = imdilate(im0,se);toc
im1(im1(:)>5000) = 5000;
im2(im2(:)>30000) = 30000;

fprintf('initializing new volume...'); tic
im1_masked = zeros(size(im0(:,:,3)), 'uint16');
im2_masked = zeros(size(im0(:,:,1)), 'uint16');
toc

fprintf('masking data...'); tic
im1_masked(im0(:,:,6)) = uint16(scaleContrast(im1(im0(:,:,6)), [0, 5000], [0, 32760]));
im2_masked(im0(:,:,6)) = uint16(scaleContrast(im2(im0(:,:,6)), [0, 30000], [0, 32760]));
toc

fprintf('rescaling data...'); tic

im1_masked(~im0(:,:,6) & im1(:,:,6)>0 ) = uint16(scaleContrast(im1(~im0(:,:,6) & im1(:,:,6)>0 ), [0, 5000], [32768, 65535]));toc
im2_masked(~im0(:,:,6) & im2(:,:,6)>0 ) = uint16(scaleContrast(im2(~im0(:,:,6) & im2(:,:,6)>0 ), [0, 30000], [32768, 65535]));toc


rt1save = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/slice-tiff_chroma_cleaned_rotated/ch1/YFPMaskedSpectralMixed/';
rt2save = '/groups/betzig/betziglab/4Stephan/171123_Mousevisualcortex/slice-tiff_chroma_cleaned_rotated/ch2/YFPMaskedSpectralMixed/';

if ~exist(rt1save, 'dir')
   mkdir(rt1save)
end

if ~exist(rt2save, 'dir')
   mkdir(rt2save)
end


fprintf('saving data...'); tic
writetiff(im1_masked, [rt1save num2str(ex) '.tif']);
writetiff(im2_masked, [rt2save num2str(ex) '.tif']);toc
end