function GU_ExM_JaneliaCluster_MergeMouseSomato_stitchedSlices(ex, minVol)

% Gokul Upadhyayula, Dec 2017

rt_mem = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis/1008/RotatedStacks_cubicInterp/nonEmpty/slice-tiff/ch0/';
rt_Homer = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246/RotatedStacks_cubicInterp/nonEmpty/slice-tiff/ch0/';
rt_HomerV2Ass = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch2/Analysis/246_V2Ass/RotatedStacks/nonEmpty/slice-tiff/ch0/';
rt_save = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed_cubicInterp_Clean/';

mx = 10629;
my = 66496;

if ~exist(rt_save, 'dir')
    mkdir(rt_save)
end

s = 15; % slab half thickness

RangeIn = {[0 15000], [0 30000], [0 30000]};
RangeOut = {[0 32767], [32768 49151], [49152 2^16-1]};

if ischar(ex)
    ex = str2double(ex);
end

if ischar(minVol)
    minVol = str2double(minVol);
end


% read tiff
fprintf('reading volume...'); tic
im = zeros(my,mx,s*2, 'uint16');
imAss = zeros(my,mx,s*2, 'uint16');
nf = [(ex-s):(ex+s)]+1;
for i = 1:numel(nf)
    im(:,:,i) = readtiff([rt_mem [num2str(nf(i)) '.tif']]);
    imAss(:,:,i) = readtiff([rt_HomerV2Ass [num2str(nf(i)) '.tif']]);
    
end
toc

im_all{1} = im(:,:,s); % mem

% clean membrane
if minVol > 0
        fprintf('discarding small noise voxels for membrane channel...');
    CC = bwconncomp(logical(im), 26);
    csize = cellfun(@numel, CC.PixelIdxList);
    idx = csize>=minVol;
    CC.NumObjects = sum(idx);
    CC.PixelIdxList = CC.PixelIdxList(idx);
    mask = labelmatrix(CC)~=0;
    im = uint16(im) .* uint16(mask);
end

fprintf('geneating labels...'); tic

lholesHomer = bwlabeln(logical(imAss));
overlap = lholesHomer(logical(im));
idx = ismember(lholesHomer,unique(overlap));

clear lholesHomer

imCleanHomer = zeros(size(im), 'uint16');
imCleanHomer(idx) = imAss(idx);

im_all{2} = imCleanHomer(:,:,s); % associated Homer

clear im idx mask overlap imAss

imNonAssHomer = readtiff([rt_Homer [num2str(ex) '.tif']]);
if sum(im_all{2}(:)) > 0
    imNonAssHomer(im_all{2}(:) > 0) = 0;
end

im_all{3} = imNonAssHomer;

clear imNonAssHomer imCleanHomer


% rescale data
fprintf('Rescaling volumes...')
tic
for i = 1:numel(im_all)
   im_all{i}(im_all{i}(:)> RangeIn{i}(2)) = RangeIn{i}(2); 
   im_all{i}(im_all{i}(:)> 0) = uint16(scaleContrast(im_all{i}(im_all{i}(:)> 0), RangeIn{i}, RangeOut{i}));
end
toc

fprintf('Max Projecting volumes...')
tic
% max project volume
maxProj = zeros(size(im_all{1}),'uint16');
for i = 1:numel(im_all)
   maxProj = max(maxProj, im_all{i});
end
toc

fprintf('Saving Mixed Volumes...')
tic
writetiff(maxProj, [rt_save num2str(ex) '.tif']);
toc