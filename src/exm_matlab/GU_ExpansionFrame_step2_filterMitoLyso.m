function GU_ExpansionFrame_step2_filterMitoLyso(data)
se = strel('disk',7);
seM = strel('disk',12);
se2 = strel('disk',12);
i = 1;
fn{i} = [data(i).source 'Analysis' filesep 'holes.tif'];
md{i} = [data(i).source 'Analysis' filesep 'MaskedData.tif'];
m{i} = [data(i).source 'Analysis' filesep 'Mask.tif'];
mfn{i} = data(i).framePaths{1}{2};
lfn{i} = data(i).framePaths{1}{3};

MData = logical(readtiff(md{i}));
maskdata = logical(readtiff(m{i}));
%     mask = logical(readtiff(fn{i}));

% keep only the largest object
[lholes, ~] = bwlabeln(MData);
CC = bwconncomp(lholes, 26);
csize = cellfun(@numel, CC.PixelIdxList);
idx = csize>=max(csize);
CC.NumObjects = sum(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
MData = labelmatrix(CC)~=0;

maskdata(MData==0) = 0;
mask = MData - maskdata;
mask(mask<0) = 0;
mask = logical(mask);
clear maskdata MData

% sorts labels by ascending sizes
[lholes, ~] = bwlabeln(mask);
[smask] = GU_LabelsSortBySize(lholes);
writetiff(uint8(mask), [data(i).source 'Analysis' filesep 'holes.tif']);
writetiff(uint16(smask), [data(i).source 'Analysis' filesep 'SortedHoles.tif']);

CC = bwconncomp(lholes, 26);
csize = cellfun(@numel, CC.PixelIdxList);
holeVol = csize*data(i).pixelSize^2*data(i).dz; %%%%%%% calc volume of all holes
%  rp3 = regionprops3(mask);
%mito
mitoIm = readtiff(mfn{i});
mito = zeros(size(mask), 'uint16');
mito(mask>0) = mitoIm(mask>0);
mito = imdilate(mito, seM);
cc = bwconncomp(logical(mito));
numPixels = cellfun(@numel,cc.PixelIdxList);
mitoSig = zeros(1,numel(numPixels));
for j = 1:numel(numPixels)
    mitoSig(j) = sum(mito(cc.PixelIdxList{j}));
end
nMS = mitoSig./numPixels;% normalized mito signal
T = thresholdOtsu(nMS(nMS < prctile(nMS, 95)));
idx = nMS > T;
mitoMask = zeros(size(mask), 'uint8');
cc.PixelIdxList = cc.PixelIdxList(idx);

for j = 1:numel(cc.PixelIdxList)
    mitoMask(cc.PixelIdxList{j}) = 1;
end
mitoMask(mask==0) = 0;
[lholes, ~] = bwlabeln(mitoMask);
CC = bwconncomp(lholes, 26);
csize = cellfun(@numel, CC.PixelIdxList);
mitoholeVol = csize*data(i).pixelSize^2*data(i).dz;%%%%%%% calc volume of mito holes
mitorp3 = regionprops3(mitoMask);
writetiff(mitoMask, [data(i).source 'Analysis' filesep 'MitoMask_new.tif']);
clear mitoIm mito
% lyso
if ~isempty(lfn{i})
    lmask = imerode(mask, se);
    lysoIm = readtiff(lfn{i});
    lyso = zeros(size(lmask), 'uint16');
    lyso(lmask>0) = lysoIm(lmask>0);
    
    cc = bwconncomp(lmask);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    lysoSig = zeros(1,numel(numPixels));
    for j = 1:numel(numPixels)
        lysoSig(j) = sum(lyso(cc.PixelIdxList{j}));
    end
    nLS = lysoSig./numPixels;% normalized mito signal
    T = thresholdOtsu(nLS(nLS<prctile(nLS, 99.7)));
    idx = nLS > T; %| (lysoSig > prctile(lysoSig,99.5) & lysoSig < prctile(lysoSig,99.9));
    lysoMask = zeros(size(mask), 'uint8');
    cc.PixelIdxList = cc.PixelIdxList(idx);
    for j = 1:numel(cc.PixelIdxList)
        lysoMask(cc.PixelIdxList{j}) = 1;
    end
    lysoMask = imdilate(lysoMask, se2);
    lysoMask(mask==0) = 0;
    lysoMask = lysoMask - mitoMask;
    lysoMask(lysoMask<0) = 0;
    lysoMask = logical(lysoMask);
    [lholes, ~] = bwlabeln(lysoMask);
    CC = bwconncomp(lholes, 26);
    csize = cellfun(@numel, CC.PixelIdxList);
    lysoholeVol = csize*data(i).pixelSize^2*data(i).dz;%%%%%%% calc volume of mito holes
    lysorp3 = regionprops3(lysoMask);
    writetiff(uint8(lysoMask), [data(i).source 'Analysis' filesep 'LysoMask_new.tif']);
    clear lysoIm lyso
end
save([data(i).source 'parameters_rp3_lyso_mito_all'], 'lysoholeVol', 'lysorp3', 'mitorp3', 'mitoholeVol','holeVol')

