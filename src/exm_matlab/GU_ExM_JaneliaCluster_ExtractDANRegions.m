function GU_ExM_JaneliaCluster_ExtractDANRegions(ex)

% Gokul Upadhyayula, Oct 2018
% downsample fly brain ; generated a perspective downsampled z-stack
danrt = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/DAN/';
nc82rt = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/clean_DANAssociated_masked_rad7/';

fn = dir([danrt '*.tif']);
fn = {fn.name};


% stitched data dimensions
danMRrt_perspective = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/DAN_maskedRegionsPerspective/';
for i = 1:15
    if ~exist([danMRrt_perspective num2str(i)], 'dir')
        mkdir([danMRrt_perspective num2str(i)])
    end
end

nc82MRrt_perspective = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/nc82_DANAssociated_masked_rad7_maskedPers/';
for i = 1:15
    if ~exist([nc82MRrt_perspective num2str(i)], 'dir')
        mkdir([nc82MRrt_perspective num2str(i)])
    end
end
maskFrames = '/groups/betzig/betziglab/Gokul/ExLLSM_GokulWorkspace/ExM_FlyBrain/maskframes/';

mx = 15055;
my = 27964;
mz = 6647;

% mask dimensions
sx = 940;
sy = 1747;
sz = 1661;

% mask downsampled
fx = 16.032;%my/sy;
fy = 16.0155;%my/sy;
fz = mz/sz;

MPR = 0.4;% max perspective reduction over depth
PRS = MPR/6646; % perspective reduction step size

z = str2double(ex); % the z-slice to separate
ii = 6646-z;


n = z;
im = readtiff([danrt num2str(n) '.tif']); %dan
im2 = readtiff([nc82rt num2str(n) '.tif']); %nc82
maskslice = readtiff([maskFrames num2str(ceil(n/fz)) '.tif']);
so = size(im);
mask = imresize(maskslice, so, 'nearest');

% separate data based on masks
for j = 1:15
    %dan
    mim = zeros(size(im), 'uint16');
    mim(mask==j) = im(mask==j);
    mim = imresize(mim, 1-PRS*ii);
    pim = zeros(so, 'uint16');
    sr = size(mim);
    dr = so-sr;
    dr2 = round((dr)/2);
    pim(dr2(1)+1:end-(dr(1)-dr2(1)),dr2(2)+1:end-(dr(2)-dr2(2))) = mim;
    writetiff(pim, [danMRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif'], 'Compression', 'lzw');
    
    % nc82
    mim = zeros(size(im2), 'uint16');
    mim(mask==j) = im2(mask==j);
    mim = imresize(mim, 1-PRS*ii);
    pim = zeros(so, 'uint16');
    sr = size(mim);
    dr = so-sr;
    dr2 = round((dr)/2);
    pim(dr2(1)+1:end-(dr(1)-dr2(1)),dr2(2)+1:end-(dr(2)-dr2(2))) = mim;
    writetiff(pim, [nc82MRrt_perspective filesep num2str(j) filesep  num2str(n) '.tif'], 'Compression', 'lzw');
end



