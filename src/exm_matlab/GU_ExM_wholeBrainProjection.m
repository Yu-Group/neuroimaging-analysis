clear rt
% rt{1} = 'M:\TIFF\ch0_DANmem_ch1_nc82ALL_masked_rescaled_v3\';
rt{1} = 'M:\TIFF\ch0_DANmem_ch1_nc82DANassociated_masked_rescaled_v3\';
srt{1} = [rt{1} 'proj' filesep];
% srt{2} = [rt{2} 'crop' filesep];
mkdir(srt{1});
% mkdir(srt{2});
d = dir([rt{1} '*.tif']);
d = {d.name};
% crop coordinates
danim = zeros(27964,15055, 'uint16');
nc82im = zeros(27964,15055, 'uint16');

tic
parfor_progress(numel(d)/2);
for i = 1:numel(d)/2
    im0 = readtiff([rt{1} 'z' num2str(i-1,'%04d') '_c0_t0.tif']);
    im0(im0>=53247) = 0;
    im1 = readtiff([rt{1} 'z' num2str(i-1,'%04d') '_c1_t0.tif']);
    im1(im1>=53247) = 0;
    danim = max(danim, im0);
    nc82im = max(nc82im, im0);
    parfor_progress;
end
parfor_progress(0);
toc
writetiff(danim, [srt{1} 'DANmaxproj.tif'])
writetiff(nc82im, [srt{1} 'nc82maxproj.tif'])
