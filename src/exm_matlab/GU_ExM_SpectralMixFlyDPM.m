function GU_ExM_SpectralMixFlyDPM(ex)

% this function combines volumes in different regions of the 16bit volume
% this assumes that the volumes are spatially non-overlapping and max projected
% Gokul Upadhyayula, Oct 2017
ex = str2double(ex);
RangeIn = {[0 32767], [0 32767]};
RangeOut = {[0 32767], [32768 2^16-1]};


% rt{1} = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch0/Analysis/1008/RotatedStacks/allRenamedTiles/slice-tiff/ch0/';
rt{1} = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch1/Analysis/246_V2Ass/RotatedStacks/allRenamedTiles/slice-tiff/ch0/';
rt{2} = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/ch1/Analysis/246_V2nonAss/RotatedStacks/allRenamedTiles/slice-tiff/ch0/';
saveFolder = '/groups/betzig/betziglab/4Stephan/160919_DPM3/Processing/SpectrallyMixed/';

for i = 1:numel(rt)
    fn{i} = [rt{i} num2str(ex) '.tif'];
end

SavePath = [saveFolder num2str(ex) '_SpectralMixed.tif'];
tic
GU_CombineVolumesIntoSingleBitDepth(fn, RangeIn, RangeOut, SavePath);
toc

