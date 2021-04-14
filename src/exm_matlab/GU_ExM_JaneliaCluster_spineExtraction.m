function GU_ExM_JaneliaCluster_spineExtraction(ex)

% Gokul Upadhyayula, Nov 2017
if ischar(ex)
    ex = str2double(ex);
end

ef = 4.04;
px = .097; %in um
pz = 0.25;
zAniso = pz/px;

rt = '/groups/betzig/betziglab/Rui/160126_Sample4_Y10_YFP_HB/ch0_Analysis_1ch/1008/';
cd(rt)
fn = dir(['*part*.tif']);
fn = {fn.name}';
% fn = 'ch0_1008_clean.tif';


GU_ExM_extractSpineParameters([rt fn{ex}], 'zAniso', zAniso, ...
    'PixelSize', 0.097, 'MinBranchLength', 10, 'RemoveOverlap', 0, 'CropScalingFactor', 5,'MinVoxelVolume', 0)