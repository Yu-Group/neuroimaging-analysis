% GU_ExM_GRatioWorkspace
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/G-RatioData/Raw/AxonSegements for statistics_newstitch/';
cd(rt)
load('dataLoad.mat')
%%

% write mha
for i = 1:numel(data)
    mbp = readtiff(data(i).framePaths{1}{2});
   mhaWriter([data(i).source 'mbp.mha'], mbp, [0.097, 0.097, 0.18], 'short'); 
    i
end

%% calc local GRatio
parfor i = 2:numel(data)
    mask = mhaReader([data(i).source 'mask.mha']);
    GU_ExM_calcAxonMyelinGRatio_frame(data(i), 'Mask', mask, 'TValue', [0.1, 0.5], 'UseCanny', false, 'LargestObject',true, 'MyelinOffset', 3,'OTSUMaxPer', 99.9);
end

%% 
for i = 2:numel(data)
    mask = mhaReader([data(i).source 'mask.mha']);
    mbp = readtiff(data(i).framePaths{1}{2});
    axon = readtiff(data(i).framePaths{1}{1});
    
    mbp(~logical(mask)) = 0;
    axon(~logical(mask)) = 0;
    writetiff(mbp,[data(i).source 'masked_MBP.tif']);
    writetiff(axon,[data(i).source 'masked_Axon.tif']);
end