function GU_ExM_JaneliaCluster_PunctaDensity_LARGEVOL(fnv1, fnv2, thresh1, thresh2)
% requires precalculation of the OTSU threshold
% Gokul Upadhyayula, 2017
if ischar(thresh1)
    thresh1 = str2double(thresh1);
end

if ischar(thresh2)
    thresh2 = str2double(thresh2);
end

GU_ExM_calcPunctaDensityVol(fnv1, fnv2,...
        'FindLocalMaxima', false, 'Threshold',[thresh1, thresh2] ,'MinThreshold', [thresh1, thresh2],'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [246, 1008],'ExpansionFactor', 3.99);