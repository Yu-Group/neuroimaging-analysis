function GU_ExM_JaneliaCluster_Mosuse_PunctaDensity_LARGEVOL(fnv1, fnv2)



GU_ExM_calcPunctaDensityVol(fnv1, fnv2,...
        'FindLocalMaxima', false, 'Threshold', [1450,2200], 'MinThreshold', [1450,2200],'Verbose', true,'OTSUGaussKernel', 0.5,...
        'MinVoxelVolume', [516, 1008],'ExpansionFactor', 3.62);