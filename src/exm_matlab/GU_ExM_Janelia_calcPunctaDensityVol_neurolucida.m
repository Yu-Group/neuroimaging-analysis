function GU_ExM_Janelia_calcPunctaDensityVol_neurolucida(vol1, vol2)




GU_ExM_calcPunctaDensityVol(vol1, vol2,...
            'FindLocalMaxima', false, 'MinThreshold', [1450,2200],'Verbose', true,'OTSUGaussKernel', 0.5,...
            'MinVoxelVolume', [246, 1008],'ExpansionFactor', 3.62); % sigma 3.9