% GU_ExM_MBON_workspace

rt = 'D:\GoogleDrive_DataTransfer\JaneliaSharedwithGokul\ExpansionRelated\Wholebrainsynapse\Ex23_MBON_dedritesvsaxons_z0p18';
cd(rt)
axonfn{1} = [rt filesep 'Ex1_MBON-a1-axon-boutons\Raw\a1-axon-boutons-ch1-nc82.tif'];
axonfn{2} = [rt filesep 'Ex1_MBON-a1-axon-boutons\Raw\a1-axon-boutons-ch0-membrane.tif'];
labelMask{1} = [rt filesep 'Ex1_MBON-a1-axon-boutons\Mask\a1-axon-boutons-ch0-mask.tif'];
CleanIm{1} = [rt filesep 'Ex1_MBON-a1-axon-boutons\Raw\Analysis\1008\a1-axon-boutons-ch0-membrane_clean.tif'];
denfn{1} = [rt filesep 'Ex2_MBON-a1-dendrites\Raw\a1-dendrites-ch1-nc82.tif'];
denfn{2} = [rt filesep 'Ex2_MBON-a1-dendrites\Raw\a1-dendrites-ch0-membrane.tif'];
labelMask{2} = [rt filesep 'Ex2_MBON-a1-dendrites\Mask\a1-dendrites-ch0-mask.tif'];
CleanIm{2} = [rt filesep 'Ex2_MBON-a1-dendrites\Raw\Analysis\1008\a1-dendrites-ch0-membrane_clean.tif'];
%%
GU_ExM_calcPunctaDensityVol(axonfn{1}, axonfn{2},...
    'FindLocalMaxima', true,'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', [246, 1008],'ExpansionFactor', 3.99,'minAZdia', 0.1);

GU_ExM_calcPunctaDensityVol(denfn{1}, denfn{2},...
    'FindLocalMaxima', true,'Verbose', true,'OTSUGaussKernel', 0.5,...
    'MinVoxelVolume', [246, 1008],'ExpansionFactor', 3.99,'minAZdia', 0.1);

%% get surface area and volume
px = 0.097;
zAniso = 0.18/px;
ef = 3.99;

for k = 1:numel(CleanIm)
    tic
    if exist([labelMask{k}(1:end-4) '_Seg_filled_masked.tif'], 'file')
        mask = readtiff([labelMask{k}(1:end-4) '_Seg_filled_masked.tif']);
    else
        if exist( [CleanIm{k}(1:end-4) '_filled.tif'], 'file')
            mask = logical(readtiff( [CleanIm{k}(1:end-4) '_filled.tif']));
        else
            mask = readtiff(CleanIm{k}); toc
            mask = logical(mask);
            
            % fill gaps
            tic
            for z = 1:size(mask,3)
                mask(:,:,z) = imfill(mask(:,:,z), 8, 'holes');
            end
            for x = 1:size(mask,2)
                mask(:,x,:) = imfill(squeeze(mask(:,x,:)), 8, 'holes');
            end
            for y = 1:size(mask,1)
                mask(y,:,:) = imfill(squeeze(mask(y,:,:)), 8, 'holes');
            end
            for z = 1:size(mask,3)
                mask(:,:,z) = imfill(mask(:,:,z), 8, 'holes');
            end
            for x = 1:size(mask,2)
                mask(:,x,:) = imfill(squeeze(mask(:,x,:)), 8, 'holes');
            end
            for y = 1:size(mask,1)
                mask(y,:,:) = imfill(squeeze(mask(y,:,:)), 8, 'holes');
            end
            toc
            
            % cm = imclose(cm, binarySphere(HSIZE));
            writetiff(uint8(mask), [CleanIm{k}(1:end-4) '_filled.tif']);
        end
        
        
        % limit content of the mask to just the hand labeled regions
        
        [tmpmask, ~] = imread_big(labelMask{k});
        mask(~logical(tmpmask)) = 0;
        clear tmpmask
        writetiff(uint8(mask), [labelMask{k}(1:end-4) '_Seg_filled_masked.tif']);
    end
    
    % crop to avoid long interpolation times
    tic
    idx = find(mask);
    [ny,nx,nz] = size(mask);
    [yi,xi,zi] = ind2sub([ny,nx,nz], idx);
    b = 5; % boundary pixels for cropping
    xa = max(min(xi)-b,1):min(max(xi)+b,nx);
    ya = max(min(yi)-b,1):min(max(yi)+b,ny);
    za = max(min(zi)-b,1):min(max(zi)+b,nz);
    mask = mask(ya,xa,za);
    toc
    
    tic
    % interpolate
    [ny,nx,nz] = size(mask);
    [y,x,z] = ndgrid(1:ny,1:nx,1:nz);
    [Y,X,Z] = ndgrid(1:ny,1:nx,1:1/zAniso:nz);
%     imask = interp3(x,y,z,double(mask),X,Y,Z,'nearest');
    imask = interp3((mask),X,Y,Z,'nearest');
    toc
    
    % write out cropped cells
    % writetiff(uint8(imask), [WP filesep name ext]);
    % calculate vol & area
    perim = bwperim(imask);
    vol(k) = sum(imask(:)) * (px/ef)^3;
    area(k) = sum(perim(:)) * (px/ef)^2;
end
axonVol = vol(1);
DenVol = vol(2);
axonArea = area(1);
DenArea = area(2);