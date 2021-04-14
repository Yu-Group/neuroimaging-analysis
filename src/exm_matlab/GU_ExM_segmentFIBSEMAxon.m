% GU - for expansion - segment FIBSEM Axons

% load raw cropped data
rt = 'E:\Gokul\ShuData_HPF\';
cd(rt)
im = readtiff([rt '2017-1019-16x16x16-center-crop.view.tif']);
lim = readtiff([rt 'axonSegmentedLabels.tif']);

%%
nA = unique(lim(lim(:)> 0)); % # of axons
ed = strel('sphere', 5); % erosion and dilation factor

% loop over the segmented regions to segment axon

% for i = 1:numel(nA)
    cm = zeros(size(lim), 'uint8');
    cm(ismember(lim(:),[1,3,4,6])) = 1;% == nA(i);
    cm = imdilate(cm, ed);

    for z = 1:size(cm,3)
        cm(:,:,z) = imfill(cm(:,:,z), 8, 'holes');
    end
    for x = 1:size(cm,2)
        cm(:,x,:) = imfill(squeeze(cm(:,x,:)), 8, 'holes');
    end
    for y = 1:size(cm,1)
        cm(y,:,:) = imfill(squeeze(cm(y,:,:)), 8, 'holes');
    end
    
    cm = imerode(cm, ed);
    maskedIm = zeros(size(cm), 'uint8');
    maskedIm(logical(cm)) = im(logical(cm));
    
    
    writetiff(maskedIm, 'maskedIm.tif');
    ecm = imerode(cm,ed);
    emaskedIm = zeros(size(ecm), 'uint8');
    emaskedIm(logical(ecm)) = im(logical(ecm));
     writetiff(emaskedIm, 'erodedmaskedIm.tif');
     
      mhaWriter('data.mha', im, [16,16,16], 'uint8');
     mhaWriter('maskedData.mha', maskedIm, [16,16,16], 'uint8');
     mhaWriter('mask.mha', cm, [16,16,16], 'uint8');
% end

%% read curated mask

[curatedMask, ~] = mhaReader('Curated_mask.mha');
maskL = uint8(bwlabeln(curatedMask));
 mhaWriter('labeledCuratedMask.mha', maskL, [16,16,16], 'uint8');
%%

