function GU_ExM_CropWholeFlyBrainSlice(fp, xmax, ymax)

% Gokul Upadhyayula, Oct 2017

if ischar(xmax)
    xmax = str2double(xmax);
end
if ischar(ymax)
    ymax = str2double(ymax);
end
[p, fn, ~] = fileparts(fp);

% xmax = 15055;
% ymax = 27964;

im = readtiff(fp);
im = im(1:ymax, 1:xmax);
writetiff(uint16(im), fp);

maxVoxOccupancy = sum(im(:)>0)/numel(im(:));
idx = find(logical(im));
[ny,nx,nz] = size(logical(im));
[yi,xi,zi] = ind2sub([ny,nx,nz], idx);

yi_min = min(yi);
yi_max = max(yi);
xi_min = min(xi);
xi_max = max(xi);
zi_min = min(zi);
zi_max = max(zi);

if ~exist([p filesep 'SliceProperties'], 'dir')
    mkdir([p filesep 'SliceProperties']);
end

save([p filesep 'SliceProperties' filesep fn '.mat'], 'maxVoxOccupancy', 'yi_min', 'yi_max', 'xi_min', 'xi_max','zi_min', 'zi_max');