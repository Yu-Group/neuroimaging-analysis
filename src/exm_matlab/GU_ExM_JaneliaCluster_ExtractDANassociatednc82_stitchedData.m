function GU_ExM_JaneliaCluster_ExtractDANassociatednc82_stitchedData(ex)

% Gokul Upadhyayula, Nov 2017
se = strel('sphere', 7);
ex = str2double(ex);
fprintf('calculating file names...'); tic
% get file names
mx = 15055;
my = 27964;
% mz = 6647;


matRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/';
% nc82RT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/slice-tiff/ch0/';
nc82RT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/nc82ALL_masked_rescaled_v3/';
nc82DANsaveRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_DANAssociated_masked/';
nc82nonDANsaveRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_nonDANAssociated_masked/';
fn_mat = ['20171025_CombinedSynapses_withParams.mat'];
fn_nc82 = [num2str(ex-1) '.tif'];
toc

% read tiff
fprintf('reading volume...'); tic
im = zeros(my,mx,11, 'uint16');
nf = [(ex-5):(ex+5)]-1;
for i = 1:numel(nf)
    im(:,:,i) = readtiff([nc82RT [num2str(nf(i)) '.tif']]);
end
toc

fprintf('initializing new volume...'); tic
imDAN = zeros(size(im(:,:,1)), 'uint16');
imNONDAN = zeros(size(im(:,:,1)), 'uint16');
DM = zeros(size(im), 'logical');
toc

if sum(sum(im(:,:,6))) > 0
    
    fprintf('loading nc82 associated DAN index, and corrdinates...'); tic
    load([matRT fn_mat], 'allData', 'allDANidx');
    allData = allData(allDANidx,:);
    
    zidx = allData(:,3) > (ex-5) & allData(:,3) < (ex+5);
    DM(sub2ind(size(im),allData(zidx,2),allData(zidx,1),allData(zidx,3)-ex+5)) = 1;
    DM = imdilate(DM,se);
    imDAN(DM(:,:,6)) = im(DM(:,:,6));
    imNONDAN(~DM(:,:,6)) = im(~DM(:,:,6));
    toc
    
end

fprintf('writing DAN / nonDAN associated nc82 puncta...'); tic
writetiff(uint16(imDAN), [nc82DANsaveRT fn_nc82]);
writetiff(uint16(imNONDAN), [nc82nonDANsaveRT fn_nc82]);
toc