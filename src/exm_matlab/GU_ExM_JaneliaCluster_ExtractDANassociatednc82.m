function GU_ExM_JaneliaCluster_ExtractDANassociatednc82(ex)

% Gokul Upadhyayula, Nov 2017
se = strel('sphere', 8);
ex = str2double(ex);
fprintf('calculating file names...'); tic
% get file names
mx = 15055;
my = 27964;
mz = 6647;
OL = 50;
s = repmat(750, [1,3]);

nx = ceil(mx/(s(1)-OL));
ny = ceil(my/(s(2)-OL));
nz = ceil(mz/(s(3)-OL));

% generate x y z cropping coordinates
clear xmin xmax ymin zmin ymax zmax

xmin = zeros(nx*ny*nz,1);
xmax = zeros(nx*ny*nz,1);
ymin = zeros(nx*ny*nz,1);
ymax = zeros(nx*ny*nz,1);
zmin = zeros(nx*ny*nz,1);
zmax = zeros(nx*ny*nz,1);

nn = 1;
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            xmin(nn) = (1+((i-1)*(s(1)-OL)));
            xmax(nn) = (s(1)*i-(OL*(i-1)));
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            ymax(nn) = (s(2)*j-(OL*(j-1)));
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            zmax(nn) = (s(3)*k-(OL*(k-1)));
            nn = nn + 1;
        end
    end
end


matRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch0/Analysis/DataStructures_Recalc/';
nc82RT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean/';
nc82DANsaveRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_DANAssociated/';
nc82nonDANsaveRT = '/groups/betzig/betziglab/4Stephan/171016_Flybrainsynapse/171016_100nmcluster/ch1/Analysis/clean_nonDANAssociated/';
fn_mat = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '_clean_densityDataRC.mat'];
fn_nc82 = [num2str(xmin(ex)) ',' num2str(ymin(ex)) ',' num2str(zmin(ex)) '_' num2str(xmax(ex)) ',' num2str(ymax(ex)) ',' num2str(zmax(ex)) '_clean.tif'];
toc

% read tiff
fprintf('reading volume...'); tic
im = readtiff([nc82RT fn_nc82]);
toc

fprintf('initializing new volume...'); tic
imDAN = zeros(size(im), 'uint16');
imNONDAN = zeros(size(im), 'uint16');
DM = zeros(size(im), 'logical');
toc

if sum(im(:) > 0)
    variableInfo = who('-file', [matRT fn_mat]);
    PF = ismember('DANSynapseIdx', variableInfo) & ismember('filteredData', variableInfo);
    if PF
        fprintf('loading nc82 associated DAN index, and corrdinates...'); tic
        load([matRT fn_mat], 'filteredData', 'DANSynapseIdx');
        filteredData = filteredData(DANSynapseIdx,:);
        toc
        
%         fprintf('geneating labels...'); tic
%         lholes = bwlabeln(logical(im));
%         toc
%         
%         fprintf('looping through coordinates...'); tic
%         danPpuncatIdx = zeros(numel(filteredData(:,1)),1);
%         for n = 1:numel(filteredData(:,1))
%             danPpuncatIdx(n) = lholes(filteredData(n,2),filteredData(n,1),filteredData(n,3));
%         end
%         idx = ismember(lholes,danPpuncatIdx);
%         imDAN(idx) = im(idx);
%         imNONDAN(~idx) = im(~idx);
            DM(filteredData(:,2),filteredData(:,1),filteredData(:,3)) = 1;
            DM = imdilate(DM,se);
            imDAN(DM) = im(DM);
            imNONDAN(~DM) = im(~DM);
        toc
    end
end

fprintf('writing DAN / nonDAN associated nc82 puncta...'); tic
writetiff(uint16(imDAN), [nc82DANsaveRT fn_nc82]);
writetiff(uint16(imNONDAN), [nc82nonDANsaveRT fn_nc82]);
toc