 
rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/Mitoclustersegmentation_test/171129_YFPsample4/Ex1_nosoma_0p18/';
cd(rt)
fn = dir([rt '*.tif']);
fn = {fn.name};

px = 0.112;
py = 0.097;
pz = 0.18;
%%
GU_ExM_calcPunctaDensityVol_1ch([rt fn{2}], ...
        'FindLocalMaxima', true, 'Verbose', true,'OTSUGaussKernel', 0.5,...
        'OTSUMaxPer', 99,'MinVoxelVolume', [246],'ExpansionFactor', 3.99,...
        'minAZdia', 0)
    
    %%
    load('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/Mitoclustersegmentation_test/171129_YFPsample4/Ex1_nosoma_0p18/Analysis_1ch/246/DataStructures/ch1 [9100, 4090, 2510]_Ex1_nosomaneurites_246_densityData.mat')
    %%
    
    figure,
    
    scatter
    
    
    %%
    imcX = filteredData(:,1) * px;
    imcY = filteredData(:,2) * py;
    imcZ = filteredData(:,3) * pz;


    %%
    clusterSizeMin = 20;
    bandwidth = 2;
    % first round of clustering to remove duplicate detections
    fprintf('clustering to remove duplicates (diameter pre expanded ~ 200)*expansion factor...')
    tic
    ptData(:,1) = imcX;
    ptData(:,2) = imcY;
    ptData(:,3) = imcZ;
    [clusterInfo, pointToClusterMap, pointTraj] = MeanShiftClustering( ptData, bandwidth);
    [idx2] = find([clusterInfo.numPoints] > clusterSizeMin);
    toc
    %%
    [idx2] = find([clusterInfo.numPoints] > clusterSizeMin);
figure

% extract the cluster points
if numel(idx2) > 0
for i = 1:numel(idx2)
    clus(i).x2= imcX(clusterInfo(idx2(i)).ptIdData);
    clus(i).y2= imcY(clusterInfo(idx2(i)).ptIdData);
    clus(i).z2= imcZ(clusterInfo(idx2(i)).ptIdData);
    clus(i).center = clusterInfo(idx2(i)).ptClusterCenter;
    scatter3(clus(i).x2,clus(i).y2,clus(i).z2,'.')
    hold on
    pause
end
end  
    
    ClusCenters = cat(1, clus.center);
    scatter3(ClusCenters(:,1),ClusCenters(:,2),ClusCenters(:,3),50),hold on;
    figure, 
    scatter3(imcX, imcY, imcZ, '.') 
    hold on
    scatter3()
    
    %%
    % second round of clustering to group
    fprintf('clustering to calculate spot density (diameter -- 1um3)*expansion factor^3...')
    tic
    ClusCenters = cat(1, clusterInfo.ptClusterCenter);
    [clusterInfo2, ~, ~] = MeanShiftClustering( ClusCenters, bandwidth2);
    toc
    
    ClusCenters2 = cat(1, clusterInfo2.ptClusterCenter);
    cm = jet(numel([clusterInfo2.numPoints]));
    cm = flip(cm,2);
    [~,idx] = sort([clusterInfo2.numPoints], 'descend');
    %%
       s  = regionprops(im,'centroid');
    imcentroids = cat(1, s.Centroid);
    imcX = imcentroids(:,1);
    imcY = imcentroids(:,2);
    imcZ = imcentroids(:,3);
    toc
    
    % first round of clustering to remove duplicate detections
    fprintf('clustering to remove duplicates (diameter pre expanded ~ 200)*expansion factor...')
    tic
    ptData(:,1) = imcX;
    ptData(:,2) = imcY;
    ptData(:,3) = imcZ;
    [clusterInfo, ~, ~] = MeanShiftClustering( ptData, bandwidth);
    toc
    
    % second round of clustering to group
    fprintf('clustering to calculate spot density (diameter -- 1um3)*expansion factor^3...')
    tic
    ClusCenters = cat(1, clusterInfo.ptClusterCenter);
    [clusterInfo2, ~, ~] = MeanShiftClustering( ClusCenters, bandwidth2);
    toc
    
    ClusCenters2 = cat(1, clusterInfo2.ptClusterCenter);
    cm = jet(numel([clusterInfo2.numPoints]));
    cm = flip(cm,2);
    [~,idx] = sort([clusterInfo2.numPoints], 'descend');