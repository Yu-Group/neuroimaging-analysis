rt = 'D:\Gokul\Mousebrainsynapse\Ex5_0p18\Analysis\1008\';
fn = 'Position-_channel 0 [3458, 39079, 5821]_clean.tif';
sfn = 'filledMask.tif';

im = readtiff([rt fn]);
cm3 = readtiff([rt 'filled_3d.tif']);
skel = readtiff([rt 'Skel_filled_3d.tif']);
diaSkel =  readtiff([rt 'DistanceSkel.tif']);
%%
cm = logical(im);
HSIZE = 5;
% 2D fill
for z = 1:size(cm,3)
    cm(:,:,z) = imfill(cm(:,:,z), 8, 'holes');
end
for x = 1:size(cm,2)
    cm(:,x,:) = imfill(squeeze(cm(:,x,:)), 8, 'holes');
end
for y = 1:size(cm,1)
    cm(y,:,:) = imfill(squeeze(cm(y,:,:)), 8, 'holes');
end
for z = 1:size(cm,3)
    cm(:,:,z) = imfill(cm(:,:,z), 8, 'holes');
end
for x = 1:size(cm,2)
    cm(:,x,:) = imfill(squeeze(cm(:,x,:)), 8, 'holes');
end
for y = 1:size(cm,1)
    cm(y,:,:) = imfill(squeeze(cm(y,:,:)), 8, 'holes');
end
writetiff(uint8(cm), [rt 'filled_2d.tif']);
% 3d fill
cm3 = imclose(logical(im), binarySphere(HSIZE));
writetiff(uint8(cm3), [rt 'filled_3d.tif']);

%%

SE = strel('sphere',3);

ecm3 = imerode(cm3, SE);
writetiff(uint8(ecm3), [rt 'Erodfilled_3d.tif']);

%% mex function
%pad edges
pcm3 = padarray(cm3,[3 3 3],0,'both');
tic
skel = skeleton3D(logical(pcm3));
toc
skel = skel(4:end-3,4:end-3,4:end-3);
writetiff(uint8(skel), [rt 'Skel_filled_3d.tif']);
%% non mex
% tic 
% skel2 = Skeleton3D(cm3);
% toc
% writetiff(uint8(skel2), [rt 'Skel2_filled_3d.tif']);
%%
tic
 [A,node,link] = Skel2Graph3D(skel,10);
 toc
 
 %%
 zAniso = .18/0.097;
 dtcm3 = bwdistsc(~cm3,[1 1 zAniso]);
 writetiff(single(dtcm3), [rt 'DistanceTransfrom.tif']);
 %%
 diaSkel = zeros(size(dtcm3),'uint16');
diaSkel(skel) = uint16(dtcm3(skel));
writetiff(diaSkel, [rt 'DistanceSkel.tif']);

%% generate a cleaned up graph
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
skel2 = Graph2Skel3D(node,link,w,l,h);
writetiff(uint8(skel2), [rt 'Skel_fromGraph.tif']);

tic
 [A,node,link] = Skel2Graph3D(skel2,10);
 toc
 
 %%
 for i = 1:numel(link)
    link(i).diameter = dtcm3(link(i).point); 
 end
 
 %%
 for i = 1:numel(link)
 link(i).ep = 0;
 end
 idx = find([node.ep]);
 endLinks = [node(idx).links];
clear endlinkidx
allEndLinkDia = [];
allMiddleLinkDia = [];
 for i = endLinks
  link(i).ep = 1;
  allEndLinkDia = [allEndLinkDia min(link(i).diameter)];
 end
 
 
  idx = find([link.ep] ==0);
allMiddleLinkDia = [];
 for i = idx
  allMiddleLinkDia = [allMiddleLinkDia min(link(i).diameter)];
 end
 
 figure
histogram(allEndLinkDia*.0268*2)
hold on
histogram(allMiddleLinkDia*.0268*2)

%% get coordinates
syxz = size(skel);
for i = 1:numel(link)
   [link(i).coord(:,1),link(i).coord(:,2),link(i).coord(:,3)] = ind2sub(syxz, [link(i).point]); 
end

for i = 1:numel(node)
     [node(i).coord(:,1),node(i).coord(:,2),node(i).coord(:,3)] = ind2sub(syxz, [node(i).idx]); 
end