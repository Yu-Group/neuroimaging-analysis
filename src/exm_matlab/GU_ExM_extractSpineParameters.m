function GU_ExM_extractSpineParameters(volPath, varargin)

% identify the spine using skeletonization
% input: volume path
% Gokul Upadhyayula, Oct, 2017


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('volPath', @isstr);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('zAniso', 1.8557, @isnumeric);
ip.addParameter('PixelSize', 0.097, @isnumeric); % in um
ip.addParameter('ClosingSphereSize', 5, @isnumeric); % to make spheres from the local max
ip.addParameter('PadSize', [3 3 3], @isvector);
ip.addParameter('MinBranchLength', 10, @isnumeric); % min length of branches
ip.addParameter('RemoveOverlap', 25, @isnumeric); % # of pixels at the border
ip.addParameter('CropScalingFactor', 5, @isnumeric); % cropping spines
ip.addParameter('CSFPercentile', 50, @isnumeric); % cropping spines
ip.addParameter('nPtsTrim', 0, @isnumeric); % number of points to trim
% image filtering options
ip.addParameter('MinVoxelVolume', 2000, @isnumeric); % set to zere if no "pre-cleaning" is necessary to remove high-freq noise (168 comes from sig = 2.3)
ip.addParameter('Verbose', true, @islogical); % to make spheres from the local max
ip.addParameter('Figure', false, @islogical); % to make spheres from the local max
ip.parse(volPath, varargin{:});
p = ip.Results;

% generate directories
[rt, fn, ~] = fileparts(volPath);
mkdir ([rt filesep 'DistanceMap']);
mkdir ([rt filesep 'Skeleton']);
mkdir ([rt filesep 'DistMapSkel']);
mkdir ([rt filesep 'SkelGraph']);
mkdir ([rt filesep 'GapClosed']);
mkdir ([rt filesep 'SpineCrops_mask']);
mkdir ([rt filesep 'SpineCrops_data']);
mkdir ([rt filesep 'SkelCrop_data']);

CSF = p.CropScalingFactor;
CSFPercentile = p.CSFPercentile;
nTrimPts = p.nPtsTrim;
if ~exist([rt filesep 'SkelGraph' filesep fn '.mat'], 'file') || p.Overwrite
    
    if ~exist([rt filesep 'GapClosed' filesep fn '_Closed.tif'], 'file') || p.Overwrite
        % read volume
        im = readtiff(volPath);
        
        % clean volume
        CC = bwconncomp(logical(im), 26);
        csize = cellfun(@numel, CC.PixelIdxList);
        idx = csize>=p.MinVoxelVolume;
        CC.NumObjects = sum(idx);
        CC.PixelIdxList = CC.PixelIdxList(idx);
        imGL = labelmatrix(CC)~=0;
        im = uint16(im) .* uint16(imGL);
        clear imGL
        
        % 3d fill
        fprintf('Closing Gaps...')
        tic
        im3close = imclose(logical(im), binarySphere(p.ClosingSphereSize));
        
        toc
        writetiff(uint8(im3close), [rt filesep 'GapClosed' filesep fn '_Closed.tif']);
    else
        im3close = readtiff([rt filesep 'GapClosed' filesep fn '_Closed.tif']);
    end
    
    % mex function - skeletonization with padded edges
    if ~exist([rt filesep 'Skeleton' filesep fn '_SK.tif'], 'file') || p.Overwrite
        if p.PadSize > 0
            pcm3 = padarray(im3close, p.PadSize, 0, 'both');
        else
           pcm3 = im3close;
        end
        fprintf('Skeletonizing...')
        tic
        skel = skeleton3D(logical(pcm3));
        toc
        skel = skel(4:end-3, 4:end-3, 4:end-3);
        
        clear pcm3
        
        % Generate Graph, clean up and generate skel
        fprintf('Converting Skeleton to Graph, cleaning up and back to Skeleton...')
        tic
        [~,node,link] = Skel2Graph3D(skel,p.MinBranchLength);
        w = size(skel,1);
        l = size(skel,2);
        h = size(skel,3);
        skel = Graph2Skel3D(node,link,w,l,h);
        toc
        writetiff(uint8(skel), [rt filesep 'Skeleton' filesep fn '_SK.tif']);
    else
        skel = readtiff([rt filesep 'Skeleton' filesep fn '_SK.tif']);
    end
    
    % generate new graph
    [~,node,link] = Skel2Graph3D(skel,p.MinBranchLength);
    
    % distance map
    if ~exist([rt filesep 'DistanceMap' filesep fn '_DT.tif'], 'file') || p.Overwrite
        fprintf('Calculating Distance Map...')
        tic
        dtcm3 = bwdistsc(~im3close,[1 1 p.zAniso]);
        tic
        writetiff(single(dtcm3), [rt filesep 'DistanceMap' filesep fn '_DT.tif']);
    else
        dtcm3 = readtiff( [rt filesep 'DistanceMap' filesep fn '_DT.tif']);
    end
    
    %     clear im im3close
    
    % generate skel readout of the distancemap
    if ~exist([rt filesep 'DistMapSkel' filesep fn '_DTSK.tif'], 'file') || p.Overwrite
        fprintf('Generating Skeleton Distance Map...')
        tic
        diaSkel = zeros(size(dtcm3),'single');
        diaSkel(skel) = single(dtcm3(skel));
        toc
        writetiff(diaSkel, [rt filesep 'DistMapSkel' filesep fn '_DTSK.tif']);
    else
        diaSkel = readtiff([rt filesep 'DistMapSkel' filesep fn '_DTSK.tif']);
    end
    
    % write out link diameter
    for i = 1:numel(link)
        link(i).diameter = dtcm3(link(i).point);
    end
    
    % get coordinates
    syxz = size(skel);
    for i = 1:numel(link)
        [link(i).coord(:,1),link(i).coord(:,2),link(i).coord(:,3)] = ind2sub(syxz, [link(i).point]);
    end
    
    for i = 1:numel(node)
        [node(i).coord(:,1),node(i).coord(:,2),node(i).coord(:,3)] = ind2sub(syxz, [node(i).idx]);
    end
    
    clear skel dtcm3 diaSkel csize CC i idx
    
    save([rt filesep 'SkelGraph' filesep fn '.mat']);
else
    load([rt filesep 'SkelGraph' filesep fn '.mat'], 'link', 'node', 'syxz');
end

w = syxz(1);
l = syxz(2);
h = syxz(3);

% remove nodes and links at the volume boundary



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % total length of network
% wl = sum(cellfun('length',{node.links}));
% 
% skel2 = Graph2Skel3D(node,link,w,l,h);
% [~,node2,link2] = Skel2Graph3D(skel2,p.MinBranchLength);
% 
% % calculate new total length of network
% wl_new = sum(cellfun('length',{node2.links}));
% 
% % iterate the same steps until network length changed by less than 0.5%
% while(wl_new~=wl)
%     wl = wl_new;
%     skel2 = Graph2Skel3D(node2,link2,w,l,h);
%     [~,node2,link2] = Skel2Graph3D(skel2,0);
%     wl_new = sum(cellfun('length',{node2.links}));
% end
% 
% for i = 1:numel(link2)
%     [link2(i).coord(:,1),link2(i).coord(:,2),link2(i).coord(:,3)] = ind2sub(syxz, [link2(i).point]);
% end
% 
% for i = 1:numel(node2)
%     [node2(i).coord(:,1),node2(i).coord(:,2),node2(i).coord(:,3)] = ind2sub(syxz, [node2(i).idx]);
% end
% if ~exist('dtcm3', 'var')
%     dtcm3 = readtiff( [rt filesep 'DistanceMap' filesep fn '_DT.tif']);
% end
% for i = 1:numel(link2)
%     link2(i).diameter = dtcm3(link2(i).point)*2;
% end
% 
% %% calculate end note parameters
% 
% %check if node is at the overlap region
% olx = [1:p.RemoveOverlap, syxz(2)-p.RemoveOverlap:syxz(2)];
% oly = [1:p.RemoveOverlap, syxz(1)-p.RemoveOverlap:syxz(1)];
% olz = [1:p.RemoveOverlap, syxz(3)-p.RemoveOverlap:syxz(3)];
% for i=1:length(node2)
%     if ismember(node2(i).comx, olx) || ismember(node2(i).comy, oly) || ismember(node2(i).comz, olz)
%         node2(i).RemoveEdgeNode = 1;
%     else
%         node2(i).RemoveEdgeNode = 0;
%     end
% end
% 
% % generate cell arrays of end links that are not at the volume edge
% en = 1;
% bn = 1;
% for i=1:length(node2)
%     for j=1:length(node2(i).links)
%         d = link2(node2(i).links(j)).diameter;
%         if node2(i).RemoveEdgeNode == 0
%             if(node2(link2(node2(i).links(j)).n2).ep==1)
%                 allEndLinkDia{en} = d;
%                 en = en+1;
%             else
%                 allMiddleLinkDia{bn} = d;
%                 bn = bn+1;
%             end
%         end
%     end
% end

%% generate masked region of all identified spines



%% generate cropped region of all identified spines
% 
% if ~exist('im', 'var')
%     im = readtiff(volPath);
% end
% % diaSkel = readtiff([rt filesep 'DistMapSkel' filesep fn '_DTSK.tif']);
% en = 1;
% for i=1:length(node2)
%     for j=1:length(node2(i).links)
%         if(node2(link2(node2(i).links(j)).n2).ep==1) && node2(i).RemoveEdgeNode == 0
%             SF = prctile(link2(node2(i).links(j)).diameter*CSF, CSFPercentile);
%             
%             ys = (min(link2(node2(i).links(j)).coord(1:end-nTrimPts,1))-SF):(max(link2(node2(i).links(j)).coord(1:end-nTrimPts,1))+ SF)+1;
%             xs = (min(link2(node2(i).links(j)).coord(1:end-nTrimPts,2))-SF):(max(link2(node2(i).links(j)).coord(1:end-nTrimPts,2))+ SF)+1;
%             zs = (min(link2(node2(i).links(j)).coord(1:end-nTrimPts,3))-SF):(max(link2(node2(i).links(j)).coord(1:end-nTrimPts,3))+ SF)+1;
%             
%             %             ys = (min(link2(node2(i).links(j)).coord(nTrimPts:end,1))-SF):(max(link2(node2(i).links(j)).coord(nTrimPts:end,1))+ SF)+1;
%             %             xs = (min(link2(node2(i).links(j)).coord(nTrimPts:end,2))-SF):(max(link2(node2(i).links(j)).coord(nTrimPts:end,2))+ SF)+1;
%             %             zs = (min(link2(node2(i).links(j)).coord(nTrimPts:end,3))-SF):(max(link2(node2(i).links(j)).coord(n%TrimPts:end,3))+ SF)+1;
%             
%             % spine mask
%             if numel(uint16(ys(ys>=1)))>10 && numel(uint16(xs(xs>=1)))>10 && numel(uint16(zs(zs>=1)))>10
%                 mask = logical(dtcm3(uint16(ys(ys>=1 & ys<syxz(1))),uint16(xs(xs>=1 & xs<syxz(2))),uint16(zs(zs>=1 & zs<syxz(3)))));
%                 SpineCrop = im(uint16(ys(ys>=1 & ys<syxz(1))),uint16(xs(xs>=1 & xs<syxz(2))),uint16(zs(zs>=1 & zs<syxz(3))));
%                 SkelCrop = zeros(size(dtcm3),'single');
%                 SkelCrop(link2(node2(i).links(j)).point) = link2(node2(i).links(j)).diameter;
%                 SkelCrop = SkelCrop(uint16(ys(ys>=1 & ys<syxz(1))),uint16(xs(xs>=1 & xs<syxz(2))),uint16(zs(zs>=1 & zs<syxz(3))));
%                 idx = find(mask);
%                 [ny,nx,nz] = size(mask);
%                 [yi,xi,zi] = ind2sub([ny,nx,nz], idx);
%                 b = 5; % boundary pixels for cropping
%                 xa = max(min(xi)-b,1):min(max(xi)+b,nx);
%                 ya = max(min(yi)-b,1):min(max(yi)+b,ny);
%                 za = max(min(zi)-b,1):min(max(zi)+b,nz);
%                 mask = mask(ya,xa,za);
%                 SpineCrop = SpineCrop(ya,xa,za);
%                 SkelCrop = SkelCrop(ya,xa,za);
%                 
%                 % interpolate
%                 [ny,nx,nz] = size(mask);
%                 [y,x,z] = ndgrid(1:ny,1:nx,1:nz);
%                 [Y,X,Z] = ndgrid(1:ny,1:nx,1:1/p.zAniso:nz);
%                 imask = interp3(x,y,z,double(mask),X,Y,Z,'nearest');
%                 iSpineCrop = interp3(x,y,z,double(SpineCrop),X,Y,Z,'linear');
%                 iSkelCrop = interp3(x,y,z,double(SkelCrop),X,Y,Z,'linear');
%                 
%                 % write out cropped spines
%                 writetiff(uint8(imask), [rt filesep 'SpineCrops_mask' filesep fn '_n1_' num2str(link2(node2(i).links(j)).n1) '_n2_' num2str(link2(node2(i).links(j)).n2) '.tif']);
%                 writetiff(uint16(iSpineCrop), [rt filesep 'SpineCrops_data' filesep fn '_n1_' num2str(link2(node2(i).links(j)).n1) '_n2_' num2str(link2(node2(i).links(j)).n2) '.tif']);
%                 writetiff(uint16(iSkelCrop), [rt filesep 'SkelCrop_data' filesep fn '_n1_' num2str(link2(node2(i).links(j)).n1) '_n2_' num2str(link2(node2(i).links(j)).n2) '.tif']);
%                 
%                 % calc spine mip
%                 %                 MIPimask = max(imask,[],3);
%                 %                 MIPimask = squeeze(MIPimask(:,:,1));
%                 %                 writetiff(uint8(MIPimask), [rt filesep 'SpineCrops_mask' filesep 'zMIP' filesep fn '_n1_' num2str(link2(node2(i).links(j)).n1) '_n2_' num2str(link2(node2(i).links(j)).n2) '.tif']);
%                 %                 MIPiSpineCrop = max(iSpineCrop,[],3);
%                 %                 MIPiSpineCrop = squeeze(MIPiSpineCrop(:,:,1));
%                 %                 writetiff(uint16(MIPiSpineCrop), [rt filesep 'SpineCrops_data' filesep 'zMIP' filesep fn '_n1_' num2str(link2(node2(i).links(j)).n1) '_n2_' num2str(link2(node2(i).links(j)).n2) '.tif']);
%                 en = en+1;
%                 node2(i).croppedFile = 1;
%             end
%         end
%     end
% end
% 
% clear dtcm3
% save([rt filesep 'SkelGraph' filesep fn '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% % generate figure
% if p.Figure
%
%     figure();
%     hold on;
%     [ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
%     for i=1:length(node2)
%         x1 = node2(i).comx;
%         y1 = node2(i).comy;
%         z1 = node2(i).comz;
%
%         if(node2(i).ep==1)
%             ncol = ceB;
%         else
%             ncol = ceR;
%         end
%
%         for j=1:length(node2(i).links)    % draw all connections of each node
%             if(node2(link2(node2(i).links(j)).n2).ep==1)
%                 col= cfB; % branches are blue
%             else
%                 col= cfR; % links are red
%             end
%             if(node2(link2(i).n1).ep==1)
%                 col= cfK;
%             end
%
%
%             % draw edges as lines using voxel positions
%             for k=1:length(link2(node2(i).links(j)).point)-1
%                 x3 = link2(node2(i).links(j)).coord(k,2);
%                 x2 = link2(node2(i).links(j)).coord(k+1,2);
%                 y3 = link2(node2(i).links(j)).coord(k,1);
%                 y2 = link2(node2(i).links(j)).coord(k+1,1);
%                 z3 = link2(node2(i).links(j)).coord(k,3);
%                 z2 = link2(node2(i).links(j)).coord(k+1,3);
%
%                 %                 [x3,y3,z3]=ind2sub([w,l,h],link(node(i).links(j)).point(k));
%                 %                 [x2,y2,z2]=ind2sub([w,l,h],link(node(i).links(j)).point(k+1));
%                 line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
%             end
%         end
%
%         % draw all nodes as yellow circles
%         plot3(y1,x1,z1,'o','Markersize',5,...
%             'MarkerFaceColor',ncol);
%
%     end
%     axis equal%;axis off;
%     set(gcf,'Color','white');
%     drawnow;
%     %     view(-17,46);
%     view(3)
%
%     ha = setupFigure(1,2, 'AxesWidth', 4, 'AxesHeight', 2,'SameAxes', false,...
%         'XSpace', [1.5 0.85 0.5], 'YSpace', [1 1 1]);
%     set(ha, 'FontSize', 6);
%
%     for i=1:length(node2)
%         for j=1:length(node2(i).links)
%             d = link2(node2(i).links(j)).diameter;
%             if(node2(link2(node2(i).links(j)).n2).ep==1)
%                 axes(ha(1))
%                 plot(1:numel(d),d, 'Color', cfG)
%                 %                 hold off
%             else
%                 axes(ha(2))
%                 plot(1:numel(d),d, 'Color', ceG)
%                 %                 hold off
%             end
%         end
%     end
% end
