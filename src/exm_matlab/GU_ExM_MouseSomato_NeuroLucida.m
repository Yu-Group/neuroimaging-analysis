% GU_ExM_MouseSomato_NeuroLucida processing

rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/Neurolucida/Ex9_Spine_Homer_collocalization/';

p = recursiveDir(rt, 'maxlevel', 1);

for i = 28%[3:27,29:numel(p)]
    fnlist = dir([p{i} '*.tif']);
    fnlist = {fnlist.name}';
    
    %             GU_ExM_calcPunctaDensityVol(fnlist{3}, fnlist{1},...
    %             'FindLocalMaxima', false, 'MinThreshold', [1450,2200],'Verbose', true,'OTSUGaussKernel', 0.5,...
    %             'MinVoxelVolume', [246, 1008],'ExpansionFactor', 3.62); % sigma 3.9
    %     GU_ExM_Janelia_calcPunctaDensityVol_neurolucida
    submit = ['bsub -n 4 -R"affinity[core(1)]" -J' ' "HoCalc' num2str(i) '" -o ~/logs/HoCalc' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_ExM_Janelia_calcPunctaDensityVol_neurolucida ' [p{i} fnlist{3}] ' ' [p{i} fnlist{1}] ''''];
    [stat, res] = system(submit,'-echo');
end

%%
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/Neurolucida/Ex9_Spine_Homer_collocalization/';

load([rt 'coordinateOffsets.mat'])
rtcoord = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/Neurolucida/Ex9_Spine_Homer_collocalization/Newcoordinates/';
for i = 1:numel(Filename)
    t = readtable([rtcoord Filename{i} '.xlsx']);
    data(i).x = t.x(~isnan(t.x));
    data(i).y = t.y(~isnan(t.x));
    data(i).z = t.z(~isnan(t.x));
    if sum(strcmp('HeadDiameter__m_',t.Properties.VariableNames))
        data(i).HeadDiameter = t.HeadDiameter__m_(~isnan(t.x));
    else
        data(i).HeadDiameter = t.HeadDiameter_m_(~isnan(t.x));
    end
    data(i).xoffset = xcordinatesum(i);
    data(i).yoffset = ycordinatesum(i);
    data(i).zoffset = zcordinatesum(i);
    data(i).BackboneLength = t.BackboneLength__m_(~isnan(t.x));
    data(i).ExMNeckDiameterm = t.ExMNeckDiameter__m_(~isnan(t.x));
    data(i).SpineType =t.SpineType(~isnan(t.x));
    data(i).Tree = t.Tree(~isnan(t.x));
    
end

%% get pixel positions
px = 0.0267;
pz = 0.0503;
zAniso = 0.18/0.097;
for i = 1:numel(Filename)
    data(i).xPix = abs(data(i).x - data(i).xoffset)/px;
    data(i).yPix = abs(data(i).y - data(i).yoffset)/px;
    data(i).zPix = abs(data(i).z - data(i).zoffset)/pz*zAniso;
    data(i).HeadDiameterPix = data(i).HeadDiameter/px;
end


%% map the distances

for i = 1:numel(Filename)
    frt = [p{i+3} 'Analysis' filesep '246' filesep 'DataStructures' filesep];
    fn = dir([frt '*.mat']);
    fn = fn.name;
    d = load([frt fn], 'filteredData');
    d = d.filteredData;
    d(:,3) = d(:,3)*zAniso;
    Mdl = KDTreeSearcher(d);
    for j = 1:numel(data(i).xPix)
        IdxKDT = rangesearch(Mdl,horzcat(data(i).xPix(j),data(i).yPix(j),data(i).zPix(j)),data(i).HeadDiameterPix(j)/2);% Search radius = bandwidth2;
        data(i).nHomerHeadradius(j) = numel(IdxKDT{:});
        IdxKDT = rangesearch(Mdl,horzcat(data(i).xPix(j),data(i).yPix(j),data(i).zPix(j)),0.5/px);% Search radius = bandwidth2;
        data(i).nHomerHalfMicron(j) = numel(IdxKDT{:});
    end
    data(i).SpineHomerDist = knnsearch(Mdl, horzcat(data(i).xPix,data(i).yPix,data(i).zPix), 'k',1);
end
