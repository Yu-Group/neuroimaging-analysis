% GU_ExM_MouseVariableHomerStainingWorkspace

rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/MouseSomatoData/VariableHomerStaining2Neurons';
cd(rt)
load('seg48_HomerInfo.mat')
load('Set32_HomerInfo.mat')
load('Set32_YFPPerimCorrdinates.mat')
load('Set48_YFPPerimCorrdinates.mat')

%%
px = .114;
py = .097;
pz = .155;

% calculate distance of Homer to Neuron #32
HomerAll32(:,1) = HomerAll32_rp3.Centroid(:,1)*px;
HomerAll32(:,2) = HomerAll32_rp3.Centroid(:,2)*py;
HomerAll32(:,3) = HomerAll32_rp3.Centroid(:,3)*pz;

HomerAss32(:,1) = HomerAss32_rp3.Centroid(:,1)*px;
HomerAss32(:,2) = HomerAss32_rp3.Centroid(:,2)*py;
HomerAss32(:,3) = HomerAss32_rp3.Centroid(:,3)*pz;

YFP32_perimcoordXYZ(:,1) = YFP32_perimcoordXYZ(:,1)*px;
YFP32_perimcoordXYZ(:,2) = YFP32_perimcoordXYZ(:,2)*py;
YFP32_perimcoordXYZ(:,3) = YFP32_perimcoordXYZ(:,3)*pz;


tic
Mdl = KDTreeSearcher(YFP32_perimcoordXYZ); % all nonDAN synapse positions
toc

tic
[~,D_HomerAll32] = knnsearch(Mdl, HomerAll32, 'k',1);
toc

tic
[~,D_HomerAss32] = knnsearch(Mdl, HomerAss32, 'k',1);
toc

figure,
histogram(D_HomerAll32, 0:0.05:5);
hold on
histogram(D_HomerAss32, 0:0.05:5);

% calculate Homer Distance to Neuron #48
HomerAll48(:,1) = HomerAll48_rp3.Centroid(:,1)*px;
HomerAll48(:,2) = HomerAll48_rp3.Centroid(:,2)*py;
HomerAll48(:,3) = HomerAll48_rp3.Centroid(:,3)*pz;

HomerAss48(:,1) = HomerAss48_rp3.Centroid(:,1)*px;
HomerAss48(:,2) = HomerAss48_rp3.Centroid(:,2)*py;
HomerAss48(:,3) = HomerAss48_rp3.Centroid(:,3)*pz;

YFP48_perimcoordXYZ(:,1) = YFP48_perimcoordXYZ(:,1)*px;
YFP48_perimcoordXYZ(:,2) = YFP48_perimcoordXYZ(:,2)*py;
YFP48_perimcoordXYZ(:,3) = YFP48_perimcoordXYZ(:,3)*pz;


tic
Mdl = KDTreeSearcher(YFP48_perimcoordXYZ); % all nonDAN synapse positions
toc

tic
[~,D_HomerAll48] = knnsearch(Mdl, HomerAll48, 'k',1);
toc

tic
[~,D_HomerAss48] = knnsearch(Mdl, HomerAss48, 'k',1);
toc

figure,
histogram(D_HomerAll48, 0:0.05:5);
hold on
histogram(D_HomerAss48, 0:0.05:5);
%% calculate surface density of the homer for the two

xAniso = px/py;
zAniso = pz/py;

SA_32 = numel()
