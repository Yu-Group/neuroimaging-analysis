% % GU_ExMFigures_HomerBassoonSpines
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/HomerBassoon_Spine/dataforFigures/';
cd(rt)
load('20180529_HomerBassson_rp.mat')
ef = 4.04;
px = .097; %in um
pz = 0.25;
zAniso = pz/px;
%% plot distances - centroid or weighted centroid

% % H = H_wc;
% % B = B_wc;
% % HY = HY_wc;
% % BY = BY_wc;

%%%%%%%%%%%%%%%%%%%% centroids %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homer_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),H_rp)]*px/ef;
% Homer_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),H_rp)]*px/ef;
% Homer_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),H_rp)]*zAniso*px/ef;
%
% Bassoon_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),B_rp)]*px/ef;
% Bassoon_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),B_rp)]*px/ef;
% Bassoon_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),B_rp)]*zAniso*px/ef;
%
% HomerYFP_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),HY_rp)]*px/ef;
% HomerYFP_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),HY_rp)]*px/ef;
% HomerYFP_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),HY_rp)]*zAniso*px/ef;
%
% BassoonYFP_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),BY_rp)]*px/ef;
% BassoonYFP_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),BY_rp)]*px/ef;
% BassoonYFP_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),BY_rp)]*zAniso*px/ef;

%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%% weighted centroids %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Homer_centroid(:,1) = H_wc.WeightedCentroid(:,1)*px/ef;
Homer_centroid(:,2) = H_wc.WeightedCentroid(:,2)*px/ef;
Homer_centroid(:,3) = H_wc.WeightedCentroid(:,3)*zAniso*px/ef;

Bassoon_centroid(:,1) = B_wc.WeightedCentroid(:,1)*px/ef;
Bassoon_centroid(:,2) = B_wc.WeightedCentroid(:,2)*px/ef;
Bassoon_centroid(:,3) = B_wc.WeightedCentroid(:,3)*zAniso*px/ef;

HomerYFP_centroid(:,1) = HY_wc.WeightedCentroid(:,1)*px/ef;
HomerYFP_centroid(:,2) = HY_wc.WeightedCentroid(:,2)*px/ef;
HomerYFP_centroid(:,3) = HY_wc.WeightedCentroid(:,3)*zAniso*px/ef;

BassoonYFP_centroid(:,1) = BY_wc.WeightedCentroid(:,1)*px/ef;
BassoonYFP_centroid(:,2) = BY_wc.WeightedCentroid(:,2)*px/ef;
BassoonYFP_centroid(:,3) = BY_wc.WeightedCentroid(:,3)*zAniso*px/ef;


tic
B_Mdl = KDTreeSearcher(Bassoon_centroid);
toc

tic
H_Mdl = KDTreeSearcher(Homer_centroid);
toc

tic
BY_Mdl = KDTreeSearcher(BassoonYFP_centroid);
toc

tic
HY_Mdl = KDTreeSearcher(HomerYFP_centroid);
toc

tic
[~,D_HB] = knnsearch(H_Mdl, Bassoon_centroid, 'k',1);
toc

tic
[~,D_YHB] = knnsearch(HY_Mdl, BassoonYFP_centroid, 'k',1);
toc

%
ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
histogram(D_HB*1000,0:10:600)
set(ha, 'FontSize', 6);
xlim([0 600])
xticks(0:150:600)
ylabel('# events');
xlabel('Distance (nm)');
legend(['n=' num2str(numel(D_HB))])
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'DistanceFromHomer2ClosestBasoon.eps']);
% [mu, sigma, ~, ~] = fitGaussianModeToCDF(D_HB, 'display',true)

ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
histogram(D_YHB*1000,0:10:600)
set(ha, 'FontSize', 6);
xlim([0 600])
xticks(0:150:600)
ylabel('# events');
xlabel('Distance (nm)');
legend(['n=' num2str(numel(D_YHB))])
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'DistanceFromYFPHomer2ClosestBasoon.eps']);
% [mu, sigma, ~, ~] = fitGaussianModeToCDF(D_YHB, 'display',true)

%% plot homer and bassoom major axes

% rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/HomerBassoon_Spine/';
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/HomerBassoon_Spine/';
cd(rt)
load('20180601_HomerBassson_rp3_PrincipleAxisLength.mat')

ef = 4.04;
px = .097; %in um
pz = 0.25;
zAniso = pz/px;

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
ms = 3;

fn{1} = [date 'HomerBassoon_MajorAxisScatter_ms' num2str(ms) '_majorAxis.eps'];
fn{3} = [date 'HomerBassoon_MajorAxisScatter_ms' num2str(ms) '_minorAxis.eps'];

for i = [1,3]
    
    Homer_PA_all(:,1) = iH_PA.PrincipalAxisLength(:,i)*px/ef;
    Bassoon_PA(:,1) = iB_PA.PrincipalAxisLength(:,i)*px/ef;
    Homer_PA = Homer_PA_all(idx);
    
    ha = setupFigure(1,1, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', false,...
        'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
    axes(ha(1))
    a = jet(256);
    params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};
    X = [(Homer_PA*1000)'; (Bassoon_PA*1000)'];
    densityScatter(X, params);
    colormap(jet)
    xlim([0 3000]); ylim([0 3000]);grid on
    % xticklabels([10^-1,10^0, 10, 10^2])
    % yticks(0:5:20)
    ylabel('Bassoon Major Axis (nm)');
    xlabel('Homer Major Axis (nm)');
    legend(['n=' num2str(numel(Homer_PA))])
    % title('Soma')
    f0 = gcf();
%     print(f0, '-painters','-depsc', '-loose',fn{i});
end
%% correlation

[rho,pval] = corr(Bassoon_PA,Homer_PA);

[MIC,NLR,MAS,MEV,MCN] = maxInformationCoef([Bassoon_PA';Homer_PA'])

%% plot distributions

Homer_PA_all(:,1) = iH_PA.PrincipalAxisLength(:,1)*px/ef;
Bassoon_PA(:,1) = iB_PA.PrincipalAxisLength(:,1)*px/ef;
Homer_PA = Homer_PA_all(idx);


ha = setupFigure(1,1, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', false,...
    'XSpace', [1.5 2 1.5], 'YSpace', [1.5 1.5 1.5]);
set(ha, 'FontSize', 6);

n = 1;
axes(ha(n))

GU_stairsHistPlot({(Homer_PA*1000), (Bassoon_PA*1000)}, 'ShowECDF', true, 'BinSize',10,'CreateFigure', false,...
    'PercentileEnd', 100, 'LineWidth',0.75)
ax = gca;
ax.BoxStyle = 'full';
box on
% title(ln{n})
legend(['n = ' num2str(numel(Homer_PA))], 'Homer', ['n = ' num2str(numel(Bassoon_PA))], 'Bassoon')
legend('boxoff')
set(gca,'fontsize',6, 'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k')
xlim([0 2000]);
yyaxis left
%     ylim([0 0.008]);
xticks(0:250:2000)
grid off

ylabel('Relative Frequency', 'color', 'k');
xlabel('Major Axis (nm)', 'color', 'k');
n = n+1;
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'HomerBassoonMajorAxis.eps']);

%%
close all
opts = {'display',true,'Bandwidth', 0.01,'MeanIntialAdjustmentFactor', 1.0, 'CreateFigure', false};
figure
[mu, sigma, ~, ~] = GU_fitGaussianModeToPDF(Bassoon_PA, opts{:})
figure
[mu, sigma, ~, ~] = GU_fitGaussianModeToPDF(Homer_PA, opts{:})
%%
% rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/HomerBassoon_Spine/';
% % rt = 'W:\Users\Gokul\Dropbox\Manuscript_ExLLSM\GokulWorkspace\HomerBassoon_Spine\';
% cd(rt)
%
% ef = 4.04;
% px = .097; %in um
% pz = 0.25;
% zAniso = pz/px;
%
% %% reclean bassoon and find pairs of homer and basson
%
% % imB = readtiff(['D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch1\ch1_1008_clean.tif']);
% % imB = imB-1500;
% % writetiff(imB, 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch1\ch1_1008_clean_v2.tif');
% imB = readtiff(['D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch1_v2\Analysis_1ch\246\ch1_1008_clean_v2_246_clean.tif']);
% imH = readtiff(['D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch2\ch2_1008_clean.tif']);
% %%
% %
% % tic
% % labelB = bwlabeln(logical(imB));
% % toc
% % logHom = logical(imH);
% % dilogHom = imdilate(logHom, strel('sphere', 3));
% % JuxBassIdx = unique(labelB(dilogHom));
% % JuxBassIdx = JuxBassIdx(JuxBassIdx>0);
% % JuxBass = zeros(size(imB),'logical');
% % JuxBass(ismember(labelB,JuxBassIdx)) = 1;
% %
% % JIB = zeros(size(imB), 'uint16');
% % JIB(JuxBass) = imB(JuxBass); %% bassoon that is juxaposed with the homer - for all
% %
% % labelH = bwlabeln(logical(imH));
% % dilogBas = imdilate(JuxBass, strel('sphere',3));
% % JuxHomerIdx = unique(labelH(dilogBas));
% % JuxHomer = zeros(size(imH),'logical');
% % JuxHomer(ismember(labelH, JuxHomerIdx)) = 1;
% %
% % JIH = zeros(size(imH), 'uint16');
% % JIH(JuxHomer) = imH(JuxHomer); %% bassoon that is juxaposed with the homer - for all
% %
% % writetiff(JIB, 'BassonJuxH.tif');
% % writetiff(JIH, 'HomerJuxB.tif');
% %
% % %%
% %
% % JIB = readtiff('BassonJuxH.tif');
% % JIH = readtiff('HomerJuxB.tif');
% %     tic
% % iimB = GU_interp3DlargeVolume(JIB, 'zAniso', 0.25/0.097, 'nzvols', 7);toc
% % writetiff(iimB, 'interBassonJuxH.tif');
% %
% % tic
% % iimH = GU_interp3DlargeVolume(JIH, 'zAniso', 0.25/0.097, 'nzvols', 7);toc
% % writetiff(iimH, 'interHomerJuxB.tif');
% % tic
% % B_rp3  = regionprops3(logical(iimB));toc
% % tic
% % H_rp3  = regionprops3(logical(iimH));toc
% % save('20180523_interpolated_HomerBassson_rp3.mat', 'B_rp3', 'H_rp3');
% ef = 4.04;
% px = .097; %in um
% pz = 0.25;
% zAniso = pz/px;
%
% rt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\';
%  cd(rt)
% load('20180523_interpolated_HomerBassson_rp3.mat')
% Homer_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),H_rp3)]*px/ef;
% Homer_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),H_rp3)]*px/ef;
% Homer_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),H_rp3)]*px/ef;
%
% Bassoon_centroid(:,1) = [arrayfun(@(x) x.Centroid(1),B_rp3)]*px/ef;
% Bassoon_centroid(:,2) = [arrayfun(@(x) x.Centroid(2),B_rp3)]*px/ef;
% Bassoon_centroid(:,3) = [arrayfun(@(x) x.Centroid(3),B_rp3)]*px/ef;
%
% % yfp_h_idx = zeros(1, numel(Homer_centroid(:,1)));
% % for i = 1:numel(Homer_centroid(:,1))
% % yfp_h_idx(i) = im(round(Homer_centroid(i,2)/px*ef),round(Homer_centroid(i,1)/px*ef),round(Homer_centroid(i,3)/px*ef/zAniso));
% % end
%
% tic
% B_Mdl = KDTreeSearcher(Bassoon_centroid);
% toc
%
% tic
% H_Mdl = KDTreeSearcher(Homer_centroid);
% toc
%
% load('20180525_YFP_homer_idx.mat')
% tic
% H_Mdl_yfp = KDTreeSearcher(Homer_centroid(find(yfp_h_idx),:));
% toc
% %%
% load('Bassoon_YFP_ch0_1008_densityData.mat', 'Vol2PunctaIdx', 'filteredData', 'filteredsynapseDensity');
% B_filteredData = filteredData;
% B_YFPidx = Vol2PunctaIdx;
% B_filteredsynapseDensity = filteredsynapseDensity;
%
% load('Homer_YFP_ch0_1008_densityData.mat', 'Vol2PunctaIdx', 'filteredData', 'filteredsynapseDensity');
% H_filteredData = filteredData;
% H_YFPidx = Vol2PunctaIdx;
% H_filteredsynapseDensity= filteredsynapseDensity;
%
% V = 1;% in um3 surveyed volume
% r = (V*3/4/pi)^(1/3); % in um to give a Vum3 spherical volume
% re = r*ef; % radius in expanded samples; in um
% bandwidth2 = re/px;
%
% %%
% % % read im first
% % im = readtiff('ch0_1008_clean.tif');
% % pim = bwperim(logical(im));
% % [iy,ix,iz] = ind2sub(size(pim), find(pim>0));
% % yfp_zcorr(:,1) = ix*px/ef;
% % yfp_zcorr(:,2) = iy*px/ef;
% % yfp_zcorr(:,3) = iz*px/ef;
% %
% % tic
% % Y_Mdl = KDTreeSearcher(yfp_zcorr); % all nonDAN synapse positions
% % toc
%
% %% get local maxima
% % % JIB = readtiff('BassonJuxH.tif');
% % % JIH = readtiff('HomerJuxB.tif');
% % %
% % % tic
% % % B_LM_filteredData = GU_getLocalMaxima(JIB, 'ExpansionFactor', 4.04, 'Sigmas', [2.5, 2.5],'PixelSize', 0.097,'minDia', 0.0); toc
% % %
% % % tic
% % % H_LM_filteredData = GU_getLocalMaxima(JIH, 'ExpansionFactor', 4.04, 'Sigmas', [2.5, 2.5],'PixelSize', 0.097,'minDia', 0.0); toc
%
% % load('20180524_localMaxima_nonInterpolatedHomerBassoonJuxaposed.mat')
% %
% % B_LM_filteredData(:,1) = B_LM_filteredData(:,1)*px/ef;
% % B_LM_filteredData(:,2) = B_LM_filteredData(:,2)*px/ef;
% % B_LM_filteredData(:,3) = B_LM_filteredData(:,3)*zAniso*px/ef;
% %
% % H_LM_filteredData(:,1) = H_LM_filteredData(:,1)*px/ef;
% % H_LM_filteredData(:,2) = H_LM_filteredData(:,2)*px/ef;
% % H_LM_filteredData(:,3) = H_LM_filteredData(:,3)*zAniso*px/ef;
% %
% % tic
% % B_Mdl = KDTreeSearcher(B_LM_filteredData);
% % toc
% %
% % tic
% % H_Mdl = KDTreeSearcher(H_LM_filteredData);
% % toc
% %% recalculate densities of spots
% B_filteredData_zcorr(:,1) = B_filteredData(:,1)*px/ef;
% B_filteredData_zcorr(:,2) = B_filteredData(:,2)*px/ef;
% B_filteredData_zcorr(:,3) = B_filteredData(:,3)*zAniso*px/ef;
%
% tic
% B_Mdl = KDTreeSearcher(B_filteredData_zcorr); % all nonDAN synapse positions
% toc
%
% H_filteredData_zcorr(:,1) = H_filteredData(:,1)*px/ef;
% H_filteredData_zcorr(:,2) = H_filteredData(:,2)*px/ef;
% H_filteredData_zcorr(:,3) = H_filteredData(:,3)*zAniso*px/ef;
%
% tic
% H_Mdl = KDTreeSearcher(H_filteredData_zcorr); % all nonDAN synapse positions
% toc
%
% %% find yfp associated pairs
% %
% % tic
% % [idx,D] = knnsearch(Y_Mdl, B_filteredData_zcorr, 'k',1);
% % toc
% %
% % idx300nmDistCutoff = D<0.3;
% % perimYFPidx = idx(idx300nmDistCutoff);
% %
% % yfp_zcorr_BAss = yfp_zcorr(perimYFPidx,:);
% %
% % tic
% % [~,D_y2b] = knnsearch(B_Mdl, yfp_zcorr_BAss, 'k',1);
% % toc
% %
% % tic
% % [~,D_y2h] = knnsearch(H_Mdl, yfp_zcorr_BAss, 'k',1);
% % toc
% %
% % idxBasson2YFP_300nmDistCutoff = D_y2b<0.3;
% % AssBasson2YFP = B_filteredData_zcorr(idxBasson2YFP_300nmDistCutoff,:);
% % idxHomer2YFP_300nmDistCutoff = D_y2h<0.3;
% % AssHomer2YFP = H_filteredData_zcorr(idxHomer2YFP_300nmDistCutoff,:);
% %
% %
% % tic
% % HY_Mdl = KDTreeSearcher(AssHomer2YFP); % all nonDAN synapse positions
% % toc
% %
% % tic
% % [~,D_H2B] = knnsearch(HY_Mdl, AssBasson2YFP, 'k',1);
% % toc
% %
% % tic
% % [~,D] = knnsearch(B_Mdl, H_filteredData_zcorr, 'k',1);
% % toc
%
% %% distance - homer to closest bassoon
% tic
% [~,D] = knnsearch(H_Mdl, B_LM_filteredData, 'k',1);
% toc
%
% ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
%     'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
%  histogram(D*1000,0:10:600)
% set(ha, 'FontSize', 6);
% xlim([0 600])
% xticks(0:150:600)
% ylabel('# events');
% xlabel('Distance (nm)');
% legend(['n=' num2str(numel(D))])
% f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date 'DistanceFromHomer2ClosestBasoon.eps']);
% %%
%
% % tic
% % Mdl = KDTreeSearcher(B_filteredData_zcorr); % all nonDAN synapse positions
% % toc
% %
% % tic
% % [idx,D] = knnsearch(Mdl, H_filteredData_zcorr(H_YFPidx,:), 'k',1);
% % toc
% tic
% [idx,D] = knnsearch(B_Mdl, Homer_centroid, 'k',1);
% toc
%
% ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
%     'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
%  histogram(D*1000,0:25:600)
% set(ha, 'FontSize', 6);
% % xlim([0 600])
% xticks(0:150:1000)
% ylabel('# events');
% xlabel('Distance (nm)');
% legend(['n=' num2str(numel(D))])
% f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date 'DistanceFromYFPHomer2ClosestBasoon.eps']);
%
% %% yfp associated homer
%
% tic
% [idx,D] = knnsearch(B_Mdl, Homer_centroid(find(yfp_h_idx),:), 'k',1);
% toc
%
% ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
%     'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
%  histogram(D*1000,0:25:600)
% set(ha, 'FontSize', 6);
% xlim([0 600])
% xticks(0:150:1000)
% ylabel('# events');
% xlabel('Distance (nm)');
% legend(['n=' num2str(numel(D))])
% f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date 'CentroidDistanceFromYFPHomer2ClosestBasoon.eps']);
% [mu, sigma, ~, ~] = fitGaussianModeToCDF(D, 'Display',true)
% %% load homer bassoon 3d region props
% % imB = readtiff('ch1_1008_clean.tif');
% % iimB = GU_interp3DlargeVolume(imB, 'zAniso', 0.25/0.097, 'nzvols', 7);
% % writetiff(iimB, 'interBasson.tif');
% % clear iimB imB
% %
% % imH = readtiff('ch2_1008_clean.tif');
% % iimH = GU_interp3DlargeVolume(imH, 'zAniso', 0.25/0.097, 'nzvols', 7);
% % writetiff(iimH, 'interHomer.tif');
% % %%
% % iimB = readtiff('interBasson.tif');
% % iimH = readtiff('interHomer.tif');
% % %%
% %
% % B_rp3  = regionprops3(logical(iimB));
% % H_rp3  = regionprops3(logical(iimH));
% % save('20180522_interpolated_HomerBassson_rp3.mat', 'B_rp3', 'H_rp3');
% %%
% load('20180522_interpolated_HomerBassson_rp3.mat');
%
% %%
%
% B_MA = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), B_rp3, 'unif',0);
% B_MA = horzcat(B_MA{:})*px/ef;
% H_MA = arrayfun(@(x) max([x.FirstAxisLength, x.SecondAxisLength, x.ThirdAxisLength]), H_rp3, 'unif',0);
% H_MA = horzcat(H_MA{:})*px/ef;
%
%
% %
% % figure
% % i = 1
% % x = B_rp3(i).FirstAxis;
% % y = B_rp3(i).SecondAxis;
% % z = B_rp3(i).ThirdAxis;
% % quiver3(zeros(3,1),zeros(3,1),zeros(3,1),x',y',z');
% % hold on
% % quiver3(zeros(3,1),zeros(3,1),zeros(3,1),x',x',x');
% %% synaptic width
% ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', true,...
%     'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1.5 1.5]);
%
% % histogram(H_MA*1000,0:100:4000)
% % hold on
% % histogram(B_MA*1000,0:100:4000)
% GU_stairsHistPlot({[],[],[],B_MA,H_MA}, 'ShowECDF', true, 'BinSize',.010,'CreateFigure', false,...
%         'PercentileEnd', 99.9)
% set(ha, 'FontSize', 6);
% xlim([0 2])
% xticks(0:0.5:5)
% yyaxis left
% ylabel('Relative Freq.');
% xlabel('Major Axis (\mum)');
% legend([''],['Bassoon n=' num2str(numel(B_MA))],[''],['Homer n=' num2str(numel(H_MA))])
% f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date 'MajorAxis_paired.eps']);
%
% %%
% ms = 10
% ha = setupFigure(1,1, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
%     'XSpace', [1.5 0.85 0.85], 'YSpace', [1.5 0.85 0.5]);
% axes(ha(1))
% a = jet(256);
% params = {'markerSize', ms; 'colMap', a(end:-1:1,:)};
% X = [ H_MA;D'*1000];
% densityScatter(X, params);
%
% % figure, scatter(D*1000, H_MA)
%
% %% spine detection
% ef = 4.04;
% px = .097; %in um
% pz = 0.25;
% zAniso = pz/px;
%
% % rt = 'E:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch0_Analysis_1ch\1008\';
% rt = 'D:\GoogleDrive_DataTransfer\ExLLSMData_SharedwithGokulRui\160126_Sample4_Y10_YFP_HB\ch0_Analysis_1ch\1008\';
% cd(rt)
% fn = 'ch0_1008_clean_zpart_2.tif';
%
% GU_ExM_extractSpineParameters([rt fn], 'zAniso', zAniso, ...
%     'PixelSize', 0.097, 'MinBranchLength', 10, 'RemoveOverlap', 0, 'CropScalingFactor', 5,'MinVoxelVolume', 0)
%
% %%
% % p = 8;
% % nz = 100;
% % for i = 1:p
% %    writetiff(im(:,:,((i-1)*100+1):(i*100)), [fn(1:end-4) '_zpart_' num2str(i) '.tif']);
% %    i
% % end
%
