rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/PN_nc82Syd1_surfaceAreaVol/';
cd(rt);
load('allnc82_PNIdx.mat')
load('allSyd1_PNIdx.mat')
%%
ef = 3.99;
epx = 0.097/ef;
epz = 0.18/ef;
%%
nc82(:,1) = allnc82filteredData(:,1) * epx;
nc82(:,2) = allnc82filteredData(:,2) * epx;
nc82(:,3) = allnc82filteredData(:,3) * epz;

syd1(:,1) = allSyd1filteredData(:,1) * epx;
syd1(:,2) = allSyd1filteredData(:,2) * epx;
syd1(:,3) = allSyd1filteredData(:,3) * epz;
%% generate nc82 model

tic
Mdl = KDTreeSearcher(syd1);
toc
% search for the closest V5
tic
[~,D] = knnsearch(Mdl, nc82, 'k',1);
toc
D2D = D(:,1);

%%
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;

ha = setupFigure(1,1, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', true,...
    'XSpace', [1.1 1 1], 'YSpace', [1 1 1]);
set(ha, 'FontSize', 6);

axes(ha(1))
histogram(D2D*1000, 150, 'EdgeColor', cfK, 'FaceColor', cfK,'BinWidth',25)
    title('syd1 distance to closest nc82 ')
    legend(['n = ' num2str(numel(D2D))])
    legend('boxoff')
    set(gca,'fontsize',6, 'FontName', 'Helvetica')
    xlim([-10 1000]);
    ax = gca;
    ax.BoxStyle = 'full';
    box on
%     yyaxis left
%     ylim([0 1500]);
    xticks(0:50:1500)
%     grid on
%     if i == 1
        ylabel('# nc82 puncta');
%     end
%     if ~isEven(i)
        xlabel('Distance to V5 (nm)');
%     end

f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[date ' nc82_syd1_DistancePlots.eps']);