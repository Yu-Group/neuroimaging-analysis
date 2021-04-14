%% colocalization of the v5 and nc82

% rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/v5_nc82/';
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/v5_nc82/';
cd(rt)
load('v5_nc82_coord.mat');
nc82 = nc82(logical(maskidx_nc82),:);
v5 = v5(logical(maskidx_v5),:);

nc82(:,1) = nc82(:,1) * epx;
nc82(:,2) = nc82(:,2) * epx;
nc82(:,3) = nc82(:,3) * epz;

v5(:,1) = v5(:,1) * epx;
v5(:,2) = v5(:,2) * epx;
v5(:,3) = v5(:,3) * epz;
%% generate nc82 model

tic
Mdl = KDTreeSearcher(v5);
toc
% search for the closest V5
tic
[~,D] = knnsearch(Mdl, nc82, 'k',2);
toc
D2D = D(:,1);
%% figure plot distance between nc82 and closest v5

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;

ha = setupFigure(1,1, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', true,...
    'XSpace', [1.1 1 1], 'YSpace', [1 1 1]);
set(ha, 'FontSize', 6);

axes(ha(1))
histogram(D2D*1000, 150, 'EdgeColor', cfK, 'FaceColor', cfK,'BinWidth',25)
    title('nc82 distance to closest V5 ')
    legend(['n = ' num2str(numel(D2D))])
    legend('boxoff')
    set(gca,'fontsize',6, 'FontName', 'Helvetica')
    xlim([-10 500]);
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
print(f0, '-painters','-depsc', '-loose',[date ' nc82_v5_DistancePlots.eps']);

%%

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;

ha = setupFigure(1,1, 'AxesWidth', 5, 'AxesHeight', 5,'SameAxes', true,...
    'XSpace', [1.1 1 1], 'YSpace', [1 1 1]);
set(ha, 'FontSize', 6);

axes(ha(1))
[f,x] = ecdf(D2D*1000);
plot(x,f)
    title('nc82 distance to closest V5 ')
    legend(['n = ' num2str(numel(D2D))])
    legend('boxoff')
    set(gca,'fontsize',6, 'FontName', 'Helvetica')
    xlim([-10 500]);
    ax = gca;
    ax.BoxStyle = 'full';
    box on
%     yyaxis left
%     ylim([0 1500]);
    xticks(0:50:1500)
%     grid on
%     if i == 1
        ylabel('Cumulative Frequency');
%     end
%     if ~isEven(i)
        xlabel('Distance to V5 (nm)');
%     end

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date ' nc82_v5_DistancePlots_ecdf.eps']);