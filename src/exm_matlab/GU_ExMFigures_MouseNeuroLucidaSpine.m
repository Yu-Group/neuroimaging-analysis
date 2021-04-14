6% GU_ExMFigures_MouseNeuroLucidaSpine

rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/NeurolucidaData/';
load([rt 'data_coordinates_offsets.mat']);
cd(rt)
%%



SpineType = arrayfun(@(x) x.SpineType, data, 'unif',0);
SpineType_all  = vertcat(SpineType{:});

Homer_HeadRadius = arrayfun(@(x) x.nHomerHeadradius, data, 'unif',0);
Homer_HeadRadius_all = horzcat(Homer_HeadRadius{:});

Homer_HalfMicron = arrayfun(@(x) x.nHomerHalfMicron, data, 'unif',0);
Homer_HalfMicron_all = horzcat(Homer_HalfMicron{:});

Homer_dist = arrayfun(@(x) x.SpineHomerDist, data, 'unif',0);
Homer_dist_all = vertcat(Homer_dist{:});

BackboneLength = arrayfun(@(x) x.BackboneLength, data, 'unif',0);
BackboneLength_all = vertcat(BackboneLength{:});

NeckDiameter = arrayfun(@(x) x.ExMNeckDiameterm, data, 'unif',0);
NeckDiameter_all = vertcat(NeckDiameter{:});

nHomer_HeadRadius = arrayfun(@(x) sum(x.nHomerHeadradius), data, 'unif',1);
perHomer_HeadRadius = arrayfun(@(x) sum(x.nHomerHeadradius)/numel(x.nHomerHeadradius), data, 'unif',1);

uSpineType = unique(SpineType_all);
uSpineType = uSpineType(end:-1:1);
%
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
pos = {1:28, 1:4, 5:8, 9:12, 12:16, 17:20, 21:24, 25:28};



% legend('Search: Spine Head Radius')
allCat = {'All Positions', Filename{1:4:end}};

%% Spine Head Radius
for k = 1:2
    ha = setupFigure(3,8, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
        'XSpace', [1.5 1.5 1.5], 'YSpace', [2 2 2]);
    set(ha, 'FontSize', 6);
    tmp = jet;
    cmap = tmp(end:-1:1,:);
    
    for i = 1:numel(pos)
        if k == 1
            nHomer = horzcat(Homer_HeadRadius{pos{i}}); % head radius
        elseif k ==2
            nHomer = horzcat(Homer_HalfMicron{pos{i}}); % half micron search
        end
        idx = logical(nHomer);
        SpineType_select = vertcat(SpineType{pos{i}});
        for j = 1:numel(uSpineType)
            ST{j} = logical(strcmp(SpineType_select, uSpineType{j}));
            %         nS(j) = sum(double(strcmp(SpineType_select(idx), uSpineType{j})).*nHomer(idx)')./sum(idx)*100;
            nS(j) = sum(logical(strcmp(SpineType_select(idx), uSpineType{j})))./sum(logical(strcmp(SpineType_select, uSpineType{j})));
        end
        axes(ha(i))
        set(ha(i), 'FontSize', 8);
        bar(nS,'FaceColor',cfO,'EdgeColor',ceO);
        set(ha(i), 'XTickLabel',uSpineType, 'XTick',1:numel(uSpineType), 'FontSize',8)
        rotateticklabel(ha(i), 45)
        ylabel('Fraction Spines with Homer')
        title(allCat{i}(1:4))
        ylim([0 1])
        
        axes(ha(i+8))
        set(ha(i+8), 'FontSize', 8);
        BackboneLength_select = vertcat(BackboneLength{pos{i}});
        %     GU_stairsHistPlot({BackboneLength_select(idx)*1000}, 'BinSize',50, 'xlabel', 'Backbone Length', 'ylabel','Rel Freq.', 'CreateFigure',false)
        %     edgesBackbone = 0:0.25:5;
        %     histogram(BackboneLength_select(idx),edgesBackbone,'FaceColor',cfG,'EdgeColor',ceG,'Normalization','probability')
        %     ylabel('Freq')
        %     xlabel('Backbone Length (µm)')
        %     xlim([0 5])
        clear X
        for j = 1:numel(uSpineType)
            ST{j} = logical(strcmp(SpineType_select, uSpineType{j}));
            X(:,1) = ones(1, numel(BackboneLength_select))*(j-1)*2;
            X(:,2) = BackboneLength_select;
            f = idx & ST{j}';
            if sum(f)>1
                densityScatter(X(f,:)', {'colMap', cmap}); %, {'colMap', jet}
            elseif sum(f) == 1
                scatter(X(f,1), X(f,2),20,'b', 'filled');
            end
            X(:,1) = ones(1, numel(BackboneLength_select))*(j-1)*2+1;
            X(:,2) = BackboneLength_select;
            f = ~idx & ST{j}';
            if sum(f)>1
                densityScatter(X(f,:)', {'colMap', cmap});
            end
        end
        ylim([0 5])
        xlim([-1 9])
        set(ha(i), 'XTickLabel',uSpineType, 'XTick',1:2:2*numel(uSpineType), 'FontSize',8)
        rotateticklabel(ha(i), 45)
        ylabel('Backbone Length (µm)')
        
        axes(ha(i+16))
        set(ha(i+16), 'FontSize', 8);
        edgesNeck = 0:0.05:1;
        NeckDia_select = vertcat(NeckDiameter{pos{i}});
        %     GU_stairsHistPlot({NeckDia_select(idx)*1000}, 'BinSize',50, 'xlabel', 'Backbone Length', 'ylabel','Rel Freq.', 'CreateFigure',false)
        %     histogram(NeckDia_select(idx),edgesNeck,'FaceColor',cfB,'EdgeColor',ceB,'Normalization','probability')
        %     ylabel('Freq')
        %     xlabel('Neck Diameter (µm)')
        %     xlim([0 1])
        clear X
        for j = 1:numel(uSpineType)
            ST{j} = logical(strcmp(SpineType_select, uSpineType{j}));
            X(:,1) = ones(1, numel(NeckDia_select))*(j-1)*2;
            X(:,2) = NeckDia_select;
            f = idx & ST{j}';
            if sum(f)>1
                densityScatter(X(f,:)', {'colMap', cmap});
            end
            X(:,1) = ones(1, numel(NeckDia_select))*(j-1)*2+1;
            X(:,2) = NeckDia_select;
            f = ~idx & ST{j}';
            if sum(f)>1
                densityScatter(X(f,:)', {'colMap', cmap});
            end
        end
        ylim([0 1])
        xlim([-1 9])
        set(ha(i), 'XTickLabel',uSpineType, 'XTick',1:2:2*numel(uSpineType), 'FontSize',8)
        rotateticklabel(ha(i), 45)
        ylabel('Neck Diameter (µm)')
    end
    f0 = gcf();
    if k ==1
        print(f0, '-painters','-depsc', '-loose',['SpineHeadRadius_density.eps']);
        export_fig(['SpineHeadRadius_density.pdf']);
    elseif k ==2
        print(f0, '-painters','-depsc', '-loose',['SpineHalfMicron_density.eps']);
        export_fig(['SpineHalfMicron_density.pdf']);
    end
end
%% Spine 0.5 um search radius

ha = setupFigure(3,8, 'AxesWidth', 3, 'AxesHeight', 3,'SameAxes', false,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [2 2 2]);
set(ha, 'FontSize', 6);

for i = 1:numel(pos)
    nHomer = horzcat(Homer_HalfMicron{pos{i}});
    idx = logical(nHomer);
    SpineType_select = vertcat(SpineType{pos{i}});
    for j = 1:numel(uSpineType)
        %         nS(j) = sum(double(strcmp(SpineType_select(idx), uSpineType{j})).*nHomer(idx)')./sum(idx)*100;
        nS(j) = sum(logical(strcmp(SpineType_select(idx), uSpineType{j})))./sum(logical(strcmp(SpineType_select, uSpineType{j})));
    end
    axes(ha(i))
    set(ha(i), 'FontSize', 8);
    bar(nS,'FaceColor',cfO,'EdgeColor',ceO);
    set(ha(i), 'XTickLabel',uSpineType, 'XTick',1:numel(uSpineType), 'FontSize',8)
    rotateticklabel(ha(i), 45)
    ylabel('Fraction Spines with Homer')
    title(allCat{i}(1:4))
    ylim([0 1])
    
    axes(ha(i+8))
    set(ha(i+8), 'FontSize', 8);
    BackboneLength_select = vertcat(BackboneLength{pos{i}});
    %     GU_stairsHistPlot({BackboneLength_select(idx)*1000}, 'BinSize',50, 'xlabel', 'Backbone Length', 'ylabel','Rel Freq.', 'CreateFigure',false)
    edgesBackbone = 0:0.25:5;
    histogram(BackboneLength_select(idx),edgesBackbone,'FaceColor',cfG,'EdgeColor',ceG,'Normalization','probability')
    ylabel('Freq')
    xlabel('Backbone Length (µm)')
    xlim([0 5])
    
    axes(ha(i+16))
    set(ha(i+16), 'FontSize', 8);
    edgesNeck = 0:0.05:1;
    NeckDia_select = vertcat(NeckDiameter{pos{i}});
    %     GU_stairsHistPlot({NeckDia_select(idx)*1000}, 'BinSize',50, 'xlabel', 'Backbone Length', 'ylabel','Rel Freq.', 'CreateFigure',false)
    histogram(NeckDia_select(idx),edgesNeck,'FaceColor',cfB,'EdgeColor',ceB,'Normalization','probability')
    ylabel('Freq')
    xlabel('Neck Diameter (µm)')
    xlim([0 1])
end
f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',['SpineHalfMicron.eps']);
% export_fig(['SpineHalfMicron.pdf']);

%% re-plotting spine statistics from Neurolucida
%1) "Backbone Length(µm)" (column E) vs. "Head Diameter(µm)" (column R)  and (2) "Neck Backbone Length(µm)" (column T) vs. "ExM Neck Diameter (µm)" (column W)
rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/Spinestatistics/';
% rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/Spinestatistics/';
cd(rt)
load('20180126_SpineStatistics.mat')


allData = [PosA;PosB;PosC;PosD;PosE;PosF;PosG];
regions = {'allData', 'PosA', 'PosB', 'PosC', 'PosD', 'PosE', 'PosF', 'PosG'};
% regions = {'PosA', 'PosB', 'PosC', 'PosD', 'PosE', 'PosF', 'PosG'};

uSpineType = unique(allData.SpineType);
uSpineType = uSpineType(end:-1:1);
%%

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
edges1 = 0:0.02:1;
edges5 = 0:0.1 :5;
for i = 1:numel(regions)
    
    d = eval(regions{i});
    backboneL = [-1;-1;-1;-1; d.BackboneLengthm];
    HeadDiameter = [-1;-1;-1;-1; d.HeadDiameterm];
    NeckBackboneL = [-1;-1;-1;-1; d.NeckBackboneLengthm];
    ExMNeckDiameter = [-1;-1;-1;-1; d.ExMNeckDiameterm];
    SpineT = [uSpineType; d.SpineType];
   
%     scatter3(backboneL,HeadDiameter,ExMNeckDiameter,100,NeckBackboneL,'.')
%     xlim([0 5])
%     ylim([0 1])
%     zlim([0 1])
%     colormap jet
%     % %     Mushroom: blue
%     % % Thin: red
%     % % Stubby: green
%     % % Filopodia: yellow
%     
%     
%     %
%     ha = scatterhist(backboneL,HeadDiameter,'Group',SpineT,'Kernel','off','Location','NorthEast',...
%         'Direction','out','Color','bgry','LineStyle',{'-','-.',':','--'},...
%         'LineWidth',[1,1,1,1],'Marker','....','MarkerSize',[10,10,10,10]);
%     xlim(ha(1),[0 5]);
%     ylim(ha(1),[0 1]);
%     set(gcf, 'Position', [500, 500, 500, 500])
%     xlabel('Backbone Length (µm)')
%     ylabel('Head Diameter (µm)')
%     pause
%     histogram(ha(3), HeadDiameter(HeadDiameter<=1),edges1,'orientation','horizontal','EdgeColor', ceK, 'FaceAlpha', 0)
%     ylim(ha(3),[0 1]);
%     xlim(ha(3),[0 30]);
%     yticklabels(ha(3),[]);
%     pause
%     histogram(ha(2), backboneL,edges5,'orientation','vertical','EdgeColor', ceK, 'FaceAlpha', 0)
%     xlim(ha(2),[0 5]);
%     ylim(ha(2),[0 40]);
%     xticklabels(ha(2),[]);
%     ax = gca;
%     ax.GridColor = 'white';
%     ax.Color = 'black';
%     ax.YColor = 'white';
%     ax.XColor = 'white';
%     ax.ZColor = 'white';
%     f0 = gcf();
%     print(f0, '-painters','-depsc', '-loose',[date regions{i} '_backBoneLength_HeadDiameter.eps']);
%     close
%     %
%     
%     ha = scatterhist(NeckBackboneL,ExMNeckDiameter,'Group',SpineT,'Kernel','off','Location','NorthEast',...
%         'Direction','out','Color','bgry','LineStyle',{'-','-.',':','--'},...
%         'LineWidth',[1,1,1,1],'Marker','....','MarkerSize',[10,10,10,10]);
%     xlim(ha(1),[0 5]);
%     ylim(ha(1),[0 1]);
%     set(gcf, 'Position', [500, 500, 500, 500])
%     xlabel('Neck Backbone Length (µm)')
%     ylabel('ExM Neck Diameter (µm)')
%     
%     pause
%     histogram(ha(3), ExMNeckDiameter,edges1,'orientation','horizontal','EdgeColor', ceK, 'FaceAlpha', 0)
%     ylim(ha(3),[0 1]);
%     xlim(ha(3),[0 60]);
%     yticklabels(ha(3),[]);
%     pause
%     histogram(ha(2), NeckBackboneL,edges5,'orientation','vertical','EdgeColor', ceK, 'FaceAlpha', 0)
%     xlim(ha(2),[0 5]);
%     ylim(ha(2),[0 80]);
%     xticklabels(ha(2),[]);
%     ax = gca;
%     ax.GridColor = 'white';
%     ax.Color = 'black';
%     ax.YColor = 'white';
%     ax.XColor = 'white';
%     ax.ZColor = 'white';
%     f0 = gcf();
%     print(f0, '-painters','-depsc', '-loose',[date regions{i} '_NeckBackboneL_ExMNeckDiameter.eps']);
%     close
%     
    %     
    ha = scatterhist(NeckBackboneL,HeadDiameter,'Group',SpineT,'Kernel','off','Location','NorthEast',...
        'Direction','out','Color','bgry','LineStyle',{'-','-.',':','--'},...
        'LineWidth',[1,1,1,1],'Marker','....','MarkerSize',[10,10,10,10]);
    [rho_NKL_HD(i,1),pval_NKL_HD(i,1)] = corr(NeckBackboneL, HeadDiameter);
    
    [rho_NKL_ND(i,1),pval_NKL_ND(i,1)] = corr(NeckBackboneL, ExMNeckDiameter);
    [rho_BL_ND(i,1),pval_BL_ND(i,1)] = corr(backboneL, ExMNeckDiameter);
     [rho_HD_ND(i,1),pval_HD_ND(i,1)] = corr(HeadDiameter, ExMNeckDiameter);
     
     [rho_HD_BL(i,1),pval_HD_BL(i,1)] = corr(HeadDiameter, backboneL);
     
     backboneL, HeadDiameter, NeckBackboneL, ExMNeckDiameter
    xlim(ha(1),[0 5]);
    ylim(ha(1),[0 1]);
    set(gcf, 'Position', [500, 500, 500, 500])
    xlabel('Neck Backbone Length (µm)')
    ylabel('Head Diameter (µm)')
    
%     pause
    histogram(ha(3), HeadDiameter(HeadDiameter<=1),edges1,'orientation','horizontal','EdgeColor', ceK, 'FaceAlpha', 0)
    ylim(ha(3),[0 1]);
    xlim(ha(3),[0 30]);
    yticklabels(ha(3),[]);
    
%     pause
    histogram(ha(2), NeckBackboneL,edges5,'orientation','vertical','EdgeColor', ceK, 'FaceAlpha', 0)
    xlim(ha(2),[0 5]);
    ylim(ha(2),[0 80]);
    xticklabels(ha(2),[]);
% 
%     f0 = gcf();
%     print(f0, '-painters','-depsc', '-loose',[date regions{i} '_NeckBackboneL_HeadDiameter.eps']);
%     close
end
%% stats


for i = 1:numel(regions)
    clc
    i
    d = eval(regions{i});
    backboneL = d.BackboneLengthm;
    HeadDiameter = d.HeadDiameterm;
    NeckBackboneL = d.NeckBackboneLengthm;
    ExMNeckDiameter = d.ExMNeckDiameterm;
    SpineT = d.SpineType;
    
    [median(backboneL), mad(backboneL,1)]
    [median(HeadDiameter), mad(HeadDiameter,1)]
    [median(NeckBackboneL), mad(NeckBackboneL,1)]
    [median(ExMNeckDiameter), mad(ExMNeckDiameter,1)]
    SpineCat = categories(SpineT)
    SpineCount = countcats(SpineT)
    pause
end