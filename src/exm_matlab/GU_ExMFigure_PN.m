
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/PN_stereo/';
cd(rt)
load('PN_area_vol.mat')

% curation
PN(3).cell(3).vol(1) = PN(3).cell(3).vol(1)+PN(3).cell(3).vol(5);
PN(3).cell(3).area(1) = PN(3).cell(3).area(1)+PN(3).cell(3).area(5);


% figure


ha = setupFigure(2,1, 'AxesWidth', 20, 'AxesHeight', 10,'SameAxes', false,...
    'XSpace', [2 1.85 1.85], 'YSpace', [1.5 1.85 1.5]);
axes(ha(1))
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
c = {ceG, ceO, ceR, ceP, ceB};
m = {'s', 'd', 'o'};
for k = 1:numel(PN)
    x = k;
    for j = 1:numel(PN(k).cell)
        
        v = PN(k).cell(j).vol;
        v = v(v>4);
        a = PN(k).cell(j).area;
        a = a(v>4);
        axes(ha(1))
        plot(x, v, m{j}, 'MarkerSize',10, 'Color', c{k}, 'MarkerFaceColor', c{k});
        hold on
        
        axes(ha(2))
        plot(x, a, m{j}, 'MarkerSize',10, 'Color', c{k}, 'MarkerFaceColor', c{k});
        
        
        x = x+0.2;
    end
    
end
axes(ha(1))
xlim([0.5 5.5])
ylim([0 35])
ylabel('Volume (\mum^3)');
xlabel('Drosophila');

axes(ha(2))
xlim([0.5 5.5])
ylim([0.01 1000])
ylabel('Surface Area (\mum^2)');
xlabel('Drosophila');
set(gca, 'YScale', 'log')
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'PN_volume_area.eps']);
%%

% figure for the PN - syd1-nc82
ef = 3.99;
px = 0.097;
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/PN_nc82Syd1_surfaceAreaVol/';
cd(rt)
CAidx = {20:22, 33:35, 24:27};
LHidx = {1:19, 1:32, 1:22};

fn = dir([rt '*.mat']);
fn = {fn.name}';

ha = setupFigure(2,1, 'AxesWidth', 10, 'AxesHeight', 5,'SameAxes', false,...
    'XSpace', [2 1.85 1.85], 'YSpace', [1.5 1.85 1.5]);
axes(ha(1))
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
c = {ceB, ceO, ceG, ceP, ceB};

m = {'s', 'd', 'o'};
x = 1;
xx = 2;
for k = 1:3
    load([rt fn{k}]);
    
    v_CA = Volume(CAidx{k})/px^3*(px/ef)^3;
    v_LH = Volume(LHidx{k})/px^3*(px/ef)^3;
    
    a_CA = SurfaceArea(CAidx{k})/px^2*(px/ef)^2;
    a_LH = SurfaceArea(LHidx{k})/px^2*(px/ef)^2;
    
    axes(ha(1))
    plot(x+0.1.*rand(numel(v_CA),1), v_CA, m{k}, 'MarkerSize',5, 'Color', c{k}, 'MarkerFaceColor', c{k});
    plot(xx+0.1.*rand(numel(v_LH),1), v_LH, m{k}, 'MarkerSize',3, 'Color', c{k}, 'MarkerFaceColor', c{k});
    
    axes(ha(2))
    plot(x+0.1.*rand(numel(a_CA),1), a_CA, m{k}, 'MarkerSize',5, 'Color', c{k}, 'MarkerFaceColor', c{k});
    plot(xx+0.1.*rand(numel(a_LH),1), a_LH, m{k}, 'MarkerSize',3, 'Color', c{k}, 'MarkerFaceColor', c{k});
    
    x = x+0.2;
    xx = xx+0.2;
end

axes(ha(1))
xlim([0.75 2.75])
ylim([0.1 100])
ylabel('Volume (\mum^3)');
ax = gca;
ax.GridColor = 'white';
ax.Color = 'black';
ax.YColor = 'white';
ax.XColor = 'white';
ax.ZColor = 'white';
% xlabel('Drosophila');
set(gca, 'YScale', 'log')
axes(ha(2))
xlim([0.75 2.75])
ylim([1 100])
ylabel('Surface Area (\mum^2)');
% xlabel('Drosophila');

set(gcf,'color','black')
set(gcf,'DefaultAxesColor','k')
ax = gca;
ax.GridColor = 'white';
ax.Color = 'black';
ax.YColor = 'white';
ax.XColor = 'white';
ax.ZColor = 'white';


set(gca, 'YScale', 'log')
f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[date 'PN_nc82Syd1_BoutonsVolume_Area.eps']);

