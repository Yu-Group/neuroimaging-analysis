%GU_ExMFigures_MBON_fly
%%
rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/DataForFigures/Ex23_MBON_dedritesvsaxons_z0p18/';
cd(rt)

load('SurfaceArea_Vol.mat')
AxonDAN = load('Ex1_a1-axon-boutons-ch0-membrane_246_densityData.mat', 'Vol2PunctaIdx','filteredData');
DendriteDAN = load('Ex2_a1-dendrites-ch0-membrane_246_densityData.mat', 'Vol2PunctaIdx','filteredData');

px = 0.097;
zAniso = 0.18/px;
ef = 3.99;

nc82_AxonDANCoord = AxonDAN.filteredData(AxonDAN.Vol2PunctaIdx,:);
nc82_DendriteDANCoord = DendriteDAN.filteredData(DendriteDAN.Vol2PunctaIdx,:);

axonSize = [3690,2190,632];
DendriteSize = [952,1722,613];

[axonMask_coord(:,2),axonMask_coord(:,1),axonMask_coord(:,3)] = ind2sub(axonSize, axonMaskIdx);
[DendriteMask_coord(:,2),DendriteMask_coord(:,1),DendriteMask_coord(:,3)] = ind2sub(DendriteSize, DendriteMaskIdx);

maskedAxonDAN = ismember(nc82_AxonDANCoord,axonMask_coord,'rows');
maskedDendriteDAN = ismember(nc82_DendriteDANCoord,DendriteMask_coord,'rows');

nAxonDAN = sum(maskedAxonDAN);
nDendriteDAN = sum(maskedDendriteDAN);

DendriteMask_coord(:,3) = DendriteMask_coord(:,3) * zAniso;
axonMask_coord(:,3) = axonMask_coord(:,3) * zAniso;
nc82_AxonDANCoord(:,3) = nc82_AxonDANCoord(:,3) * zAniso;
nc82_DendriteDANCoord(:,3) = nc82_DendriteDANCoord(:,3) * zAniso;

DendriteMask_coord = DendriteMask_coord * px/ef;
axonMask_coord = axonMask_coord * px/ef;
nc82_AxonDANCoord = nc82_AxonDANCoord * px/ef;
nc82_DendriteDANCoord = nc82_DendriteDANCoord * px/ef;
%%
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;

ha = setupFigure(1,2, 'AxesWidth',10, 'AxesHeight', 10,'SameAxes', false,...
    'XSpace', [1.25 2 1.25], 'YSpace', [1.5 0.85 1]);
axes(ha(1))
scatter3(DendriteMask_coord(1:50:end,1),DendriteMask_coord(1:50:end,2),DendriteMask_coord(1:50:end,3), 1, ceB,'.', 'MarkerFaceAlpha', 0.1)
view(3)
axis equal
xlim([0 40])
ylim([0 30])
zlim([0 30])
xticks([0:10:50])
yticks([0:10:50])
zticks([0:10:50])
xlabel('\mum')
grid on
% hold on
axes(ha(2))
scatter3(nc82_DendriteDANCoord(maskedDendriteDAN,1),nc82_DendriteDANCoord(maskedDendriteDAN,2),nc82_DendriteDANCoord(maskedDendriteDAN,3), 10, cfO, 'o')
view(3)
axis equal
xlim([0 40])
ylim([0 30])
zlim([0 30])
xticks([0:10:50])
yticks([0:10:50])
zticks([0:10:50])
grid on
xlabel('\mum')
f0 = gcf();
print(f0, '-painters','-dpdf', '-loose',[date ' Dendrite_Mask_nc82.pdf']);
%%

ha = setupFigure(1,2, 'AxesWidth',10, 'AxesHeight', 10,'SameAxes', false,...
    'XSpace', [1.25 2 1.25], 'YSpace', [1.5 0.85 1]);
axes(ha(1))
scatter3(axonMask_coord(1:50:end,1),axonMask_coord(1:50:end,2),axonMask_coord(1:50:end,3), 1, ceB,'.', 'MarkerFaceAlpha', 0.1)
view(3)
axis equal
xlim([0 50])
ylim([0 90])
zlim([0 30])
xticks([0:10:50])
yticks([0:10:100])
zticks([0:10:50])
xlabel('\mum')
grid on

axes(ha(2))
scatter3(nc82_AxonDANCoord(maskedAxonDAN,1),nc82_AxonDANCoord(maskedAxonDAN,2),nc82_AxonDANCoord(maskedAxonDAN,3), 10, ceO, 'o')
view(3)
axis equal
xlim([0 50])
ylim([0 90])
zlim([0 30])
xticks([0:10:50])
yticks([0:10:100])
zticks([0:10:50])
xlabel('\mum')
grid on
f0 = gcf();
print(f0, '-painters','-dpdf', '-loose',[date ' Axon_Mask_nc82.pdf']);
