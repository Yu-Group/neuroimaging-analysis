rt = '/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/BleachingAnalysis/';
% rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/BleachingAnalysis/';
cd(rt)

LLSMOS = load([rt 'LLSM_bleachingPerTile.mat']);
Airy = load([rt 'Airy_bleachingPerTile.mat']);
SD = load([rt 'SpinningDisk_bleachingPerSlice.mat']);

% order LLSM sample sweep data
LLSM = load('LLSM_SampleSweep_BleachingParameters_20180125.mat', 'allParams');
LLSM_fileList = load('LLSM_SampleSweep_BleachingParameters_20180125.mat', 'fnlist');

mx = 11393; %1
my = 10663; %7
mz = 5301; % 41
csx = 5697;
csy = 1524;
csz = 130;
OL = 0;
s = [csx,csy,csz];

nx = ceil(mx/(s(1)-OL));
ny = ceil(my/(s(2)-OL));
nz = ceil(mz/(s(3)-OL));

%z: linear, from small to large
%y: meander, start with large to small (at the smallest z) and as z increases small to large, large to small, ...., etc.


% generate x y z cropping coordinates
xmin = zeros(nx*ny*nz,1);
xmax = zeros(nx*ny*nz,1);
ymin = zeros(nx*ny*nz,1);
ymax = zeros(nx*ny*nz,1);
zmin = zeros(nx*ny*nz,1);
zmax = zeros(nx*ny*nz,1);

nn = 1;
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            xmin(nn) = (1+((i-1)*(s(1)-OL)));
            xmax(nn) = (s(1)*i-(OL*(i-1)));
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            ymax(nn) = (s(2)*j-(OL*(j-1)));
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            zmax(nn) = (s(3)*k-(OL*(k-1)));
            
            nn = nn + 1;
        end
    end
end

clear AcqOrder
nn = 1;
for k = 1:nz
    if isEven(k)
        yr = 1:ny;
    else
        yr = ny:-1:1;
    end
    for j = yr
        for i = 1:nx
            xmin(nn) = (1+((i-1)*(s(1)-OL)));
            xmax(nn) = (s(1)*i-(OL*(i-1)));
            ymin(nn) = (1+((j-1)*(s(2)-OL)));
            ymax(nn) = (s(2)*j-(OL*(j-1)));
            zmin(nn) = (1+((k-1)*(s(3)-OL)));
            zmax(nn) = (s(3)*k-(OL*(k-1)));

            AcqOrder{nn} = [num2str(xmin(nn)) ',' num2str(ymin(nn)) ',' num2str(zmin(nn)) '_' num2str(xmax(nn)) ',' num2str(ymax(nn)) ',' num2str(zmax(nn)) '.tif'];
            nn = nn + 1;
        end
    end
end
AcqOrder = AcqOrder';

[~, Acqidx] = ismember(AcqOrder, LLSM_fileList.fnlist); 
LLSM.allParams = LLSM.allParams(Acqidx);
nonempidx = arrayfun(@(x) ~isempty(x.MaxI), LLSM.allParams);
LLSM.allParams = LLSM.allParams(nonempidx);
% %%
% figure, 
% for i = 1:nn-2
%     plot3(xmin(i:i+1), ymin(i:i+1),zmin(i:i+1))
% hold on
% pause
% end


%%

[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;
AirytileSize = [1948, 1948, 1525*9];
LLSM_OS_tileSize = [360, 704, 501*960]; % obj scan
LLSMtileSize = [5697, 1524, 130*numel(LLSM.allParams)];
SDtileSize = [3584, 3584, 1816];

AiryPx = [0.0517, 0.0517, 0.270]/1000; % in mm
SDPx = [0.160, 0.160, 0.250]/1000; % in mm
LLSM_OS_Px = [0.097, 0.097, 0.180]/1000; % in mm -- for obj scan
LLSMPx = [0.1229, 0.104, 0.16515]/1000; % in mm -- for sample sweep

minOcc = 0.005; % 0.5 percent of the voxels need to be occupied


Idx_Airy = arrayfun(@(x) (x.nVox)/(1948*1948), Airy.allParams, 'unif',0);
Idx_Airy = horzcat(Idx_Airy{:});

normInt_Airy = arrayfun(@(x) (x.SumI./x.nVox), Airy.allParams, 'unif',0);
normInt_Airy = horzcat(normInt_Airy{:});
normInt_Airy = normInt_Airy(Idx_Airy>minOcc);

sumInt_Airy = arrayfun(@(x) (x.SumI), Airy.allParams, 'unif',0);
sumInt_Airy = horzcat(sumInt_Airy{:});
sumInt_Airy = sumInt_Airy(Idx_Airy>minOcc);

% llsm - sample sweep
Idx_LLSM = arrayfun(@(x) (x.nVox)/(LLSMtileSize(1)*LLSMtileSize(2)), LLSM.allParams, 'unif',0);
Idx_LLSM = horzcat(Idx_LLSM{:});

normInt_LLSM = arrayfun(@(x) (x.SumI./x.nVox), LLSM.allParams, 'unif',0);
normInt_LLSM = horzcat(normInt_LLSM{:});
normInt_LLSM = normInt_LLSM(Idx_LLSM>minOcc);

sumInt_LLSM = arrayfun(@(x) (x.SumI), LLSM.allParams, 'unif',0);
sumInt_LLSM = horzcat(sumInt_LLSM{:});
sumInt_LLSM = sumInt_LLSM(Idx_LLSM>minOcc);

% llsm - obj scan
Idx_LLSMOS = arrayfun(@(x) (x.nVox)/(LLSM_OS_tileSize(1)*LLSM_OS_tileSize(2)), LLSMOS.allParams, 'unif',0);
Idx_LLSMOS = horzcat(Idx_LLSMOS{:});

normInt_LLSMOS = arrayfun(@(x) (x.SumI./x.nVox), LLSMOS.allParams, 'unif',0);
normInt_LLSMOS = horzcat(normInt_LLSMOS{:});
normInt_LLSMOS = normInt_LLSMOS(Idx_LLSMOS>minOcc);

sumInt_LLSMOS = arrayfun(@(x) (x.SumI), LLSMOS.allParams, 'unif',0);
sumInt_LLSMOS = horzcat(sumInt_LLSMOS{:});
sumInt_LLSMOS = sumInt_LLSMOS(Idx_LLSMOS>minOcc);

normInt_SD = SD.SumI./SD.nVox;
sumInt_SD = SD.SumI;




%% figures
bf = 2;

ha = setupFigure(1,2, 'AxesWidth', 4, 'AxesHeight', 4,'SameAxes', false,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1 1]);
set(ha, 'FontSize', 6);


% Airy

vals = normInt_Airy(~isnan(normInt_Airy));
x = [0:numel(vals)-1]*(AirytileSize(1) * AirytileSize (2) * prod(AiryPx(:))); 
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_Airy(~isnan(Int_Airy)), sf),99);
maxV = max(smooth(normInt_Airy(~isnan(normInt_Airy)), sf));
vals = smooth(normInt_Airy(~isnan(normInt_Airy)), sf)/maxV;
% 
axes(ha(1))
title('Sum / # Vox')
xlabel('Imaged Volume (mm^3)');
ylabel('Normalized Intensity');
plot(x, vals, 'color', ceP, 'Linewidth',2);

maxV = max(smooth(sumInt_Airy(~isnan(sumInt_Airy)), sf));
vals = smooth(sumInt_Airy(~isnan(sumInt_Airy)), sf)/maxV;
axes(ha(2))
plot(x, vals, 'color', ceP, 'Linewidth',2);

title('Sum')
xlabel('Imaged Volume (mm^3)');
ylabel('Normalized Intensity');


% LLSM - sample sweep
x = [0:LLSMtileSize(3)-1]*(LLSMtileSize(1) * LLSMtileSize(2) * prod(LLSMPx(:))); 
normInt_LLSM =normInt_LLSM(3500:end);
% normInt_LLSM = normInt_LLSM(find(vals == max(vals)):end);
% vals = normInt_LLSM(~isnan(normInt_LLSM));
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_LLSM(~isnan(Int_LLSM)), sf),99);
axes(ha(1))
maxV = max(smooth(normInt_LLSM(~isnan(normInt_LLSM)), sf));
vals = smooth(normInt_LLSM(~isnan(normInt_LLSM)), sf)/maxV;
vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSMtileSize(1) * LLSMtileSize(2) * prod(LLSMPx(:))); 
plot(x, vals, 'color', ceB, 'Linewidth',2)


axes(ha(2))
sumInt_LLSM = sumInt_LLSM(4000:end);
maxV = max(smooth(sumInt_LLSM(~isnan(sumInt_LLSM)), sf));
vals = smooth(sumInt_LLSM(~isnan(sumInt_LLSM)), sf)/maxV;
% vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSMtileSize(1) * LLSMtileSize(2) * prod(LLSMPx(:))); 
plot(x(1:end-1000), vals(1001:end), 'color', ceB, 'Linewidth',2)
hold on

% LLSM - objective scan
x = [0:LLSM_OS_tileSize(3)-1]*(LLSM_OS_tileSize(1) * LLSM_OS_tileSize(2) * prod(LLSM_OS_Px(:))); 
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_LLSM(~isnan(Int_LLSM)), sf),99);
axes(ha(1))
maxV = max(smooth(normInt_LLSMOS(~isnan(normInt_LLSMOS)), sf));
vals = smooth(normInt_LLSMOS(~isnan(normInt_LLSMOS)), sf)/maxV;
% vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSM_OS_tileSize(1) * LLSM_OS_tileSize(2) * prod(LLSM_OS_Px(:))); 
plot(x(1:end-400), vals(1:end-400), 'color', cfB, 'Linewidth',2)


axes(ha(2))
maxV = max(smooth(sumInt_LLSMOS(~isnan(sumInt_LLSMOS)), sf));
vals = smooth(sumInt_LLSMOS(~isnan(sumInt_LLSMOS)), sf)/maxV;
% vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSM_OS_tileSize(1) * LLSM_OS_tileSize (2) * prod(LLSMPx(:))); 
plot(x(1:end-400), vals(1:end-400), 'color', cfB, 'Linewidth',2)
hold on


% SD
x = [0:SDtileSize(3)-1]*(SDtileSize(1) * SDtileSize (2) * prod(SDPx(:))); 
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_SD, sf'),99);
maxV = max(smooth(normInt_SD, sf'));
vals = smooth(normInt_SD, sf)/maxV;

axes(ha(1))
plot(x, vals, 'color', ceG, 'Linewidth',2)

legend('Airyscan', 'LLSM-SS', 'LLSM-OS','SD')

ylim([0 1])
xlim([0 0.25])
maxV = max(smooth(sumInt_SD, sf));
vals = smooth(sumInt_SD, sf)/maxV;

axes(ha(2))
plot(x, vals, 'color', ceG, 'Linewidth',2)

ylim([0 1])
xlim([0 0.25])

f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[rt date 'withLLSM_OS_BleachingComparison_0_0p25.eps']);
%% 3d plot
% time to image in minutes
t_LLSM_O = 175;
t_LLSM_SS = 130;
t_SD = 237;
t_AS = 441;



% figures
bf = 2;

ha = setupFigure(1,2, 'AxesWidth', 4, 'AxesHeight', 4,'SameAxes', false,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1 1]);
set(ha, 'FontSize', 6);


% Airy

vals = normInt_Airy(~isnan(normInt_Airy));
x = [0:numel(vals)-1]*(AirytileSize(1) * AirytileSize (2) * prod(AiryPx(:))); 
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_Airy(~isnan(Int_Airy)), sf),99);
maxV = max(smooth(normInt_Airy(~isnan(normInt_Airy)), sf));
vals = smooth(normInt_Airy(~isnan(normInt_Airy)), sf)/maxV;
% 
axes(ha(1))
title('Sum / # Vox')
xlabel('Imaged Volume (mm^3)');
ylabel('time (min)');
zlabel('Normalized Intensity');
t = 0:t_AS/numel(x):t_AS;
plot3(x,t(1:end-1) ,vals, 'color', ceP, 'Linewidth',2);
view([-30, 34])
maxV = max(smooth(sumInt_Airy(~isnan(sumInt_Airy)), sf));
vals = smooth(sumInt_Airy(~isnan(sumInt_Airy)), sf)/maxV;
axes(ha(2))
t = 0:t_AS/numel(x):t_AS;
plot3(x,t(1:end-1) ,vals,  'color', ceP, 'Linewidth',2);

title('Sum')
xlabel('Imaged Volume (mm^3)');
ylabel('time (min)');
zlabel('Normalized Intensity');
view([-30, 34])
%

% LLSM - sample sweep
x = [0:LLSMtileSize(3)-1]*(LLSMtileSize(1) * LLSMtileSize(2) * prod(LLSMPx(:))); 
normInt_LLSM =normInt_LLSM(3500:end);
% normInt_LLSM = normInt_LLSM(find(vals == max(vals)):end);
% vals = normInt_LLSM(~isnan(normInt_LLSM));
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_LLSM(~isnan(Int_LLSM)), sf),99);
axes(ha(1))
maxV = max(smooth(normInt_LLSM(~isnan(normInt_LLSM)), sf));
vals = smooth(normInt_LLSM(~isnan(normInt_LLSM)), sf)/maxV;
vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSMtileSize(1) * LLSMtileSize(2) * prod(LLSMPx(:))); 
t = 0:t_LLSM_SS/numel(x):t_LLSM_SS;
plot3(x,t(1:end-1) ,vals,  'color', ceB, 'Linewidth',2)


axes(ha(2))
sumInt_LLSM = sumInt_LLSM(4000:end);
maxV = max(smooth(sumInt_LLSM(~isnan(sumInt_LLSM)), sf));
vals = smooth(sumInt_LLSM(~isnan(sumInt_LLSM)), sf)/maxV;
% vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSMtileSize(1) * LLSMtileSize(2) * prod(LLSMPx(:))); 
t = 0:t_LLSM_SS/numel(x):t_LLSM_SS;
plot3(x(1:end-1000), t(1:end-1001),vals(1001:end), 'color', ceB, 'Linewidth',2)
hold on

% LLSM - objective scan
x = [0:LLSM_OS_tileSize(3)-1]*(LLSM_OS_tileSize(1) * LLSM_OS_tileSize(2) * prod(LLSM_OS_Px(:))); 
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_LLSM(~isnan(Int_LLSM)), sf),99);
axes(ha(1))
maxV = max(smooth(normInt_LLSMOS(~isnan(normInt_LLSMOS)), sf));
vals = smooth(normInt_LLSMOS(~isnan(normInt_LLSMOS)), sf)/maxV;
% vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSM_OS_tileSize(1) * LLSM_OS_tileSize(2) * prod(LLSM_OS_Px(:))); 
t = 0:t_LLSM_O/numel(x):t_LLSM_O;
plot3(x(1:end-400), t(1:end-401),vals(1:end-400), 'color', cfB, 'Linewidth',2)


axes(ha(2))
maxV = max(smooth(sumInt_LLSMOS(~isnan(sumInt_LLSMOS)), sf));
vals = smooth(sumInt_LLSMOS(~isnan(sumInt_LLSMOS)), sf)/maxV;
% vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSM_OS_tileSize(1) * LLSM_OS_tileSize (2) * prod(LLSMPx(:))); 
t = 0:t_LLSM_O/numel(x):t_LLSM_O;
plot3(x(1:end-400), t(1:end-401),vals(1:end-400), 'color', cfB, 'Linewidth',2)
hold on


% SD
x = [0:SDtileSize(3)-1]*(SDtileSize(1) * SDtileSize (2) * prod(SDPx(:))); 
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_SD, sf'),99);
maxV = max(smooth(normInt_SD, sf'));
vals = smooth(normInt_SD, sf)/maxV;

axes(ha(1))
t = 0:t_SD/numel(x):t_SD;
plot3(x, t(1:end-1), vals, 'color', ceG, 'Linewidth',2)
% view(3)
legend('Airyscan', 'LLSM-SS', 'LLSM-OS','SD')

ylim([0 1])
xlim([0 0.25])
maxV = max(smooth(sumInt_SD, sf));
vals = smooth(sumInt_SD, sf)/maxV;

axes(ha(2))
t = 0:t_SD/numel(x):t_SD;
plot3(x, t(1:end-1),vals, 'color', ceG, 'Linewidth',2)
% view(3)
ylim([0 1])
xlim([0 0.25])

f0 = gcf();
% print(f0, '-painters','-depsc', '-loose',[rt date 'withLLSM_OS_BleachingComparison_0_0p25.eps']);


%% 2 yaxis plot
% time to image in minutes
t_LLSM_O = 175;
t_LLSM_SS = 130;
t_SD = 237;
t_AS = 441;



% figures
bf = 2;

ha = setupFigure(1,2, 'AxesWidth', 4, 'AxesHeight', 4,'SameAxes', false,...
    'XSpace', [1.5 1.5 1.5], 'YSpace', [1.5 1 1]);
set(ha, 'FontSize', 6);


% Airy

vals = normInt_Airy(~isnan(normInt_Airy));
x = [0:numel(vals)-1]*(AirytileSize(1) * AirytileSize (2) * prod(AiryPx(:))); 
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_Airy(~isnan(Int_Airy)), sf),99);
maxV = max(smooth(normInt_Airy(~isnan(normInt_Airy)), sf));
vals = smooth(normInt_Airy(~isnan(normInt_Airy)), sf)/maxV;
% 
axes(ha(1))
title('Sum / # Vox')
xlabel('Imaged Volume (mm^3)');
yyaxis right
ylabel('time (min)');
yyaxis right
ylabel('Normalized Intensity');
t = 0:t_AS/numel(x):t_AS;
yyaxis left
plot(x, vals, 'color', ceP, 'Linewidth',2);
yyaxis right
plot(x,t(1:end-1), '--','color', ceP, 'Linewidth',2);

maxV = max(smooth(sumInt_Airy(~isnan(sumInt_Airy)), sf));
vals = smooth(sumInt_Airy(~isnan(sumInt_Airy)), sf)/maxV;
axes(ha(2))
t = 0:t_AS/numel(x):t_AS;
yyaxis left
plot(x, vals, '-', 'color', ceP, 'Linewidth',2);
yyaxis right
plot(x,t(1:end-1), '--', 'color', ceP, 'Linewidth',2);

title('Sum')
xlabel('Imaged Volume (mm^3)');
yyaxis right
ylabel('time (min)');
yyaxis right
ylabel('Normalized Intensity');

% LLSM - sample sweep
x = [0:LLSMtileSize(3)-1]*(LLSMtileSize(1) * LLSMtileSize(2) * prod(LLSMPx(:))); 
normInt_LLSM =normInt_LLSM(3500:end);
% normInt_LLSM = normInt_LLSM(find(vals == max(vals)):end);
% vals = normInt_LLSM(~isnan(normInt_LLSM));
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_LLSM(~isnan(Int_LLSM)), sf),99);
axes(ha(1))
maxV = max(smooth(normInt_LLSM(~isnan(normInt_LLSM)), sf));
vals = smooth(normInt_LLSM(~isnan(normInt_LLSM)), sf)/maxV;
vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSMtileSize(1) * LLSMtileSize(2) * prod(LLSMPx(:))); 
t = 0:t_LLSM_SS/numel(x):t_LLSM_SS;
yyaxis left
plot(x, vals, '-', 'color', ceB, 'Linewidth',2)
yyaxis right
plot(x,t(1:end-1), '--', 'color', ceB, 'Linewidth',2)

axes(ha(2))
sumInt_LLSM = sumInt_LLSM(4000:end);
maxV = max(smooth(sumInt_LLSM(~isnan(sumInt_LLSM)), sf));
vals = smooth(sumInt_LLSM(~isnan(sumInt_LLSM)), sf)/maxV;
% vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSMtileSize(1) * LLSMtileSize(2) * prod(LLSMPx(:))); 
t = 0:t_LLSM_SS/numel(x):t_LLSM_SS;
yyaxis left
plot(x(1:end-1000), vals(1001:end), '-', 'color', ceB, 'Linewidth',2)
yyaxis right
plot(x(1:end-1000), t(1:end-1001), '--', 'color', ceB, 'Linewidth',2)
hold on
%
% LLSM - objective scan
x = [0:LLSM_OS_tileSize(3)-1]*(LLSM_OS_tileSize(1) * LLSM_OS_tileSize(2) * prod(LLSM_OS_Px(:))); 
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_LLSM(~isnan(Int_LLSM)), sf),99);
axes(ha(1))
maxV = max(smooth(normInt_LLSMOS(~isnan(normInt_LLSMOS)), sf));
vals = smooth(normInt_LLSMOS(~isnan(normInt_LLSMOS)), sf)/maxV;
% vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSM_OS_tileSize(1) * LLSM_OS_tileSize(2) * prod(LLSM_OS_Px(:))); 
t = 0:t_LLSM_O/numel(x):t_LLSM_O;
yyaxis left
plot(x(1:end-400), vals(1:end-400), '-','color', cfB, 'Linewidth',2)
yyaxis right
plot(x(1:end-400), t(1:end-401), '--','color', cfB, 'Linewidth',2)


axes(ha(2))
maxV = max(smooth(sumInt_LLSMOS(~isnan(sumInt_LLSMOS)), sf));
vals = smooth(sumInt_LLSMOS(~isnan(sumInt_LLSMOS)), sf)/maxV;
% vals = vals(find(vals == max(vals)):end);
x = [0:numel(vals)-1]*(LLSM_OS_tileSize(1) * LLSM_OS_tileSize (2) * prod(LLSMPx(:))); 
t = 0:t_LLSM_O/numel(x):t_LLSM_O;
yyaxis left
plot(x(1:end-400), vals(1:end-400), '-', 'color', cfB, 'Linewidth',2)
yyaxis right
plot(x(1:end-400), t(1:end-401), '--', 'color', cfB, 'Linewidth',2)

hold on


% SD
x = [0:SDtileSize(3)-1]*(SDtileSize(1) * SDtileSize (2) * prod(SDPx(:))); 
sf = round(numel(x)/bf);
% maxV = prctile(smooth(Int_SD, sf'),99);
maxV = max(smooth(normInt_SD, sf'));
vals = smooth(normInt_SD, sf)/maxV;

axes(ha(1))
t = 0:t_SD/numel(x):t_SD;
yyaxis left
plot(x, vals, '-','color', ceG, 'Linewidth',2)
ylim([0 1])
yyaxis right
plot(x, t(1:end-1), '--', 'color', ceG, 'Linewidth',2)
% ylim([0 1])
% legend('Airyscan', 'LLSM-SS', 'LLSM-OS','SD')


xlim([0 0.25])
maxV = max(smooth(sumInt_SD, sf));
vals = smooth(sumInt_SD, sf)/maxV;

axes(ha(2))
t = 0:t_SD/numel(x):t_SD;
yyaxis left
plot(x, vals, '-','color', ceG, 'Linewidth',2)
ylim([0 1])
yyaxis right
plot(x, t(1:end-1), '--','color', ceG, 'Linewidth',2)
% ylim([0 1])

xlim([0 0.25])

f0 = gcf();
print(f0, '-painters','-depsc', '-loose',[rt date 'withLLSM_OS_BleachingComparison_0_0p25_withTIME.eps']);
