cd('/Users/GU/Dropbox/Manuscript_ExLLSM/GokulWorkspace/unmixing');
% im0 = readtiff('UnmixingEx3_channel 0 [4647, 2041, 2633]_300cube.tif');
% im1 = readtiff('UnmixingEx3_channel 1 [4647, 2041, 2633]_300cube.tif');

im0 = readtiff('UnmixingEx2_channel 0 [4831, 1842, 2740].tif');
im1 = readtiff('UnmixingEx2_channel 1 [4831, 1842, 2740].tif');
%% excitation
ch0 = 41.70;
ch1 = 26.76;
rch = ch1/ch0;
%%
Bk0 = 460 + 12*2;%thresholdOtsu(im0(:));%
Bk1 = 3080 + 80*2;%thresholdOtsu(im1(:));%
im0_bk = im0 - Bk0;
im1_bk = im1 - Bk1;

%%
r = im1_bk(im1_bk>0 & im0_bk >0)./im0_bk(im1_bk>0 & im0_bk >0);
figure,
histogram(r,3000)
xlim([0 40])
%%
cim0 = im0_bk-(im1_bk*rch/3) ;

writetiff(cim0, 'unmixed_UnmixingEx2_channel 0 [4831, 1842, 2740].tif')

%%
rt0 = '/groups/betzig/betziglab/4Stephan/160404_Sample3_PM_Austin/Position2_Tile/Stitch_Igor/restitching/restitching-covariance/slice-tiff/ch0-final-flatfield-decon/';
rt1 = '/groups/betzig/betziglab/4Stephan/160404_Sample3_PM_Austin/Position2_Tile/Stitch_Igor/restitching/restitching-covariance/slice-tiff/ch1-final-flatfield-decon/';

fn0 = dir([rt0 '*.tif']);
fn0 = {fn0.name};

fn1 = dir([rt1 '*.tif']);
fn1 = {fn1.name};
% GU_LinearUnmixing(im0p, im1p, ch0exMax, ch1exMax, Bk0, Bk1, SF)
ch0exMax = 42;
ch1exMax = 27;
Bk0 = 460 + 12*2;%thresholdOtsu(im0(:));%
Bk1 = 3080 + 80*2;%thresholdOtsu(im1(:));%
SF = 3;
for i = 1:numel(fn0)
    im0p = [rt0 fn0{i}];
    im1p = [rt1 fn1{i}];
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
    
end

%% for figure 1
%% LLSM RAW data

rt0 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/raw_export/slice-tiff/ch0_mirroredX-final-raw/';
rt1 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/raw_export/slice-tiff/ch1_mirroredX-final-raw/';

fn0 = dir([rt0 '*.tif']);
fn0 = {fn0.name};

fn1 = dir([rt1 '*.tif']);
fn1 = {fn1.name};

ch0exMax = 38; % eyfp
ch1exMax = 5; % alexa 568
Bk0 = 104 + 4*2;%thresholdOtsu(im0(:));%
Bk1 = 112 + 4*2;%thresholdOtsu(im1(:));%
SF = 1;

for i = 1:numel(fn0)
    im0p = [rt0 fn0{i}];
    im1p = [rt1 fn1{i}];
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
end
%% LLSM RAW data - resub missing jobs
rtpre = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/raw_export/slice-tiff/ch0_mirroredX-final-raw/';
rt = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/raw_export/slice-tiff/ch0_mirroredX-final-raw/unmixed/';
cd(rt)
fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name};

fnlistpre = dir([rtpre '*.tif']);
fnlistpre = {fnlistpre.name};
CompletedJobsIdx = zeros(1, numel(fnlistpre));
for ex = 1:numel(fnlistpre)
    fn = ['unmixed_' num2str(ex-1) '.tif'];
    CompletedJobsIdx(ex) = ismember(fn,fnlist);%exist(fn(2:end), 'file');
end

job2resubmit = find(CompletedJobsIdx==0);


rt0 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/raw_export/slice-tiff/ch0_mirroredX-final-raw/';
rt1 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/raw_export/slice-tiff/ch1_mirroredX-final-raw/';

fn0 = dir([rt0 '*.tif']);
fn0 = {fn0.name};

fn1 = dir([rt1 '*.tif']);
fn1 = {fn1.name};

ch0exMax = 38; % eyfp
ch1exMax = 5; % alexa 568
Bk0 = 104 + 4*2;%thresholdOtsu(im0(:));%
Bk1 = 112 + 4*2;%thresholdOtsu(im1(:));%
SF = 1;

for i = job2resubmit
    im0p = [rt0 num2str(i-1) '.tif'];
    im1p = [rt1 num2str(i-1) '.tif'];
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
end
%% LLSM DECON data

rt0 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/slice-tiff/ch0/';
rt1 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/slice-tiff/ch1/';

fn0 = dir([rt0 '*.tif']);
fn0 = {fn0.name};

fn1 = dir([rt1 '*.tif']);
fn1 = {fn1.name};

ch0exMax = 38; % eyfp
ch1exMax = 5; % alexa 568
Bk0 = 1500 + 43*2;%thresholdOtsu(im0(:));%
Bk1 = 217 + 13*2;%thresholdOtsu(im1(:));%
SF = 0.1;

for i = 1:numel(fn0)
    im0p = [rt0 fn0{i}];
    im1p = [rt1 fn1{i}];
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
end

%% LLSM DECON data - resub missing jobs
rtpre = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/slice-tiff/ch0/';
rt = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/slice-tiff/ch0/unmixed/';

fnlist = dir([rt '*.tif']);
fnlist = {fnlist.name};

fnlistpre = dir([rtpre '*.tif']);
fnlistpre = {fnlistpre.name};
CompletedJobsIdx = zeros(1, numel(fnlistpre));
for ex = 1:numel(fnlistpre)
    fn = ['unmixed_' num2str(ex-1) '.tif'];
    CompletedJobsIdx(ex) = ismember(fn,fnlist);%exist(fn(2:end), 'file');
end

job2resubmit = find(CompletedJobsIdx==0);

rt0 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/slice-tiff/ch0/';
rt1 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/Stitch_Igor/slice-tiff/ch1/';

fn0 = dir([rt0 '*.tif']);
fn0 = {fn0.name};

fn1 = dir([rt1 '*.tif']);
fn1 = {fn1.name};

ch0exMax = 38; % eyfp
ch1exMax = 5; % alexa 568
Bk0 = 1500 + 43*2;%thresholdOtsu(im0(:));%
Bk1 = 217 + 13*2;%thresholdOtsu(im1(:));%
SF = 0.1;

for i = job2resubmit
    im0p = [rt0 num2str(i-1) '.tif'];
    im1p = [rt1 num2str(i-1) '.tif'];
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
end
%% Spinning disc - decon

rt0 = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Spinningdisk/Flatfield+deconed/488/';
rt1 = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Spinningdisk/Flatfield+deconed/560/';

fn0 = dir([rt0 '*.tif']);
fn0 = {fn0.name};

fn1 = dir([rt1 '*.tif']);
fn1 = {fn1.name};

ch0exMax = 38; % eyfp
ch1exMax = 5; % alexa 568
Bk0 = 12  + 32*2;%thresholdOtsu(im0(:));%
Bk1 = 100 + 63*2;%thresholdOtsu(im1(:));%
SF = 0.5;

for i = 1:numel(fn0)
    im0p = [rt0 fn0{i}];
    im1p = [rt1 fn1{i}];
    submit = ['bsub -n 1 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
end

%% Spinning disc - raw

rt0 = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Spinningdisk/Raw/488/';
rt1 = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Spinningdisk/Raw/560/';

fn0 = dir([rt0 '*.tif']);
fn0 = {fn0.name};

fn1 = dir([rt1 '*.tif']);
fn1 = {fn1.name};

ch0exMax = 38; % eyfp
ch1exMax = 5; % alexa 568
Bk0 = 107  + 4*2;%thresholdOtsu(im0(:));%
Bk1 = 126 + 9*2;%thresholdOtsu(im1(:));%
SF = 1;

for i = 1:numel(fn0)
    im0p = [rt0 fn0{i}];
    im1p = [rt1 fn1{i}];
    submit = ['bsub -n 1 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
end

%% Airy Scan - decon

rt0 = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/488/';
rt1 = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/560/';

fn0 = dir([rt0 '*.tif']);
fn0 = {fn0.name};

fn1 = dir([rt1 '*_c1']);
fn1 = {fn1.name};

ch0exMax = 38; % eyfp
ch1exMax = 5; % alexa 568
Bk0 = 0  + 1*2;%thresholdOtsu(im0(:));%
Bk1 = 0 + 1*2;%thresholdOtsu(im1(:));%
SF = 0.5;

for i = 1:numel(fn0)
    im0p = [rt0 fn0{i}];
    im1p = [rt1 fn1{i}];
    submit = ['bsub -n 1 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
end


%% Airy Scan - decon - inividual tiles

rt0 = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/Individualtiles/488/';
rt1 = '/groups/betzig/betziglab/4Stephan/171105_unmixing/Airyscan/Individualtiles/560/';

fn0 = dir([rt0 '*.tif']);
fn0 = {fn0.name};

fn1 = dir([rt1 '*.tif']);
fn1 = {fn1.name};

ch0exMax = 38; % eyfp
ch1exMax = 5; % alexa 568
Bk0 = 0  + 1*2;%thresholdOtsu(im0(:));%
Bk1 = 0 + 1*2;%thresholdOtsu(im1(:));%
SF = 0.5;

for i = 1:numel(fn0)
    im0p = [rt0 fn0{i}];
    im1p = [rt1 fn1{i}];
    submit = ['bsub -n 16 -R"affinity[core(1)]"  -R"select[broadwell]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
end

%% LLSM - decon only tiles


rt0 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/matlab_decon/';
rt1 = '/groups/betzig/betziglab/Ved/Data/20161211_ExM/YFP2_comparison2/matlab_decon/';

fn0 = dir([rt0 '*488nm*.tif']);
fn0 = {fn0.name}';

fn1 = dir([rt1 '*560nm*.tif']);
fn1 = {fn1.name}';

ch0exMax = 38; % eyfp
ch1exMax = 5; % alexa 568
Bk0 = 1500 + 43*2;%thresholdOtsu(im0(:));%
Bk1 = 217 + 13*2;%thresholdOtsu(im1(:));%
SF = 0.1;

for i = 1:numel(fn0)
    im0p = [rt0 fn0{i}];
    im1p = [rt1 fn1{i}];
    submit = ['bsub -n 1 -R"affinity[core(1)]" -J' ' "unmix' num2str(i) '" -o ~/logs/unmix' num2str(i) ' ''export MCR_CACHE_ROOT=/scratch/upadhyayulas/mcr_cache_root.'  num2str(i) ' ; /groups/betzig/home/upadhyayulas/GU_LinearUnmixing ' im0p ' ' im1p ' ' num2str(ch0exMax) ' ' num2str(ch1exMax) ' ' num2str(Bk0) ' ' num2str(Bk1) ' ' num2str(SF) ''''];
    [~, ~] = system(submit,'-echo');
end