rt = '/Volumes/D_36win/GoogleDrive_DataTransfer/JaneliaSharedwithGokul/ExpansionRelated/20161211_ExM_YFP2_threeColor_z0p18/';
fn{1} = 'c1.tif';
fn{2} = 'c2.tif';

im1 = readtiff([rt fn{1}]);
im2 = readtiff([rt fn{2}]);

%%

T1 = thresholdOtsu(im1(im1>0));
T2 = thresholdOtsu(im2(im2>0));
im1t = im1-T1;
im2t = im2-T2;
%%
r = 500:700;
imtool(im2t(r,r,200)-im1t(r,r,200)*1000)
% imtool(im2t(r,r,200))
%%
im2c = im2t-im1t*1000;

writetiff(im2c, 'c2_noc1.tif');