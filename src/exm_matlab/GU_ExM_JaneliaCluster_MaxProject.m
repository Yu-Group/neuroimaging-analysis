function GU_ExM_JaneliaCluster_MaxProject(minfn, maxfn, rt, srt, sfn)
% Gokul Upadhyayula, Oct 2018

if ~exist(srt, 'dir')
    mkdir(srt)
end

minfn = str2double(minfn);
maxfn = str2double(maxfn);

im = readtiff([rt  num2str(minfn) '.tif']);
maxim = im;

for n = (minfn+1):maxfn
    im = readtiff([rt num2str(n) '.tif']);
    maxim = max(im, maxim);
end
writetiff(maxim, [srt  sfn]);