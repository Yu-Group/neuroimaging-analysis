function GU_ExM_JaneliaCluster_SpectralMixMouseSomato(ex)

% this function combines volumes in different regions of the 16bit volume
% this assumes that the volumes are spatially non-overlapping and max projected
% Gokul Upadhyayula, Oct 2017

if ischar(ex)
    ex = str2double(ex);
end

RangeIn = [0 32760];
HomerRange = 32768;

rtSM = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed_cubicInterp_Clean_500_3200/';
rt = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/ch0/Analysis_1ch/1008/RotatedStacks/nonEmpty/slice-tiff/ch0/';
rtsave = '/groups/betzig/betziglab/4Stephan/171104_Mousebrainsynapse/SpectralMixed/LinearYFP_cubicHomer_SM/';

imSM = readtiff([rtSM num2str(ex) '.tif']);
imyfp = readtiff([rt num2str(ex) '.tif']);


imyfp(imyfp > RangeIn(2)) = RangeIn(2);
imyfp(imSM>HomerRange) = imSM(imSM>HomerRange);

writetiff(imyfp, [rtsave num2str(ex) '.tif']);
