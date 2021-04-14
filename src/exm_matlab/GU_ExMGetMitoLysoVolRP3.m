function GU_ExMGetMitoLysoVolRP3(data)
i = 1
%    l(i) = load([data(i).source 'parameters_rp3_lyso_mito_all.mat']); 
    lysoMask = readtiff([data(i).source 'Analysis' filesep 'LysoMask_new.tif']);
    mitoMask = readtiff([data(i).source  'Analysis' filesep 'MitoMask_new.tif']);

    [lholes, N] = bwlabeln(mitoMask);
    CC = bwconncomp(lholes, 26);
    csize = cellfun(@numel, CC.PixelIdxList);
    mitoholeVol = csize*data(i).pixelSize^2*data(i).dz;%%%%%%% calc volume of mito holes
    mitorp3 = regionprops3(logical(mitoMask));
    
    clear mitoMask lholes

    [lholes, N]  = bwlabeln(lysoMask);
    ucl = unique(lholes);
    ucl = ucl(ucl>0);
    CC = bwconncomp(lholes, 26);
    csize = cellfun(@numel, CC.PixelIdxList);
    lysoholeVol = csize*data(i).pixelSize^2*data(i).dz;%%%%%%% calc volume of lyso holes
    lysorp3 = regionprops3(logical(lysoMask));
    
    clear lysoMask lholes
    
    save([data(i).source  'parameters_rp3_lyso_mito_all'], 'lysoholeVol', 'lysorp3', 'mitorp3', 'mitoholeVol')

