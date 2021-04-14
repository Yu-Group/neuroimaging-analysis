    function im2 = GU_ExM_zOffset3Dframe(Path, Top,nF)
        im = readtiff(Path);
        sizeIM = size(im);
        im2 = zeros([sizeIM]+[0 0 nF]);
        if Top
            im2(:,:,nF+1:end) = im;
        else
            im2(:,:,1:end-nF) = im;
        end
        im2 = im2(:,:,1:end-nF); % crop to the orginal vol by removing end frames
        writetiff(uint16(im2),Path);
    end