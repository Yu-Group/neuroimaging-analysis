function coord = GU_ExM_getPNMaskCoord(rt)


mz = 5186;

coord(1:mz) = struct('x',[], 'y', [], 'z', [], 'label', []);
FE = zeros(1, mz);
parfor_progress(mz);
parfor j = 1:mz
    if exist([rt num2str(j) '.nrrd.tif'],'file')
        [mask] = readtiff([rt num2str(j) '.nrrd.tif']);
        u = unique(mask(mask>0));
        if ~isempty(u)
            for i = 1:numel(u)
                [y,x,~] = ind2sub(size(mask), find(mask == u(i)));
                coord(j).x{i} = x;
                coord(j).y{i} = y;
                coord(j).z{i} = repmat(j, numel(x),1);
                coord(j).label{i} = repmat(u(i),numel(x),1);
            end
        end
        FE(j) = 1;
    end
    parfor_progress;
end
parfor_progress(0);
save([rt 'data_cluster.mat'], 'coord', 'FE');
