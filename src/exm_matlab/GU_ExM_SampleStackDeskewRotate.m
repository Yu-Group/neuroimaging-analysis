function GU_ExM_SampleStackDeskewRotate(volpath, dz, varargin)

% sample scan deskew, rotate, deconvolve pipeline
% rotated data are resampled
% Gokul Upadhyayula, Nov 2017

[rt, fn, ext] = fileparts(volpath);
% opts = load([rt filesep 'deskewSettings.mat'], 'opts');
% varargin = {opts.opts{:}};
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('volpath'); % data structure from loadConditionData
ip.addRequired('dz');
ip.addParameter('PSFpath', '');
ip.addParameter('Crop', 1);
ip.addParameter('xyPixelSize', 104); % in nm
% deskew options
ip.addParameter('Deskew', 1);
ip.addParameter('SkewAngle', 32.2);
ip.addParameter('Rotate', 1);
ip.addParameter('Reversed', 1);
% deconvovle
ip.addParameter('Deconvovle', 1);
ip.addParameter('dzPSF', 0.1);
ip.addParameter('Background', 100);
ip.addParameter('nIter', 15);
ip.addParameter('Save16bit', 0); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('EdgeArtifacts', 0); % sets 20 pixels to zero to remove
ip.parse(volpath, dz, varargin{:});
p = ip.Results;

if ischar(p.xyPixelSize)
    p.xyPixelSize = str2double(p.xyPixelSize);
end

if ischar(p.dz)
    p.dz = str2double(p.dz);
end

if ischar(p.SkewAngle)
    p.SkewAngle = str2double(p.SkewAngle);
end

if ischar(p.Deskew)
    p.Deskew = logical(strcmp(p.Deskew, '1'));
end

if ischar(p.Crop)
    p.Crop = logical(strcmp(p.Crop, '1'));
end

if ischar(p.Rotate)
    p.Rotate = logical(strcmp(p.Rotate, '1'));
end

if ischar(p.Save16bit)
    p.Save16bit = logical(strcmp(p.Save16bit, '1'));
end


if ischar(p.Reversed)
    p.Reversed = logical(strcmp(p.Reversed, '1'));
end

if ischar(p.Deconvovle)
    p.Deconvovle = logical(strcmp(p.Deconvovle, '1'));
end

if ischar(p.dzPSF)
    p.dzPSF = str2double(p.dzPSF);
end
if ischar(p.Background)
    p.Background = str2double(p.Background);
end
if ischar(p.nIter)
    p.nIter = str2double(p.nIter);
end

theta = p.SkewAngle;
volout = readtiff(volpath);
zAniso = p.dz*sind(theta)/p.xyPixelSize;

zf = cotd(abs(theta));
yf = 1;
xf = cosd(abs(theta)) + tand(abs(theta))*sind(abs(theta));


% deskew
if p.Deskew
    deskewpath = [rt filesep 'DS'];
    if ~exist(deskewpath,'dir')
        mkdir(deskewpath);
    end
    
    fprintf('Deskewing Data...')
    tic
    volout = deskewFrame3D(volout, theta, p.dz, p.xyPixelSize,'Crop', logical(p.Crop),'reverse', logical(p.Reversed)); toc
    if p.Save16bit
        writetiff(uint16(volout), [deskewpath filesep fn ext]);
    else
        writetiff(volout, [deskewpath filesep fn ext]);
    end
end
mask = logical(volout);

SV = sum(volout(:));

% deconvovle
if p.Deconvovle
    deconpath = [rt filesep 'matlab_decon'];
    if ~exist(deconpath, 'dir')
        mkdir(deconpath);
    end
    
    fprintf('Deconvolving Data...')
    tic
    if SV > 0
        tic
        volout = RLdecon_for_ExMPipeline(volout, volpath, p.PSFpath, p.Background, p.nIter, p.dzPSF, zAniso, 0, [], 0, ...
            p.xyPixelSize, 0, 0, 0, 0, [0,0,1],0, []) ;toc
        
        % remove edge artifacts
        fprintf('Removing Edges...'); tic
        r = p.EdgeArtifacts;
        if r > 0
            se = strel('disk', r);
            mask = imerode(mask,se);
            mask(:,:,1:r) = 0;
            mask(:,:,end-r:end) = 0;
            volout(~mask) = 0;
            clear mask
            toc
        end
    end
    toc
    if p.Save16bit
        writetiff(uint16(volout), [deconpath filesep fn  ext]);
    else
        writetiff(volout, [deconpath filesep fn ext]);
    end
end

% rotate
if p.Rotate
    dsrpath = [rt filesep 'DSR_dx' num2str(p.xyPixelSize*xf) '_dz' num2str(p.xyPixelSize*zf)];
    if ~exist(dsrpath, 'dir')
        mkdir(dsrpath);
    end
    
    fprintf('Rotating Data...')
    tic
    volout = rotateFrame3D(volout, theta, zAniso, 'reverse', logical(p.Reversed),...
        'Crop', logical(p.Crop), 'ObjectiveScan', false); toc
    
    fprintf('resampling Rotated Data...')
    tic
    volout = GU_resampleStack3D(volout, xf, yf, zf,'Interp', 'linear');toc
end

if p.Save16bit
    writetiff(uint16(volout), [dsrpath filesep fn  ext]);
else
    writetiff(volout, [dsrpath filesep fn ext]);
end


