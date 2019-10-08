function [rgb_out, mask_out, stats] = IHC_GFAP_seg(fname, varargin) 

%% Parse inputs
p = inputParser;
% defaultRadiusList = [14,22];
% defaultNucSaturationThresh = 0.65;
% defaultCircThresh = 0.8;
defaultSmoothingKernelSize = 10;
defaultDownsampleFactor = 4;
defaultRefImage = nan;
addRequired(p,'fname');
addOptional(p,'SmoothingKernelSize',defaultSmoothingKernelSize);
addOptional(p,'DownsampleFactor',defaultDownsampleFactor);
addOptional(p,'RefImage',defaultRefImage);
parse(p,'fname',varargin{:});
SmoothingKernelSize = p.Results.SmoothingKernelSize;
downsample_factor = p.Results.DownsampleFactor;
refimg = p.Results.RefImage;

%% Load input file
rgb_in = imread(fname); 

%% Downsample
rgb = imresize(rgb_in, (1/downsample_factor));
rgbsm = imfilter(rgb,fspecial('gaussian',8*(4/downsample_factor),8*(4/downsample_factor)),'conv');

% RGB = mat2gray(rgbsmn);
RGB = mat2gray(rgbsm);
LAB = rgb2lab(RGB);

a = LAB(:,:,2); %an = (a - mean(a(:)))/std(a(:));
b = LAB(:,:,3); %bn = (b - mean(b(:)))/std(b(:));

A_center = 20;
B_center = 28;

d = ((1-(a - A_center)).^2 + (b - B_center).^2).^0.5;

Thresh = 25;

mask_TriB = d < Thresh;

rgb_out = rgb;
mask_out = mask_TriB;

stats.GFAPArea = bwarea(mask_out)*downsample_factor;
% stats.VesselLumenArea = bwarea(mask_Lumen_out)*downsample_factor;
% stats.VesselLumenPerim = sum(sum(bwperim(mask_Lumen_out)))*downsample_factor;


end

function totallen = contourlen(C)
k = 0;
    for i = 1:size(C,1)-1
        dx = C(i+1,1) - C(i,1);
        dy = C(i+1,2) - C(i,2);
        k = k + sqrt(dx^2 + dy^2);
    end
totallen = k;
end

function theta_t = contourangle(C)
k = 0;
    for i = 1:size(C,1)-1
        dx = C(i+1,1) - C(i,1);
        dy = C(i+1,2) - C(i,2);
        k = k + abs(atan(dy/dx));
    end
theta_t = k;
end

function rough_t = roughness(C)

    ddx = diff(diff(C(:,1)));
    ddy = diff(diff(C(:,2)));
    rough_t = sum(ddx.^2 + ddy.^2);
    
end

function rgbn = colornorm(rgb)
fprintf('wmspaceseg.m: Vector normalize each color channel\n')
r = double(rgb(:,:,1)); g = double(rgb(:,:,2)); b = double(rgb(:,:,3));
rn = r / sqrt((r(:))' * (r(:)));
gn = g / sqrt((g(:))' * (g(:)));
bn = b / sqrt((b(:))' * (b(:)));
rgbn = cat(3,rn,gn,bn);
end

function maskedImage = applyMask(binaryMask, image)
    maskedImage = uint8(repmat(binaryMask,[1,1,3])).*image;
end

