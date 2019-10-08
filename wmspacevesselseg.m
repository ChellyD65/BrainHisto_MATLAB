function [rgb_out, mask_WM_out, smoothC_WM_out1, mask_Vessel_out, stats] = wmspacevesselseg(fname, varargin) 

%% wnspacevesselseg.m
%  Marcello DiStasio, 2018
%  Segments the non-white matter space on an LFB stained image of a white matter
%  vessel (i.e. perivascular Virchow-Robin space), and computes some statistics on this region.

%% Parse inputs
p = inputParser;
defaultSmoothingKernelSize = 10;
defaultDownsampleFactor = 4;
addRequired(p,'fname');
addOptional(p,'SmoothingKernelSize',defaultSmoothingKernelSize);
addOptional(p,'DownsampleFactor',defaultDownsampleFactor);
parse(p,'fname',varargin{:});
SmoothingKernelSize = p.Results.SmoothingKernelSize;
downsample_factor = p.Results.DownsampleFactor;

%% Load input file
fname
rgb_in = imread(fname); 

%% Downsample
rgb = imresize(rgb_in, (1/downsample_factor));

[mask_LFB, maskedImage_LFB] = createLFBMask(rgb);
mask_LFB_ex = 1-imclearborder(1-imclose(mask_LFB,strel('disk',4)), 4); %A little cleanup
BW = mask_LFB_ex < 1;
BW2 = bwfill(bwpropfilt(BW,'area',1), 'holes'); % Select only the object with largest area, and fill in its holes

mask_WM = BW2;
maskedImage_WM = applyMask(BW2,rgb);

maskedImage_WM_lab = rgb2lab(maskedImage_WM);
L = maskedImage_WM_lab(:,:,1)/100;
L = adapthisteq(L,'NumTiles',[8 8],'ClipLimit',0.005);
maskedImage_WM_lab(:,:,1) = L*100;
maskedImage_WM_contrast = lab2rgb(maskedImage_WM_lab);

I = rgb2gray(maskedImage_WM_contrast);
I_smooth = 1-adapthisteq(imfilter(I,fspecial('gaussian',16,16),'conv')).*imerode(mask_WM, strel('disk',3));

I_diff = (I - I_smooth);
I_diff = I_diff - min(I_diff(:)); I_diff = I_diff/max(I_diff(:));
I_diff = (1-I_diff).*imerode(mask_WM, strel('disk',5));

thh = graythresh(I_diff);
BW = imclearborder(I_diff>2*thh);
BW = imclose(BW,strel('disk',4));

mask_Vessel = bwfill(bwpropfilt(BW,'area',1),'holes');

%% Jaggedness

C=bwboundaries(mask_WM);
C=C{1};
l_orig = contourlen(C);
rough_orig = roughness(C);
theta_orig = contourangle(C);

cc = 1;

for w1 = [3,5,7,9,11,13,15,17,19,21]
    
windowWidth1 = min(w1, 2*floor((l_orig/4)/2)+1 );
smoothX1 = smooth(C(:,2), windowWidth1);
smoothY1 = smooth(C(:,1), windowWidth1);
smoothC_WM1 = [smoothX1, smoothY1];
l_smooth1 = contourlen(smoothC_WM1);
rough_smooth1 = roughness(smoothC_WM1);
theta_smooth1 = contourangle(smoothC_WM1);

%% Return variables (do not rescale back up by downsampling factor)
rgb_out = rgb;
mask_WM_out = mask_WM;
mask_Vessel_out = mask_Vessel;
smoothC_WM_out1 = smoothC_WM1;
J(cc).smoothingWidth1 = w1;
J(cc).JaggedMeasure_len_orig = l_orig;
J(cc).JaggedMeasure_len_smooth1 = l_smooth1;
J(cc).JaggedMeasure_rough_orig = rough_orig;
J(cc).JaggedMeasure_rough_smooth1 = rough_smooth1;
J(cc).JaggedMeasure_theta = theta_orig;
J(cc).JaggedMeasure_theta_smooth1 = theta_smooth1;

cc = cc+1;
end

stats.NonWMarea = bwarea(mask_WM)*downsample_factor;
stats.Vesselarea = bwarea(mask_Vessel)*downsample_factor;
stats.JaggedIndex = J;

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

function maskedImage = applyMask(binaryMask, image)
    maskedImage = uint8(repmat(binaryMask,[1,1,3])).*image;
end

