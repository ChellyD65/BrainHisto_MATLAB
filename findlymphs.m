function [labels, stats, stats1] = findlymphs(fname, varargin) 

%% Parse inputs
p = inputParser;
defaultRadiusList = [14,22];
defaultNucSaturationThresh = 0.65;
defaultCircThresh = 0.8;
addRequired(p,'fname');
addOptional(p,'RadiusLimits',defaultRadiusList);
addOptional(p,'NucSaturationThresh',defaultNucSaturationThresh);
addOptional(p,'CircThresh',defaultCircThresh);
parse(p,'fname',varargin{:});
RadiusList = p.Results.RadiusLimits;
NucSaturationThresh = p.Results.NucSaturationThresh;
CircThresh = p.Results.CircThresh;

%% Load input file
rgb = imread(fname); 

%% Color Separation and thresholding;
rgbn = colornorm(rgb); % Normalize each color channel
hh = rgb2hsv(rgbn);
% Get Intensity (actually 'value')
% I = rgb2gray(rgbn);
I = hh(:,:,3);

% Smooth
fprintf('findlymphs.m: Smooth\n')
h = fspecial('gaussian',100,5);
I=imfilter(I,h,'conv');

[WhiteMatterMask] = createLFBMask(rgb);
WhiteMatterMask = bwareaopen(WhiteMatterMask,2e2);
for i = 1:8
    WhiteMatterMask = imdilate(WhiteMatterMask,strel('disk',8,8));
end
WhiteMatterMask = 1-bwareaopen(1-WhiteMatterMask,2e4);

% Threshold on intensity
fprintf('findlymphs.m: Threshold on HSV(:,:,3)\n')
level = graythresh(I); %Otsu threshold intensity
I_t = I; I_t(I_t > level) = 0;

fprintf('findlymphs.m: Apply (negative) white matter mask\n');
I_t = I_t.*(1-WhiteMatterMask);

%% Morphology
fprintf('findlymphs.m: Opening operation (remove noise specks)\n')
I_th = imhmin(I_t,max(I_t(:))*0.1);
I_tho = I_th.*bwareaopen(I_th,32); % Gets rid of noise (dots w/ area < 16px)
I_tho(I_t==0) = -Inf;

LL = watershed(I_tho); 

fprintf('findlymphs.m: Erode\n')
se = strel('disk',4,4); LL = bwlabel(imerode(LL,se)>0); %erode
stats1 = regionprops(LL,'Area','Centroid','Eccentricity', 'Perimeter'); 

fprintf('findlymphs.m: Delete any blobs with area above 1.5*2*pi*max_radius (%g) or below 0.6*2*pi*min_radius (%g)\n', ...
    pi*max(1.5*RadiusList(:)).^2, pi*min(0.6*RadiusList(:)).^2)
areas = cat(1,stats1(2:end).Area); % Ignore the background (stats(1))
for i = find(areas > pi*max(1.5*RadiusList(:)).^2 | areas < pi*min(0.6*RadiusList(:)).^2)'
   LL(LL==i+1) = 0; 
end

% Area / perimeter will approach a maximum of 0.5*r for circle
fprintf('findlymphs.m: Delete any blobs with circularity below %g\n', CircThresh)
cmm = 4*pi*cat(1,stats1(:).Area)./cat(1,stats1(:).Perimeter).^2;
cmm(isinf(cmm)) = 0;
for i = find(cmm < CircThresh)'
   LL(LL==i) = 0; 
end

%% Segmentation of circles
LL_circ = uint16(zeros(size(LL)));
for mcr = 1:size(RadiusList,1)
mcrl = RadiusList(mcr,1);
mcrh = RadiusList(mcr,2);
    fprintf('findlymphs.m: Find circles (r in range [%g,%g])\n', mcrl,mcrh)
    [centers, radii] = imfindcircles(LL>1,[mcrl, mcrh],'ObjectPolarity','bright', ...
        'Sensitivity',0.99, 'EdgeThreshold',0.95);
    fprintf('findlymphs.m: Found %g circles (r in range [%g,%g])\n', size(centers,1), mcrl,mcrh)  
    % Select regions that share a center with the circles.
    if (size(centers,1)>0)
        radius_median = median(sqrt(cat(1,stats1(:).Area)/pi));
        for i = 1:size(stats1)
            c = repmat(stats1(i).Centroid,size(centers,1),1);
            mindist = min(sqrt(sum((centers - c).^2,2)));
            stats1(i).mindist = mindist;
            if mindist < radius_median
                LL_circ(LL==i) = i;
            end
        end
    end
end

stats2 = regionprops(LL_circ,hh(:,:,2),'Area','Centroid','Eccentricity','MeanIntensity','Perimeter');

fprintf('findlymphs.m: Get mean value in circles, select those > %g\n',NucSaturationThresh)
I_m = cat(1,stats2.MeanIntensity);
I_metric = I_m;
%I_metric = (I_m - min(I_m)); I_metric = I_metric./max(I_metric); %Normalize metric to [0,1]

Hits = uint16(find(I_metric>NucSaturationThresh));
statsOut = stats2(Hits);
LLOut = LL_circ.*uint16(ismember(LL_circ,Hits));
fprintf('findlymphs.m: Found %g circles (saturation > %g)\n', size(Hits,1), NucSaturationThresh)  

%% Return variables
stats = statsOut;
labels = LLOut;

end

function rgbn = colornorm(rgb)
fprintf('findlymphs.m: Normalize each color channel\n')
r = double(rgb(:,:,1)); g = double(rgb(:,:,2)); b = double(rgb(:,:,3));
rn = r - min(r(:)); rn = rn ./ max(rn(:)); 
gn = g - min(g(:)); gn = gn ./ max(gn(:)); 
bn = b - min(b(:)); bn = bn ./ max(bn(:)); 
rgbn = cat(3,rn,gn,bn);
end

function rgbnt = purplethresh(rgbn)
    hsvv = rgb2hsv(rgbn);
    hsvv_purple = hsvv(:,:,1) > 0.74 & hsvv(:,:,1) < 0.82;
    rgbnt=bsxfun(@times,rgbn,hsvv_purple);
end
