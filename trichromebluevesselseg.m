function [rgb_out, mask_TriB_out, mask_Lumen_out, stats] = trichromeblueseg(fname, varargin) 

%addpath(fullfile('~','Dropbox','Workspace','MATLAB_TrichromeBlue_Seg','dijkstra'));
addpath(fullfile(pwd(),'dijkstra'))

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
% rgbsmn = histeq(rgbsm,imhist(refimg));

% RGB = mat2gray(rgbsmn);
RGB = mat2gray(rgbsm);
LAB = rgb2lab(RGB);

%%%%
% A good background color (pink) is: b05685 (RGB: 176, 86, 133)

% Thresh = 38;
% Thresh = 42;

a = LAB(:,:,2); an = (a - mean(a(:)))/std(a(:));
b = LAB(:,:,3); bn = (b - mean(b(:)))/std(b(:));

% Thresh = 1250;
% Threshn = 175;
% 
% % AB_blue = -65; %[44, -65];
% % AB_red = 70;
AB_blue = -6; 
AB_red = 1;



% d = ((1-(a - 70)).^2 + 1e3*(b + 65).^2).^0.5;
d = ((1-(an - AB_red)).^2 + 1e3*(bn - AB_blue).^2).^0.5;
dn = (d-mean(d(:)))/std(d(:));

Thresh = -1.75;
% 
% mask_TriB = d < Thresh;
% mask_TriB = dn < Threshn;
mask_TriB = dn < Thresh;
                  

% % % % % LAB = rgb2lab(rgbsm);
% % % % % L = LAB(:,:,1)/100;
% % % % % LAB(:,:,1) = adapthisteq(L,'NumTiles',...
% % % % %                          [8 8],'ClipLimit',0.005)*100;
% % % % % a = LAB(:,:,2);
% % % % % b = LAB(:,:,3);
% % % % %                      
% % % % % nColors = 5;
% % % % % color_markers = zeros([nColors, 2]);
% % % % % lab_color1 = [35.7214 -6.3387]; % blue
% % % % % lab_color2 = [24.4681 -5.8097]; % pink
% % % % % lab_color3 = [31.3658 -5.9495]; % purple
% % % % % lab_color4 = [20.8614 -3.0364]; % pink2
% % % % % lab_color5 = [37.5137  -5.4532];% white
% % % % % 
% % % % % color_markers = [lab_color1; lab_color2; lab_color3; lab_color4; lab_color5];
% % % % % color_labels = 0:nColors-1;
% % % % % 
% % % % % a = double(a);
% % % % % b = double(b);
% % % % % distance = zeros([size(a), nColors]);
% % % % % 
% % % % % for count = 1:nColors
% % % % %   distance(:,:,count) = ( (a - color_markers(count,1)).^2 + ...
% % % % %                       (b - color_markers(count,2)).^2 ).^0.5;
% % % % % end
% % % % % [~, label] = min(distance,[],3);
% % % % % label = color_labels(label);
% % % % % clear distance;
% % % % % 
% % % % % 
% % % % % % rgb_label = repmat(label,[1 1 3]);
% % % % % 
% % % % % mask_TriB = zeros(size(label));
% % % % % mask_TriB(label == 0) = 1;
% % % % % mask_TriB(LAB(:,:,1)>70) = 0;
% % % % % 
% % % % % % segmented_images = zeros([size(rgb), nColors],'uint8');
% % % % % % for count = 1:nColors
% % % % % %   color = rgb;
% % % % % %   color(rgb_label ~= color_labels(count)) = 0;
% % % % % %   segmented_images(:,:,:,count) = color;
% % % % % % end 
% % % % % 
% % % % % % figure; imshow(segmented_images(:,:,:,1)), title('blue objects');
% % % % % % figure; imshow(segmented_images(:,:,:,2)), title('pink objects');
% % % % % % figure; imshow(segmented_images(:,:,:,3)), title('purple objects');
% % % % % % figure; imshow(segmented_images(:,:,:,4)), title('purple2 objects');
% % % % % % figure; imshow(segmented_images(:,:,:,5)), title('white objects');

% [mask_TriB, maskedImage_TriB] = createTrichromeBlueMask2(rgbsm);
% mask_TriB_ex = imopen(mask_TriB,strel('disk',2)); %A little cleanup
mask_TriB_ex = imdilate(mask_TriB,strel('disk',2*(4/downsample_factor)));
mask_TriB_ex = bwareaopen(mask_TriB_ex,600*(4/downsample_factor));


mask_TriB_exF = imfill(mask_TriB_ex,'holes');
boundaries = bwboundaries(mask_TriB_exF);

if ~isempty(boundaries)
    for i = 1:length(boundaries)
        DistFromCenter(i) = min(sqrt(sum((boundaries{i} - [size(rgb,1),size(rgb,2)]/2).^2,2)));
    end
    CenterObject = find(DistFromCenter == min(DistFromCenter));
%     numberOfBoundaries = size(boundaries, 1);
%     for k = 1 : numberOfBoundaries
%         thisBoundary = boundaries{k};
%     end
    % Define object boundaries
    numberOfBoundaries = size(boundaries, 1);
    % Find minimum distance between each pair of boundaries
    Distances = zeros(numberOfBoundaries);
    for b1 = 1 : numberOfBoundaries
        for b2 = 1 : numberOfBoundaries
            if b1 == b2
                % Can't find distance between the region and itself
                continue;
            end
            boundary1 = boundaries{b1};
            boundary2 = boundaries{b2};
            boundary1x = boundary1(:, 2);
            boundary1y = boundary1(:, 1);
            x1=1;
            y1=1;
            x2=1;
            y2=1;
            overallMinDistance = inf; % Initialize.
            % For every point in boundary 2, find the distance to every point in boundary 1.
            for k = 1 : size(boundary2, 1)
                % Pick the next point on boundary 2.
                boundary2x = boundary2(k, 2);
                boundary2y = boundary2(k, 1);
                % For this point, compute distances from it to all points in boundary 1.
                allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
                % Find closest point, min distance.
                [minDistance(k), indexOfMin] = min(allDistances);
                if minDistance(k) < overallMinDistance
                    x1 = boundary1x(indexOfMin);
                    y1 = boundary1y(indexOfMin);
                    x2 = boundary2x;
                    y2 = boundary2y;
                    overallMinDistance = minDistance(k);
                end
            end
            % Find the overall min distance
            minDistance = min(minDistance);
            Distances(b1,b2) = minDistance;
        end
    end

    maxJump = .02 * size(rgb,1);
    m = size(mask_TriB_exF,1);
    n = size(mask_TriB_exF,2);
    mask_TriB_Selected = zeros(m,n);
    for i = 1:numberOfBoundaries
        [cost, route] = dijkstra(Distances,CenterObject,i);
        if length(route) > 1
            if (max(Distances(sub2ind(size(Distances),route(1:end-1), route(2:end)))) < maxJump)
                mask_TriB_Selected = mask_TriB_Selected + poly2mask(boundaries{i}(:,2), boundaries{i}(:,1),m,n).*mask_TriB_ex;
            end
        else
            mask_TriB_Selected = mask_TriB_Selected + poly2mask(boundaries{i}(:,2), boundaries{i}(:,1),m,n).*mask_TriB_ex;
        end
    end









%     mask_TriB_exC =  imerode(imfill(mask_TriB_Selected,'holes'),strel('disk',12));
    dsize = 14*(4/downsample_factor);
    mask_TriB_exC = imerode(imfill(imdilate(mask_TriB_Selected,strel('disk',dsize)),'holes'),strel('disk',dsize+(10/downsample_factor)));
    mask_TriB_coarse = mask_TriB_Selected;
    
else
%     mask_TriB_exC =  imerode(imfill(mask_TriB_ex,'holes'),strel('disk',12));
    dsize = 14*(4/downsample_factor);
    mask_TriB_exC = imerode(imfill(imdilate(mask_TriB_ex,strel('disk',dsize)),'holes'),strel('disk',dsize+(10/downsample_factor)));
    mask_TriB_coarse = mask_TriB_ex;
    
end




RGB = mat2gray(rgb);
LAB = rgb2lab(RGB);
a = mask_TriB_coarse.*LAB(:,:,2); an = (a - mean(a(:)))/std(a(:)); AB_red = 7;
b = mask_TriB_coarse.*LAB(:,:,3); bn = (b - mean(b(:)))/std(b(:)); AB_blue = -10; 
d = ((1-(an - AB_red)).^2 + 1e3*(bn - AB_blue).^2).^0.5;
dn_fine = (d-mean(d(:)))/std(d(:));

mask_TriB_fine = dn_fine < Thresh;
mask_TriB_out = mask_TriB_fine;

BW = mask_TriB_exC.*(dn>Thresh) > 0;
mask_Lumen = bwfill(bwpropfilt(BW,'area',1), 'holes'); % Select only the object with largest area, and fill in its holes

rgb_out = rgb;
mask_Lumen_out = mask_Lumen;

stats.CollagenArea = bwarea(mask_TriB_out)*downsample_factor;
stats.VesselLumenArea = bwarea(mask_Lumen_out)*downsample_factor;
stats.VesselLumenPerim = sum(sum(bwperim(mask_Lumen_out)))*downsample_factor;


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
% rn = r - min(r(:)); rn = rn ./ max(rn(:)); 
% gn = g - min(g(:)); gn = gn ./ max(gn(:)); 
% bn = b - min(b(:)); bn = bn ./ max(bn(:)); 
rn = r / sqrt((r(:))' * (r(:)));
gn = g / sqrt((g(:))' * (g(:)));
bn = b / sqrt((b(:))' * (b(:)));
rgbn = cat(3,rn,gn,bn);
end

% 
% function rgbnt = purplethresh(rgbn)
%     hsvv = rgb2hsv(rgbn);
%     hsvv_purple = hsvv(:,:,1) > 0.74 & hsvv(:,:,1) < 0.82;
%     rgbnt=bsxfun(@times,rgbn,hsvv_purple);
% end

function maskedImage = applyMask(binaryMask, image)
    maskedImage = uint8(repmat(binaryMask,[1,1,3])).*image;
end

