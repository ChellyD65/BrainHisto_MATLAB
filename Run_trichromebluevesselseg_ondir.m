addpath(fullfile(pwd()))
clear; close all;
tic;

% fbase = fullfile('~','Dropbox','AndersonLab','Lymphocytes_Study','Data','Trichrome_Photos_VESSELS','SecondSet','good')
fbase = fullfile('~','AndersonLab','Lymphocytes_Study','Data','Trichrome_Photos_VESSELS','NewControlCases')

imfiles = dir(fullfile(fbase, '*.jpg')) 
ilist = 1:size(imfiles,1); % Pick all files
% ilist = 1:4;

tic;
outdir = fullfile(fbase,'out')
if (exist(outdir,'dir') ~= 7)
    mkdir(outdir)
end

DownsampleFactor = 1

% refimg = fullfile(fbase,'B6184_17_VN_020.jpg')
% refimg_rgb = imread(refimg);
% refimg_rgb = imresize(refimg_rgb, (1/DownsampleFactor));
% refimg_rgbsm = imfilter(refimg_rgb,fspecial('gaussian',8,8),'conv');
refimg_rgbsm = 0;

ii = 1;
for f = ilist
    fname = imfiles(f).name;
    fprintf('\nFile (%i of %i): %s\n',ii,length(ilist),fname);
    % Check to see if an edited version is present
    ff = regexp(fname,'\.','split');
    ffed = fullfile(fbase, 'edited', strcat(ff{1},'_ed.jpg'));
    if (exist(ffed,'file'))
        fprintf('\nFound edited version: %s\n', ffed);
        fname_seg = ffed;
    else
        fname_seg = fullfile(fbase,imfiles(f).name);
    end
    
%     [rgb, mask_WM, smoothC_WM1, mask_Vessel, stats] = wmspacevesselseg(fname_seg, 'DownsampleFactor', DownsampleFactor);
%     [rgb, mask_TriB, stats] = trichromebluevesselseg(fname_seg, 'DownsampleFactor', DownsampleFactor);
    [rgb, mask_TriB, mask_Lumen, stats] = trichromebluevesselseg(fname_seg, 'DownsampleFactor', DownsampleFactor, 'RefImage', refimg_rgbsm);
    if (exist(ffed,'file'))
        % Display original image even if edited is present
        rgb = imresize(imread(fullfile(fbase,imfiles(f).name)),1/DownsampleFactor);
    end
    
    h = figure('visible','off'); 
    imshow(rgb, 'InitialMag', 'fit'); hold on;
    
    C = bwboundaries(mask_TriB); 
    for cc = 1:size(C,1)
        B=C{cc};
        plot(B(:,2), B(:,1),'-g','LineWidth',1);
    end
    C = bwboundaries(mask_Lumen); 
    for cc = 1:size(C,1)
        B=C{cc};
        plot(B(:,2), B(:,1),'-y','LineWidth',1);
    end
    
    
    saveas(h, fullfile(outdir,strcat(fname,'_out.jpg')))
    
    close all hidden;

    %% Save some metrics in the imfiles struct
    imfiles(f).CollagenArea = stats.CollagenArea;
    imfiles(f).VesselLumenArea= stats.VesselLumenArea;
    imfiles(f).VesselLumenPerim = stats.VesselLumenPerim;

    ii = ii+1;
end

        
%% Write out a CSV file with the areas
outfilenamebase = strcat('trichromeblueseg_Run_-_',datestr(now, 'yyyy-mm-dd_HH-MM-SS'))

outfile = fullfile(outdir, strcat(outfilenamebase,'.csv'))
fileID = fopen(outfile,'w');


fprintf(fileID, 'Case,Block,SulcusPhotoNum,CollagenArea,LumenArea,LumenPerim\r\n');


for i = ilist
    fss = imfiles(i).name;
    % [p,fss,e] = fileparts(filepathh);
    expression = '_|\.';
    splitStr = regexp(fss,expression,'split');
    imfiles(i).Case = splitStr{1};
    imfiles(i).Block = splitStr{2};
    imfiles(i).SulcusPhotoNum = str2num(splitStr{4});
    
    fprintf(fileID, '%s,%s,%g,%g,%g,%g\n', ...
    imfiles(i).Case, ...
    imfiles(i).Block, ...
    imfiles(i).SulcusPhotoNum, ...
    imfiles(i).CollagenArea,...
    imfiles(i).VesselLumenArea, ...
    imfiles(i).VesselLumenPerim);

end

fclose(fileID);
fprintf('\n\nSaved File: %s\n', outfile);


toc;

