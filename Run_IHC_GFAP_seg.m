% Run_IHC_GFAP_seg.m
% Marcello DiStasio, June 2019

addpath(fullfile('.'));
clear; close all;
tic;

fbase = fullfile('.','images')

imfiles = dir(fullfile(fbase, '*.jpg')) 
ilist = 1:size(imfiles,1); % Pick all files

tic;
outdir = fullfile(fbase,'out')
if (exist(outdir,'dir') ~= 7)
    mkdir(outdir)
end

DownsampleFactor = 1

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
    
    [rgb, mask_GFAP, stats] = IHC_GFAP_seg(fname_seg, 'DownsampleFactor', DownsampleFactor, 'RefImage', refimg_rgbsm);
    if (exist(ffed,'file'))
        % Display original image even if edited is present
        rgb = imresize(imread(fullfile(fbase,imfiles(f).name)),1/DownsampleFactor);
    end
    
    h = figure('visible','off'); 
    imshow(rgb, 'InitialMag', 'fit'); hold on;
    
    C = bwboundaries(mask_GFAP); 
    for cc = 1:size(C,1)
        B=C{cc};
        plot(B(:,2), B(:,1),'-g','LineWidth',1);
    end
    saveas(h, fullfile(outdir,strcat(fname,'_out.jpg')))
    
    close all hidden;

    %% Save some metrics in the imfiles struct
    imfiles(f).GFAPArea = stats.GFAPArea;
%     imfiles(f).VesselLumenArea= stats.VesselLumenArea;
%     imfiles(f).VesselLumenPerim = stats.VesselLumenPerim;

    ii = ii+1;
end
        
%% Write out a CSV file with the areas
outfilenamebase = strcat('IHC_GFAP_seg_Run_-_',datestr(now, 'yyyy-mm-dd_HH-MM-SS'))

outfile = fullfile(outdir, strcat(outfilenamebase,'.csv'))
fileID = fopen(outfile,'w');

fprintf(fileID, 'Case,Block,PhotoNum,GFAParea\r\n');

for i = ilist
    fss = imfiles(i).name;
    % [p,fss,e] = fileparts(filepathh);
    expression = '_|\.';
    splitStr = regexp(fss,expression,'split');
    imfiles(i).Case = splitStr{1};
    imfiles(i).Block = splitStr{2};
    imfiles(i).PhotoNum = str2num(splitStr{4});
    
    fprintf(fileID, '%s,%s,%g,%g,\n', ...
    imfiles(i).Case, ...
    imfiles(i).Block, ...
    imfiles(i).PhotoNum, ...
    imfiles(i).GFAPArea);

end

fclose(fileID);
fprintf('\n\nSaved File: %s\n', outfile);

toc;

