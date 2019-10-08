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

DownsampleFactor = 4
ii = 1;
for f = ilist
    fname = imfiles(f).name;
    
    % Check to see if an edited version is present -- this allows the user
    % to manually create a 'mask' restricting analysis to a certain part of
    % the image.
    
    ff = regexp(fname,'\.','split');
    ffed = fullfile(fbase, 'edited', strcat(ff{1},'_ed.jpg'));
    if (exist(ffed,'file'))
        fprintf('\n\nFound edited version: %s\n', ffed);
        fname_seg = ffed;
    else
        fname_seg = fullfile(fbase,imfiles(f).name);
    end
    fprintf('\nFile (%i of %i): %s\n',ii,length(ilist),fname_seg);
    [rgb, mask_WM, smoothC_WM1, mask_Vessel, stats] = wmspacevesselseg(fname_seg, 'DownsampleFactor', DownsampleFactor);
    if (exist(ffed,'file'))
        % Display original image even if edited is present
        rgb = imresize(imread(fullfile(fbase,imfiles(f).name)),1/DownsampleFactor);
    end
    
    h = figure('visible','off'); 
    imshow(rgb, 'InitialMag', 'fit'); hold on;
    
    C = bwboundaries(mask_WM); C=C{1};
    plot(C(:,2), C(:,1),'-y','LineWidth',1);
    saveas(h, fullfile(outdir,strcat(fname,'_out.jpg')))
    
    
    %% Save some metrics in the imfiles struct
    imfiles(f).NonWMarea = stats.NonWMarea;
    imfiles(f).Vesselarea = stats.Vesselarea;
    imfiles(f).JaggedIndex = stats.JaggedIndex;

    ii = ii+1;
end

%% Write out a CSV file with the areas
outfilenamebase = strcat('wmspaceseg_Run_-_',datestr(now, 'yyyy-mm-dd_HH-MM-SS'))

outfile = fullfile(outdir, strcat(outfilenamebase,'.csv'))
fileID = fopen(outfile,'w');

cnn = fieldnames(imfiles(1).JaggedIndex);
A = cnn(2:7)';
B = string(cat(1,imfiles(1).JaggedIndex.smoothingWidth1));
[ii,jj]=ndgrid(1:numel(A),1:numel(B));
out=arrayfun(@(x,y) [strcat(A(y),"_",B(x))],jj(:),ii(:),'un',0);

% fprintf(fileID, 'Case,Block,VesselNum,NonWMarea,Vesselarea,JaggedIndex,JaggedIndex2,JaggedIndex3,JaggedIndex4\r\n');
fprintf(fileID, strjoin(['Case,Block,VesselNum,NonWMarea,Vesselarea,',strjoin(strcat([out{:}]),","),'\r\n']));


for i = ilist
    fss = imfiles(i).name;
    % [p,fss,e] = fileparts(filepathh);
    expression = '_|\.';
    splitStr = regexp(fss,expression,'split');
    imfiles(i).Case = splitStr{1};
    imfiles(i).Block = splitStr{2};
    imfiles(i).VesselNum = splitStr{4};
    
    fprintf(fileID, '%s,%s,%s,%g,%g', ...
    imfiles(i).Case, ...
    imfiles(i).Block, ...
    imfiles(i).VesselNum, ...
    imfiles(i).NonWMarea, ...
    imfiles(i).Vesselarea);

    J = imfiles(i).JaggedIndex;
%     for j = [1:size(J,2)]
%         if j < size(J,2)
%         fprintf(fileID,',%g,%g,%g', ...
%             J(j).JaggedMeasure_len, ...
%             J(j).JaggedMeasure_theta, ...
%             J(j).JaggedMeasure_theta_smooth);
%         else
%         fprintf(fileID,',%g,%g,%g\r\n', ...
%             J(j).JaggedMeasure_len, ...
%             J(j).JaggedMeasure_theta, ...
%             J(j).JaggedMeasure_theta_smooth);
%         end
%     end
    for j = [1:size(J,2)]
        fprintf(fileID,',%g,%g,%g,%g,%g,%g', ...
            J(j).JaggedMeasure_len_orig, ...
            J(j).JaggedMeasure_len_smooth1, ...
            J(j).JaggedMeasure_rough_orig, ...
            J(j).JaggedMeasure_rough_smooth1, ...
            J(j).JaggedMeasure_theta, ...
            J(j).JaggedMeasure_theta_smooth1);
    end    
    fprintf(fileID,'\r\n');



%     % Write a line to the output table
%     fprintf(fileID, '%s,%s,%s,%g,%g,%g,%g,%g,%g\r\n', ...
%         imfiles(i).Case, ...
%         imfiles(i).Block, ...
%         imfiles(i).VesselNum, ...
%         imfiles(i).NonWMarea, ...
%         imfiles(i).Vesselarea, ...
%         imfiles(i).JaggedIndex, ...
%         imfiles(i).JaggedIndex2, ...
%         imfiles(i).JaggedIndex3, ...
%         imfiles(i).JaggedIndex4);
end

fclose(fileID);
fprintf('\n\nSaved File: %s\n', outfile);


toc;

