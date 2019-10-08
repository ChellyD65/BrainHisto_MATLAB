%% Run_findlymphs_ondir.m
% Marcello DiStasio, 2019

clear; close all;
addpath(pwd);
tic;

fbase =  fullfile('.','images')
imfiles = [dir(fullfile(fbase,'*.tif')); 
           dir(fullfile(fbase,'*.jpg')); 
           dir(fullfile(fbase,'*.jpeg'))]

outdir = fullfile(fbase,'out')
if (exist(outdir,'dir') ~= 7)
    mkdir(outdir)
end

%% Pick 50 random files
% ilist = ceil(size(imfiles,1)*rand(1,50));


%% Pick one file
% imfiles = dir(fullfile(fbase,'MSSM14-2_014.jpg'));
% imfiles = dir(fullfile(fbase,'Slide_1_009.jpg'));
% ilist = 1;

%% Pick all files
ilist = 1:size(imfiles,1);

%% Loop through files in ilist

ii = 1;
for f = ilist
    fname = fullfile(fbase,imfiles(f).name);
    fprintf('\nFile (%i of %i): %s\n',ii,length(ilist),fname)
    % Tweak the morphologic parameters here
    [l,s,s1] = findlymphs(fname, 'NucSaturationThresh', 0.25, 'RadiusLimits', [9, 18], 'CircThresh', 0.5);
    imfiles(f).s = s;
    imfiles(f).l = l;
    imfiles(f).s1 = s1;
    
    hh = figure; set(gca,'position',[0 0 1 1],'units','normalized'); set(hh, 'Visible', 'off');
    imshow(imread(fname)); hold on;
    h = viscircles(cat(1,imfiles(f).s.Centroid),sqrt(cat(1,imfiles(f).s.Area)/pi),'LineWidth',0.5);

    saveas(hh, strcat(fullfile(outdir,imfiles(f).name),'_-_LymphocyteSearch_Run_-_',datestr(now, 'yyyy-mm-dd_HH-MM-SS'),'.jpg'),'jpeg')
    close(hh);
    ii = ii+1;
end

%% Save summary (i.e. updated imfiles)
save(fullfile(outdir,strcat('LymphocyteSearch_Run_-_',datestr(now, 'yyyy-mm-dd_HH-MM-SS'),'.mat')),'imfiles');

toc;

