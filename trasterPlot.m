clear all

%% worm trace raster plot %%%%
%%%% outout format %%%%
    %%%% traces in raster colour %%%%

%% open an image first %%%%

%% user input %%%%
start = 480;
finish = 510;
frameRate = 0.2;
choreP = 0.027;
%%%% user input end %%%%

dirDat = dir('*.dat');
hold on;
%% loop through plate %%%%
for i = 1:numel(dirDat)
    wormTrace = [];
    wormDat = dlmread(dirDat(i).name);
    if wormDat(1,1) > finish || wormDat(end,1) < start
        continue;
    end
    for j = 1:size(wormDat,1)
        if wormDat(j,1) > start && wormDat(j,1) <= finish
            wormTrace = [wormTrace; wormDat(j,1) wormDat(j,4).*wormDat(j,6) wormDat(j,9)/choreP wormDat(j,10)/choreP]
        end
    end
    wormFor = [(wormTrace(:,2)>0).*wormTrace(:,4) (wormTrace(:,2)>0).*wormTrace(:,3)];
    wormPause = [(wormTrace(:,2)==0).*wormTrace(:,4) (wormTrace(:,2)==0).*wormTrace(:,3)];
    wormRev = [(wormTrace(:,2)<0).*wormTrace(:,4) (wormTrace(:,2)<0).*wormTrace(:,3)];
    wormFor( ~any(wormFor,2), : ) = [];
    wormPause( ~any(wormPause,2), : ) = [];
    wormRev( ~any(wormRev,2), : ) = [];
    plot(wormFor(:,1), wormFor(:,2),'r', wormPause(:,1), wormPause(:,2),'y', wormRev(:,1), wormRev(:,2),'b', 'LineWidth',1.5);
end
