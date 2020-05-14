p = '/Users/connylin/Dropbox/RL/PhD/Chapters/2-STH N2/Data/Recovery/Groups/Recovery 10sISI/Dance_Glee_Showmance/Dance_Glee_Showmance.mat';
load(p);
pMWT = MWTSet.PATHS.pMWT;

% run Dance_Glee_Showmance
addpath('/Users/connylin/Dropbox/RL/Code/Modules/Dance/Dance_Glee');
MWTSet = Dance_Glee_Acafellas(pMWT,'analysisNLimit',3);

%% import swanlake.dat (speed for Reversal)
% import_swanlake

% oshanespark = '-O shanespark -o nNss*b12M'; % standard for Beethoven
% oswanlake = '-O swanlake -o tnNemMawlkcspbd1e#m#M#a#w#l#k#c#s#p#b#d#e-m-M-a-w-l-k-c-s-p-b-d-e*m*M*a*w*l*kvc*s*p*b*d*';
[~,p] = dircontent(cd,'*swanlake2.dat');
a = dlmread(char(p));
%% tap is at col 16, speed at 12, 25, 37, good N at 3, time at 1
% reversal speed is defined as negative speed during a duration of
% consecutive reversals started within 1 second after the tap occurs

t1 = 319.9710;
t2 = t1+1.5;
i = a(:,1) >= t1 & a(:,1) <=t2;
a(i,3)