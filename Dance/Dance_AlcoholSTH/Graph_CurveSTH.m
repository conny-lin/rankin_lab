function Graph_CurveSTH(X,Y,E,gn,msr,pSave,varargin)

%% defaults
graphpack = 'brnfsc';
gp = graphsetpack(graphpack);
graphname = 'HabCurve';
cStrain = 'N2';
%% varargin processor 
vararginProcessor;


%% make sure gn matches ------------------------------
gnT = regexprep(gn,' ','_');
% s = parse_groupname(gnT,{'strain'});
% strains = unique(s);
% mutant = strains(~ismember(strains,cStrain));
% r = parse_groupname(gnT,{'rx'});
% rx =  unique(r);
% 
% if numel(mutant)==1
%     mutant = char(mutant); 
% elseif isempty(mutant)
%     fprintf('no mutants passed valid N\n');
% else
%     error('graph can not accomodate multiple mutants');
% end
% ----------------------------------------------------

% create standard matches -----------------------------
% if ~isempty(mutant)
%     expectedGn = {[cStrain rx{1}] [cStrain ' ' rx{2}] [mutant rx{1}] [mutant ' ' rx{2}]};
%     gp.DisplayName = expectedGn;
% else
%     expectedGn = {[cStrain rx{1}] [cStrain ' ' rx{2}] '' ''};
%     gp.DisplayName = expectedGn;   
% end
gp.DisplayName = gnT;
% ----------------------------------------------------


%% plot -------------------------------------------
fig = figure('Visible','off','PaperSize',[4 2.5]); 
ax1 = axes('Parent',fig); hold(ax1,'on');
e1 = errorbar(X,Y,E,'Marker','o','MarkerSize',3);
% -----------------------------------------------

% get settings -----------------------------------
for gi = 1:numel(gn)
    
%     gpn = fieldnames(gp);
    % find expected group name
%     k = find(ismember(gp.DisplayName,regexprep(gn{gi},'_',' ')));
%     
%     if numel(k) ~= 1
%         warning('multiple group name per graph setting assignment'); 
%         k = k(1); 
%     end
    
%     for gpi = 1:numel(gpn)
%         e1(gi).(gpn{gpi}) = gp.(gpn{gpi}){k};
%     end
    e1(gi).DisplayName = gp.DisplayName{gi};
    
end
% -----------------------------------------

% legend ------------------------
lg = legend('show');
lg.Location = 'northeastoutside';
lg.Box = 'off';
lg.Color = 'none';
lg.EdgeColor = 'none';
lg.FontSize = 8;
% -------------------------------

% axis -----------------------------------------------
xlabel('stimuli')
xlim([0 max(max(X))+0.5]);
if strcmp(msr,'RevFreq')
    yname = 'RevProb';
%     ylim([0 1]);
else
    yname = msr;
%     ylim([0 max(max(Y))*1.1]);
end
ylabel(yname);
axs = get(ax1);
ax1.TickLength = [0 0];
ax1.FontSize = 10;
%% ----------------------------------------------------

% save -------------------------------------------
savename = sprintf('%s/%s %s',pSave,graphname,msr);
printfig(savename,pSave,'w',4,'h',2.5,'closefig',1);
% -------------------------------------------

