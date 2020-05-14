function fig1 = starFish(wormid,X,Y,Bias,Tap,varargin)


%% varargin default
visibleopt = 1;
wormidshow = 0;
tapcirclesize = 3;
pathlinewidth = 0.2;
alignAtTap = 0;
respondantsOnly = 0;

%% VARARGIN PROCESSOR
A = varargin;
if numel(A) > 1
    if isinteger(numel(A)/2) == 1
        error('variable entries incorrect')
    end
    callName = A(1:2:numel(A));
    % check if call names are strings
    for x = 1:numel(callName)
       if ischar(callName{x}) == 0
           error('some var names not strings');
       end
    end
    setValue = A(2:2:numel(A)); % get input values
    % assign input to var by input type
    for x = 1:numel(callName)
        % if is char
        if eval(sprintf('ischar(%s)',callName{x})) == true
            eval(sprintf('%s = ''%s'';',callName{x},setValue{x}));
        % if is numeric
        elseif eval(sprintf('isnumeric(%s)',callName{x})) == true || ...
            eval(sprintf('islogical(%s)',callName{x})) == true 
            eval(sprintf('%s = %d;',callName{x},cell2mat(setValue(x))));
        % if is table or cell
        elseif eval(sprintf('iscell(%s)',callName{x})) == true ||...
               eval(sprintf('istable(%s)',callName{x})) == true
            a = setValue{x};
            eval(sprintf('%s = a;',callName{x}));
        % other types not coded here
        else
            error('need coding for %s',callName{x});
        end
    end
end


%% get wormid
wormidu = unique(wormid);
wormN = numel(wormidu);
% loop thru each worm
Biascolorlegend = [-1,[0 0 1];0,[.5 .5 .5];1,[1 0 0]];
close;
if visibleopt
    fig1 = figure('Visible','on','Color',[1 1 1]);
else
    fig1 = figure('Visible','off','Color',[1 1 1]);
end
axes1 = axes('Parent',fig1,'ZColor','w','YColor','w','XColor','w');
hold(axes1,'all');

for wrmi = 1:wormN
   iworm = find(wormid == wormidu(wrmi));
   % respondants only
   if respondantsOnly
       itap = find(Tap(iworm));
       bias = Bias(iworm);
       biaschange = unique(bias(itap:itap+2));
       if numel(biaschange) == 1
          if biaschange == bias(itap-1); continue; end
       end
   end
   if alignAtTap % start with tap
      iworm(1:find(Tap(iworm))-1) = [];
   end
   xp = X(iworm);
   yp = Y(iworm);
   bias = Bias(iworm)';
   tap = Tap(iworm);

   shiftpt = [1 find(diff(bias) ~= 0)+1 numel(bias)+1];

   % get position, set first value to zero
   x = xp - xp(1);
   y = yp - yp(1);
   % draw tap circle
   itap = tap == 1;
   plot(x(itap),y(itap),'Color','k','MarkerFaceColor','none','Marker','o','MarkerEdgeColor',[0 0 0],'MarkerSize',tapcirclesize,'LineStyle','none','LineWidth',.5);

   % draw path
   for sfi = 1:numel(shiftpt)-1
       t1 = shiftpt(sfi);
       t2 = shiftpt(sfi+1)-1;
       b = unique(bias(t1:t2));
       if isnan(b) == 1 && numel(b)== 1
           c = [ 0 1 0];
       else
            c = Biascolorlegend(ismember(Biascolorlegend(:,1),b),2:4);
       end
       if sfi == 1
          xx = x(t1:t2);
          yy = y(t1:t2);
       else
          xx = x(t1-1:t2);
          yy = y(t1-1:t2);
       end
       plot(xx,yy,'Color',c,'MarkerFaceColor','none','Marker','.','MarkerEdgeColor',c,'MarkerSize',2,'LineWidth',pathlinewidth);
   end
   % plot wormid 
   if wormidshow
      text(x(end),y(end),num2str(wormidu(wrmi)),'FontSize',2)
   end
   xlim([-1.1 1.1]);
   ylim([-1.1 1.1]);
end
