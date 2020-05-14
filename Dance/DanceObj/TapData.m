classdef TapData
    properties
        mwtpath
        mwtid
        groupname
        tap
        Response
        Measures
    end
    properties (Dependent)
        MWTDB
        gnu
        gnu_genotype
    end
    
    methods
        % get class properties
        function MWTDB = get.MWTDB(obj)
            p = obj.mwtpath;
            MWTDB = parseMWTinfo(unique(p));
        end
        function mwtid = get.mwtid(obj)
            p = obj.mwtpath;
            M = obj.MWTDB;
            [~,j] = ismember(p,M.mwtpath);
            mwtid = j;
        end
        function gnu = get.gnu(obj)
            gnu = unique(obj.groupname);
        end
        function gnu_genotype = get.gnu_genotype(obj)
            g = obj.gnu;
            a = regexpcellout(g,'^[A-Z]{1,}\d{1,}','match');
            strainNames = DanceM_load_strainInfo(a);
            strain = strainNames.strain;
            genotype = strainNames.genotype;
            for si = 1:numel(strain)
               g = regexprep(g,strain{si}, genotype{si});
            end
            gnu_genotype = g;
            gnu_genotype = regexprep(gnu_genotype,'N2','wildtype');
        end
        function rmANOVAData = rmANOVA(obj)
            % get inputs
            rm = obj.tap;
            dv = obj.Response;
            id = obj.mwtpath;
            msrlist = obj.Measures;
            % get tap (repeated measures variable)
            rmu = unique(rm);
            idu = unique(id);
            nid = numel(idu);
            % get output
            A = cell(nid,1);
            G = cell(nid,1);
            mwtdbu = parseMWTinfo(idu);
            [~,mwtnames] = cellfun(@fileparts,id,'UniformOutput',0);
            % reorganize data
            [~,rowind] = ismember(mwtnames,mwtdbu.mwtname);
            [~,colind] = ismember(rm,rmu);
            A = nan(numel(idu),numel(rmu));
            arraySize = size(A);
            sind = sub2ind(arraySize,rowind,colind);
            if numel(unique(sind)) > numel(A) % if duplicated data found
                error('found duplicated data');
            end
            S = struct;
            for msri = 1:numel(msrlist)
                msr = msrlist{msri};
                B = A;
                d = dv(:,msri);
                B(sind) = d;
                B(isnan(B)) = 0; % !!! fill in nan data with zero to avoid no data issue in rmanova
                S.(msr) = B;
            end
            % find dose
            a = mwtdbu.groupname;
            b = regexpcellout(a,'(?<=(_))\d{1,}mM','match');
            b(cellfun(@isempty,b)) = {'0mM'};
            dose = b;
            % find strain
            a = mwtdbu.groupname;
            strain = regexpcellout(a,'^[A-Z]{1,}\d{1,}','match');
            % construct factors
            a = [mwtdbu.mwtname, mwtdbu.groupname, strain, dose];
            F = cell2table(a,'VariableNames',{'mwtname','groupname','strain','dose'});
            % =========================================================================

            % rmanova settings ++++++++++++++++++++++++
            alpha = 0.05;
            pvlimit = 0.001;
            compName = 'groupname';
            rmName = 'tap';
            % determine factorName
            strainN = numel(unique(F.strain));
            doseN = numel(unique(F.dose));
            if strainN > 1 && doseN > 1
                factorName = 'strain*dose';
            elseif strainN == 1
                factorName = 'dose';
            elseif doseN ==1
                factorName = 'strain';
            end
            % make rmtable 
            rmtable = array2table(rmu,'VariableNames',{rmName});
            % make rmtable col names
            n = cellfun(@num2str,num2cell(rmu),'UniformOutput',0);
            n2 = repmat({rmName},numel(rmu),1);
            rmColNames = strjoinrows([n2,n],'');
            % ---------------------------------------

            % cycle through msrlist =======================================
            rmANOVAData = struct;
            for msri = 1:numel(msrlist)
                msr = msrlist{msri};
                % get data ++++++++
                d = array2table(S.(msr),'VariableNames',rmColNames);
                D = [F d];
                % -----------------

                % RMANOVA ================================================================
                textout = '';

                % rmanova multi-factor +++++++
                rmTerms = sprintf('%s%d-%s%d~%s',rmName,rmu(1),rmName,rmu(end),factorName);
                rm = fitrm(D,rmTerms,'WithinDesign',rmtable);
                t = anovan_textresult(ranova(rm), 0, 'pvlimit',pvlimit);
                textout = sprintf('%s\nRMANOVA(%s:%s):\n%s\n',textout,rmName,factorName,t);
                % ----------------------------

                % rmanova by groupnames ++++++++
                rmTerms = sprintf('%s%d-%s%d~%s',rmName,rmu(1),rmName,rmu(end),compName);
                rm = fitrm(D,rmTerms,'WithinDesign',rmtable);
                % -------------------------------

                % rmanova pairwise by groupnames ++++
                t = multcompare(rm,compName);
                % text output
                textout = sprintf('%s\nPosthoc(Tukey) by %s:\n',textout,compName);
                if isempty(t); 
                    textout = sprintf('%s\nAll comparison = n.s.\n',textout);
                else
                    t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha);
                    textout = sprintf('%s\n%s\n',textout,t);
                end

                % comparison by groups between taps +++++++++++
                % prepare comparision paris 
                g = unique(D.(compName));
                gg = pairwisecomp_getpairs(g);
                gpairs = strjoinrows(gg,' x ');
                % comparison
                t = multcompare(rm,compName,'By',rmName);
                %  keep only unique comparisons
                a = strjoinrows([t.([compName,'_1']) t.([compName,'_2'])],' x ');
                t(~ismember(a,gpairs),:) =[]; 
                % record
                    textout = sprintf('%s\n\nPosthoc(Tukey)%s by %s:\n',textout,rmName,compName); 
                if isempty(t); 
                    textout = sprintf('%s\nAll comparison = n.s.\n',textout);
                else
                    t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha);
                    textout = sprintf('%s\n%s\n',textout,t);
                end
                % -----------------------------------------------

                % comparison by taps between groups +++++++++++
                % comparison
                t = multcompare(rm,rmName,'By',compName);
            %     %%  keep only unique comparisons
            %     a = strjoinrows([t.('tap_1') t.('tap_2')],' x ');
            %     t(~ismember(a,gpairs),:) =[]; 
                % record
                    textout = sprintf('%s\n\nPosthoc(Tukey)%s by %s:\n',textout,rmName,compName); 
                if isempty(t); 
                    textout = sprintf('%s\nAll comparison = n.s.\n',textout);
                else
                    t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha);
                    textout = sprintf('%s\n%s\n',textout,t);
                end
                % -----------------------------------------------
                rmANOVAData.(msr) = textout;
                % =========================================================================
            end
            % =========================================================================

        end
        function DataTable = createDataTable(obj)
            DataLeg = table;
            DataLeg.groupname = obj.groupname;
            DataLeg.tap = obj.tap;
            A = array2table(obj.Response,'VariableNames',obj.Measures);
            DataTable = [DataLeg A];
        end
        function G  = graphicData(obj)
            msrlist = [obj.Measures];
            DataTable = createDataTable(obj);
            G = statsByGroup(DataTable,msrlist,'tap','groupname');
        end
        function fig = plotHabStd(obj,msr,visibleopt)
            % setting 
            gn = obj.gnu_genotype;
            gn = regexprep(gn,'_',' ');
            % generate graphic variables 
            G  = graphicData(obj);
            GData = G.(msr);
            X = GData.X;
            Y = GData.Y;
            E = GData.E;
            % sort N2 first
            gn = sortN2first(gn,gn);
            X = sortN2first(gn,X')';
            Y = sortN2first(gn,Y')';
            E = sortN2first(gn,E')';
            % figure 
            fig = graphHabCurvePck(X,Y,E,msr,gn,'visibleopt',visibleopt);
        end
        function D = pctHabData(obj,msr,t1,t2)
            % get habituated pct
            t = obj.tap;
            g = obj.groupname;
            m = obj.mwtid;
            M = obj.MWTDB;
            k= ismember(obj.Measures,msr);
            A = obj.Response(:,k);
            
            % get initial
            i = t==t1;
            j = t==t2;
            row1 =m(i);
            row2 = m(j);
            d1 = A(i);
            d2 = A(j);
            B = nan(numel(unique(m)),3);
            B(row1,1) = d1;
            B(row2,2) = d2;
            C = B(:,2)./B(:,1);
            C = 1-C;
            D = table;
            D.mwtid = unique(m);
            D.gname = M.groupname;
            D.pct = C;

        end
    end
    methods (Static)
        % functions
        function dose = findDose(grpname)
            a = grpname;
            b = regexpcellout(a,'(?<=(_))\d{1,}mM','match');
            b(cellfun(@isempty,b)) = {'0mM'};
            dose = b;
        end
        function strain = findStrain(grpname)
            a = grpname;
            strain = regexpcellout(a,'^[A-Z]{1,}\d{1,}','match');
        end
    end

end













