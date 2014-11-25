classdef ncell
    %NCELL, class to hold neural cell info
    % eg voltage responce, position, 3D tree structure etc..
    % for bugs and comments: stamatiad.st at gmail.com
    
    properties
        mv; % Voltage responce in milivolts
        tstop; %Voltae responce duration in miliseconds
        dt; % Simulation's descrete steps per milisecond
        freq; % Overall frequency in Hz
        nspikes; % Overall number of spikes
        spikes; % Timings of spikes in miliseconds
        position; % 3d position (3d vector)
        clusterID; % Cluster ID that ncell belongs to
        stBins; % Spike train bins
        persistentActivity; % if persistent activity
        tree % Hierarchical tree (binary) structure
        soma % Array (cell) of structures for each section of the soma
        axon % Array (cell) of structures for each section of the axon
        dends % Array (cell) of structures for each dendritic section
        % (from branch point to branch point)
        apics % Array (cell) of structures for each apical section
        
    end
    
    methods %(Access = public)
        function obj = ncell(varargin)
            % Constructon of the class:
            % Call ncell(vresponce) to initialize .mv=[];
            % Call ncell(vresponce,dt) to initialize
            %             mv=[]; % Voltage responce in milivolts
            %             tstop=[]; %Voltae responce duration in miliseconds
            %             dt=[]; % Simulation's descrete steps per milisecond
            %             freq=[]; % Overall frequency in Hz
            %             nspikes=[]; % Overall number of spikes
            %             spikes=[]; % Timings of spikes in miliseconds
            % Call ncell(vresponce,dt,'swcFileName') to associate a
            % morphology also.
            
            
            % initialize:
            mv=[]; % Voltage responce in milivolts
            tstop=[]; %Voltae responce duration in miliseconds
            dt=[]; % Simulation's descrete steps per milisecond
            freq=[]; % Overall frequency in Hz
            nspikes=[]; % Overall number of spikes
            spikes=[]; % Timings of spikes in miliseconds
            position=[]; % 3d position (3d vector)
            clusterID=[]; % Cluster ID that ncell belongs to
            stBins=[]; % Spike train bins
            persistentActivity=[]; % if persistent activity
            tree=[]; % Hierarchical tree (binary) structure
            soma={}; % Array (cell) of structures for each section of the soma
            axon={}; % Array (cell) of structures for each section of the axon
            dends={}; % Array (cell) of structures for each dendritic section
            % (from branch point to branch point)
            apics={}; % Array (cell) of structures for each apical section
            
            if nargin > 0
                membrane = varargin{1};
                if(~isempty(membrane))
                    if(mod(length(varargin{1}),2))
                        membrane = varargin{1}(1:end-1);
                    end
                    obj.mv=membrane; % in milivolts
                end
            end
            
            if nargin > 1
                if(~isempty(obj.mv))
                    obj.tstop = length(membrane) / varargin{2};
                    obj.dt = varargin{2};
                    [number_of_spikes, spike_timing] = obj.spike_count(1,obj.tstop);
                    obj.freq = number_of_spikes / (obj.tstop / 1000);
                    obj.nspikes = number_of_spikes;
                    obj.spikes = spike_timing / varargin{2};
                end
            end
            if nargin  > 2
                % Load tree from swc file:
                obj.tree = load_tree(varargin{3});
                obj.tree = ncell.fixRnames(obj.tree);
                tmpRegs = unique(obj.tree.R);
                % xxx BUG TERASTEIWN DIASTASEWN: kati paizei me to soma..
                % disect soma:
                if find(tmpRegs==1)
                    obj.soma = obj.detectRegions(1);
                end
                % disect axon:
                if find(tmpRegs==2) 
                    obj.axon = obj.detectRegions(2);
                end
                % disect basal dendrites:
                if find(tmpRegs==3) 
                    obj.dends = obj.detectRegions(3);
                end
                % disect apical dendrites:
                if find(tmpRegs==4) 
                    obj.apics = obj.detectRegions(4);
                end
                
            end
            
        end
        function [number_of_spikes, spike_timing] = spike_count(obj,on,off)
            % Returns spikes by detecting zero-crossings that occur when the membrane
            % potential is positive.
            if on ~= 1
                on = on * obj.dt;
            end
            if off ~= 1
                off = off * obj.dt;
            end
            responce = obj.mv(on:off);
            if size(responce,2) > size(responce,1)
                responce = responce';
            end
            spike_timing = find( ([0;diff(sign(diff(responce)))<0;0] & [sign(responce)==1]) );
            number_of_spikes = length(spike_timing);
        end
        % Function to return branch IDXs and points
%         function branches = detectRegions(obj,dendType)
%             obj.tree = ncell.fixRnames(obj.tree);
% %             Verify trees structure ( must be binary tree!):
%             ver_tree (obj.tree);
% %             Extract BCT:
%             RAW_BTC = sum(full(obj.tree.dA),1) ;
% %             trace backwards terminal/branch points to extract a tree
% %             branch
%             
%             branches={};
%             ctr=1;
%             for i= find( [ ((RAW_BTC==0) | (RAW_BTC>1)) &...
%                     (cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == dendType) ] )
% 
% %                 Slow recursion
%                 branch = struct();
% %                 Add the terminal node to the list:
%                 branch(1).id= i;
%                 branch(1).x = obj.tree.X(i);
%                 branch(1).y = obj.tree.Y(i);
%                 branch(1).z = obj.tree.Z(i);
%                 branch(1).d = obj.tree.D(i);
% %                         Recursively add the rest of the terminal branch nodes:
%                 branch = obj.recursionIn(branch,i);
%                 branches(ctr) = {branch(:)};
%                 ctr = ctr+1;
%             end
%             
% %             if region is soma, exclude the root node!
%             if(dendType == 1)
%                 branches(1) = [];
%             end
%         end
        function branches = detectRegions(obj,dendType)
            obj.tree = ncell.fixRnames(obj.tree);
            % Verify trees structure ( must be binary tree!):
            ver_tree (obj.tree);
            % Extract BCT:
            RAW_BTC = sum(full(obj.tree.dA),1) ;
            % trace backwards terminal/branch points to extract a tree
            % branch
            
            branches={};
            ctr=1;
            Parents = ipar_tree(obj.tree);
            AL = len_tree(obj.tree);
            for i= find( [ ((RAW_BTC==0) | (RAW_BTC>1)) &...
                    (cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == dendType) ] )
                               
                iPar = unique([i,Parents(i,:)],'stable');
%                 iPar = [i,Parents(i,:)];
                iPar = iPar(iPar~=0);
                % if region contains only root, exclude it!
                if length(iPar) == 1
                    continue;
                end
                stopIdx = find(  RAW_BTC(iPar(2:end))>1 ) +1;
%                 if (stopIdx(1)==1) || (isempty(stopIdx))
%                     continue;
%                 else
                    iPar = iPar(1:stopIdx(1));
%                 end
                branch = repmat(struct('id',[],'x',[],'y',[],'z',[],'d',[],'l',[]), length(iPar),1);
                for l=1:length(iPar)
                    branch(l).id= iPar(l);
                    branch(l).x = obj.tree.X(iPar(l));
                    branch(l).y = obj.tree.Y(iPar(l));
                    branch(l).z = obj.tree.Z(iPar(l));
                    branch(l).d = obj.tree.D(iPar(l));
                    branch(l).l = AL(iPar(l));
                end
                branches(ctr) = {branch(:)};
                ctr = ctr+1;
            end
        end
        function branch = recursionIn(obj,branch,idx)
            % Let recursion kick in to save lines of code: continue until hit a
            % bifurcation node.
            
            % Add the index of the parent node;
            tmpIdx = find(obj.tree.dA(idx,:) );
            tmpLoc = size(branch,2)+1;
            branch(tmpLoc).id=tmpIdx;
            branch(tmpLoc).x = obj.tree.X(tmpIdx);
            branch(tmpLoc).y = obj.tree.Y(tmpIdx);
            branch(tmpLoc).z = obj.tree.Z(tmpIdx);
            branch(tmpLoc).d = obj.tree.D(tmpIdx);
            % If the parent is continuoum point, continue the recursion:
            if sum(obj.tree.dA(:,branch(end).id))==1
                branch = obj.recursionIn(branch,branch(end).id);
            else
                % If the parent is bifurcation, STOP
                return;
            end
        end
        function subtree = getSubtree(obj,strID)
            obj.tree = ncell.fixRnames(obj.tree);
            % decode names (SWC naming convention here)
            switch strID
                case 'soma'
                    sID = 1;
                case 'axon'
                    sID = 2;
                case 'dend'
                    sID = 3;
                case 'apic'
                    sID = 4;
                otherwise
                    error('Not known subregion name! Exiting getSubtree()');
            end
            
            % Keep branches indicated by subtreeID (ex basal == 3):
            keep = find(cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == sID);
            
            % initialize the temp tree:
            subtree.dA = obj.tree.dA(keep,keep);
            subtree.X = obj.tree.X(keep);
            subtree.Y = obj.tree.Y(keep);
            subtree.Z = obj.tree.Z(keep);
            subtree.D = obj.tree.D(keep);
            subtree.R = obj.tree.R(keep);
            subtree.rnames=obj.tree.rnames;
            subtree.name = [obj.tree.name,'_',strID];
            subtree = ncell.fixRoot(subtree);
        end
        function obj = swapSubtree(obj,subtree,strID)
            
            %given an empty subtree, exit:
            if( isempty(subtree.dA) )
                warning('EMPTY subtree: no swapping done.');
                return;
            end
            
            % Temprorary: deleting the root (ie soma) is not implemented
            if(strID == 'soma')
                warning('For swaping the root (ie soma) use swapSoma().');
                return;
            end
            
            % decode names (SWC naming convention here)
            switch strID
                case 'axon'
                    sID = 2;
                case 'dend'
                    sID = 3;
                case 'apic'
                    sID = 4;
                otherwise
                    error('Not known subregion name! Exiting swapSubtree()');
            end
            
            obj.tree = ncell.fixRnames(obj.tree);
            % Verify subtree structure ( must be binary tree!):
            ver_tree (subtree);
            
            % remove old subtree:
            % TREESTOOLBOX HAS BUG IN delete_tree()
            obj.tree = delete_tree(obj.tree,...
                find(cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == sID) );
            
            % Concatenate trees:
            somaXYZ = [obj.tree.X(1), obj.tree.Y(1), obj.tree.Z(1)];
            subtreeXYZ = [subtree.X(1), subtree.Y(1), subtree.Z(1)];
            
            obj.tree = cat_tree(obj.tree,tran_tree(subtree,somaXYZ-subtreeXYZ),1,1);
            obj.tree = elim0_tree(obj.tree);
            
            % Update properties:
            % disect soma:
            obj.soma = obj.detectRegions(1);
            % disect axon:
            obj.axon = obj.detectRegions(2);
            % disect basal dendrites:
            obj.dends = obj.detectRegions(3);
            % disect apical dendrites:
            obj.apics = obj.detectRegions(4);
            
        end
        function obj = swapSoma(obj,subtree)
            
            %given an empty subtree, exit:
            if( isempty(subtree.dA) )
                warning('EMPTY subtree: no swapping done.');
                return;
            end
            
            if( subtree.name(end-3:end) ~= 'soma' )
                warning('Inserting subtree in soma location! Exiting swapSoma()');
                return;
            end
            
            obj.tree = ncell.fixRnames(obj.tree);
            % Verify subtree structure ( must be binary tree!):
            ver_tree (subtree);
            
            % get new tree sub-component trees:
            tmpRoot = subtree;
            tmpAxon = obj.getSubtree('axon');
            tmpDend = obj.getSubtree('dend');
            tmpApic = obj.getSubtree('apic');
            
            % move axon, basals and apicals to new soma position:
            % and concatenate trees:
            tmpName = obj.tree.name;
            obj.tree = tmpRoot;
            obj.tree.name = tmpName;
            %             tmpAxon.name = tmpName;
            %             tmpDend.name = tmpName;
            %             tmpApic.name = tmpName;
            
            somaXYZ = [tmpRoot.X(1), tmpRoot.Y(1), tmpRoot.Z(1)];
            if( ~isempty(tmpAxon.dA) )
                subtreeXYZ = [tmpAxon.X(1), tmpAxon.Y(1), tmpAxon.Z(1)];
                obj.tree = cat_tree(obj.tree,tran_tree(tmpAxon,somaXYZ-subtreeXYZ),1,1);
                %                 obj.tree = ncell.fixRoot(obj.tree);
            end
            if( ~isempty(tmpDend.dA) )
                subtreeXYZ = [tmpDend.X(1), tmpDend.Y(1), tmpDend.Z(1)];
                obj.tree = cat_tree(obj.tree,tran_tree(tmpDend,somaXYZ-subtreeXYZ),1,1);
                %                 obj.tree = ncell.fixRoot(obj.tree);
                %                 obj.tree = redirect_tree (obj.tree, 1);
            end
            if( ~isempty(tmpApic.dA) )
                subtreeXYZ = [tmpApic.X(1), tmpApic.Y(1), tmpApic.Z(1)];
                obj.tree = cat_tree(obj.tree,tran_tree(tmpApic,somaXYZ-subtreeXYZ),1,1);
                %                 obj.tree = ncell.fixRoot(obj.tree);
            end
            
            % IMPORTANT: delete the intermediate branches created by
            % fixRoot()
            obj.tree = elim0_tree(obj.tree);
%             subtree = ncell.fixRoot(subtree);
            
            % Update properties:
            % disect soma:
            obj.soma = obj.detectRegions(1);
            % disect axon:
            obj.axon = obj.detectRegions(2);
            % disect basal dendrites:
            obj.dends = obj.detectRegions(3);
            % disect apical dendrites:
            obj.apics = obj.detectRegions(4);
            
        end
        function brkState = isBroken(obj)
            % If returns one if one or more nodes are missing their parrents:
            tmp = sum(full(obj.tree.dA),2);
            % exclude first node that has no parent:
            brkState = ~all(tmp(2:end)) ;
        end
        function maxfurc = maxBranching(obj)
            % Return broken state of morphology:
            % returned int indicates the max tri(and up)-furcations
            % in the morphology. So 'normal' state is brkState==1 (soma).
            % More trifurcations might indicate that something went wrong.
            % the tools to check the morphology for errors are incomplete
            % or non existent.
            maxfurc = max(sum(full(obj.tree.dA),1)) ;
        end
        function obj = hasPersistent(obj,stim,freq,dur)
            %returns true if responce similars persistent activity
            %duration in ms
            %IGNORES the stimulus duration!
            %Pontiako implementation.. Logika, prwta briskeis ta spikes mia fora kai
            %meta dokimazeis gia ka8e para8yro poia briskontai mesa (dah)
            for i=stim:obj.tstop - dur
                if( (spike_count(obj,i,i+dur)/ ((dur)/1000)) >= freq )
                    obj.persistentActivity = 1;
                    return;
                end
            end
            obj.persistentActivity = 0;
        end
        function obj = binSpikes(obj,binSize)
            % Function to return bins containing the spikes of the cell
            % The binSize to be given is in miliseconds
            
            % check if obj has freq (must have due to constructor):
            if(isempty(obj.freq))
                error('spike_count have not run forncell obj!');
            end
            obj.stBins = histc(obj.spikes,1:binSize:(obj.tstop*obj.dt));
        end
        function meanMEP = dendMEP(obj,sID,rm_s, rm_e,ra,rm_dhalf,rm_steep, varargin)
            % Function that scans for terminal nodes in adjacency matrix 
            % and recursively falls back to the soma, computing mean MEP.
            % As in:
            % Impact of Dendritic Size and Dendritic Topology on Burst
            % Firing in Pyramidal Cells, Van Ooyen, 2010
            % sID : compute mean MEP for a swc section (soma,axon,dend,etc)
            
            MEP = @(l, b, rm, ra)l/sqrt(((b/2)*rm)/(2*ra));
            rm_sig = @(rm_s, rm_e, dh, st, dist)rm_s + (rm_e - rm_s)/(1.0 + exp((dh-dist)/st));
            
            branches={};
            ctr = 1;
            iC = find(cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == sID); %Total number of compartmets of interest
            mC = length(iC);
            dist_tmp = Pvec_tree(obj.tree);
            if nargin > 7
                RMs = repmat(rm_s/1.5,[mC,1]); %xxx
            else
                RMs = arrayfun(rm_sig, repmat(rm_s,[mC,1]), repmat(rm_e,[mC,1]), repmat(rm_dhalf,[mC,1]), repmat(rm_steep,[mC,1]),dist_tmp(iC)) ;%rm_sig(rm_s, rm_e, rm_dhalf, rm_steep, dist_tmp(i));
            end
            AL = len_tree(obj.tree);
            Diam = obj.tree.D;
            Parents = ipar_tree(obj.tree);
            % Calculate MEP for all the compartments of interest (to avoid
            % recomputation for some compartments:
            MEPs = arrayfun(MEP, AL(iC), Diam(iC), RMs, repmat(ra,[mC,1]));
                
            terms = find( (sum(full(obj.tree.dA))==0) & (cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == sID) );
            for i=terms %find(sum(full(obj.tree.dA))==0)
                %         Recursively add the rest of the terminal branch nodes:
                % EDIT: avoid recursion for it is to slow in Matlab:
%                 branch = obj.recursionMEP(MEP,rm_sig,branch,i,rm_s,
%                 rm_e,ra);
                iPar = unique([i,Parents(i,:)],'first');
                iPar = iPar(iPar~=0);
                %preallocate for speed:
                branch = repmat(struct('id',[],'MEP',[]),length(iPar),1);
                for l=1:length(iPar)
                    branch(l).id = iPar(l);
                    branch(l).MEP = MEPs( iC==iPar(l) );
                end
                branches(ctr) = {branch(:)};
                ctr = ctr+1;
            end
%             meanMEP = mean(cellfun(@(x) x.MEP, branches));
              meanMEP = mean(cellfun(@(cb) mean(cell2mat(arrayfun(@(x) x.MEP, cb,'uniformoutput',false))), branches))
        end
    end
    
    methods (Static)
        function tree = fixRnames(subtree)
            insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));
            rnums = cellfun(@(x) str2num(x),subtree.rnames);
            for i=1:3
                if(~(rnums(i)==i))
                    rnums=insert(i,rnums,i-1);
                    tmpR = find(cellfun(@(x) str2num(x),subtree.rnames((subtree.R))) >= i) ;
                    subtree.R(tmpR) = subtree.R(tmpR) + 1;
                    subtree.rnames = cellfun(@(x) {num2str(x)},mat2cell(rnums,1,ones(1,length(rnums))));
                end
            end
            tree = subtree;
        end
        function outtree = fixRoot(intree)
            % fix root (connect orphan sections to the root):
            orphans = ~logical(sum(full(intree.dA),2));
            intree.dA(find(orphans(2:end))+1,1) = 1;
            outtree=intree;
        end
        function vol = regionVolume(reg)
            volum = @(l,d) l*pi*(d/2)^2 ;
            % region is a cell of structures containing segments
            for i=1:length(reg)
                for t=1:length(reg{i})-1
                    Ls = [reg{i}(t).x,reg{i}(t).y,reg{i}(t).z];
                    Le = [reg{i}(t+1).x,reg{i}(t+1).y,reg{i}(t+1).z];
                    vtmp(t) = volum(sqrt(sum((Ls-Le).^2)),reg{i}(t).d); % is norm() the same?
                end
                vol{i} = vtmp;
            end
        end
        function len = regionLength(reg)
            % region is a cell of structures containing segments
            for i=1:length(reg)
                for t=1:length(reg{i})-1
                    Ls = [reg{i}(t).x,reg{i}(t).y,reg{i}(t).z];
                    Le = [reg{i}(t+1).x,reg{i}(t+1).y,reg{i}(t+1).z];
                    ltmp(t) = sqrt(sum((Ls-Le).^2)); % is norm() the same?
                end
                len{i} = ltmp;
            end
        end
        function diam = regionDiameter(reg,nseg)
            % region is a cell of structures containing segments
            for i=1:length(reg) % No of sections
                for t=1:length(reg{i})-1 %No of pt3d
                    dtmp(t) = reg{i}(t).d; % is norm() the same?
                    ltmp(t) = reg{i}(t).l;
                end
                % As in NEURON source code:
                diam(i) = sum((dtmp(1:end-1) + dtmp(2:end)) .* ltmp(1:end-1)) * 0.5 / (sum(ltmp(1:end-1))/nseg) ;
%                 diam{i} = dtmp;
            end

%             if nvarg < 2
%                 diam=cellfun(@(as) arrayfun(@(f) f.d,as,'uniformoutput',false),reg,'uniformoutput',false);
%                 % remove the last diam as is the parent's diam
%                 
%             else
%                 
%             end
        end
        
    end
    
end

% function meanMEP = dendMEP(obj,sID,rm_s, rm_e,ra,rm_dhalf,rm_steep)
%             % Function that scans for terminal nodes in adjacency matrix 
%             % and recursively falls back to the soma, computing mean MEP.
%             % As in:
%             % Impact of Dendritic Size and Dendritic Topology on Burst
%             % Firing in Pyramidal Cells, Van Ooyen, 2010
%             % sID : compute mean MEP for a swc section (soma,axon,dend,etc)
%             
%             MEP = @(l, b, rm, ra)l/sqrt(((b/2)*rm)/(2*ra));
%             rm_sig = @(rm_s, rm_e, dh, st, dist)rm_s + (rm_e - rm_s)/(1.0 + exp((dh-dist)/st));
%             
%             branches={};
%             ctr = 1;
%             terms = find( (sum(full(obj.tree.dA))==0) & (cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == sID) );
%             for i=terms%find(sum(full(obj.tree.dA))==0)
%                 %            Add the terminal node to the list:
%                 branch = struct();
%                  %Add the terminal node to the list:
%                 branch(1).id= i;
%                 dist_tmp = Pvec_tree(obj.tree);
% %                 rm=rm_sig(rm_s, rm_e, rm_dhalf, rm_steep, dist_tmp(i));
%                 rm = 1.5/rm_s; %xxx
%                 branch(1).MEP = MEP(obj.tree.D(i),rm,ra);
%                 %         Recursively add the rest of the terminal branch nodes:
%                 branch = obj.recursionMEP(MEP,rm_sig,branch,i,rm_s, rm_e,ra);
%                 branches(ctr) = {branch(:)};
%                 ctr = ctr+1;
%             end
%             meanMEP = mean(cellfun(@(x) x.MEP, branches));
%         end
%         function branch = recursionMEP(obj,MEP,rm_sig,branch,idx,rm_s, rm_e,ra)
%             %     Add the index of the parent node;
%             tmpIdx = find(obj.tree.dA(idx,:) );
%             tmpLoc = size(branch,2)+1;
%             branch(tmpLoc).id=tmpIdx;
%             dist_tmp = Pvec_tree(obj.tree);
%             rm=rm_sig(rm_s, rm_e, 10, 5, dist_tmp(tmpIdx));
%             branch(tmpLoc).MEP = MEP(obj.tree.D(tmpIdx),rm,ra);
% 
%             %     If the parent is continuoum point, continue the recursion:
%             tmpBO = BO_tree(obj.tree);
%             if tmpBO(branch(end).id) >= sum(cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == 1)
%                 branch = obj.recursionMEP(MEP,rm_sig,branch,branch(end).id,rm_s, rm_e,ra);
%             else
%                 %     If the parent is bifurcation, STOP
%                 return;
%             end
%         end
        