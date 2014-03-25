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
                % disect soma:
                obj.soma = obj.detectRegions(1);
                % disect axon:
                obj.axon = obj.detectRegions(2);
                % disect basal dendrites:
                obj.dends = obj.detectRegions(3);
                % disect apical dendrites:
                obj.apics = obj.detectRegions(4);
                
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
            for i= find( [ ((RAW_BTC==0) | (RAW_BTC>1)) &...
                    (cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == dendType) ] )
                
                branch = struct();
                %Add the terminal node to the list:
                branch(1).id= i;
                branch(1).x = obj.tree.X(i);
                branch(1).y = obj.tree.Y(i);
                branch(1).z = obj.tree.Z(i);
                branch(1).d = obj.tree.D(i);
                %         Recursively add the rest of the terminal branch nodes:
                branch = obj.recursionIn(branch,i);
                branches(ctr) = {branch(:)};
                ctr = ctr+1;
            end
            
            % if region is soma, exclude the root node!
            if(dendType == 1)
                branches(1) = [];
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
            tmp = sum(full(obj.tree.dA),2);
            brkState = all(tmp(2:end)) ;
        end
        function obj = hasPersistent(obj,stim,freq,dur)
            %returns true if responce similars persistent activity
            %duration in ms
            %IGNORES the stimulus duration!
            for i=stim:obj.tstop - dur
                if( (spike_count(obj,i,i+dur)/ ((obj.tstop - dur)/1000)) >= freq )
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
        function meanMEP = measureMEP(obj,rm,ra)
            % Function that scans for terminal nodes in adjacency matrix 
            % and recursively falls back to the soma, computing mean MEP.
            % As in:
            % Impact of Dendritic Size and Dendritic Topology on Burst
            % Firing in Pyramidal Cells, Van Ooyen, 2010
            
            MEP = @(b, rm, ra)sqrt(((b/2)*rm)/(2*ra));
            
            
            branches={};
            ctr = 1;
            terms = find( (sum(full(obj.tree.dA))==0) & (cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) ~= 1) )
            for i=find(sum(full(obj.tree.dA))==0)
                %            Add the terminal node to the list:
                branch = struct();
                 %Add the terminal node to the list:
                branch(1).id= i;
                branch(1).MEP = MEP(obj.tree.D(i),rm,ra);
                %         Recursively add the rest of the terminal branch nodes:
                branch = obj.recursionMEP(MEP,branch,i,rm,ra);
                branches(ctr) = {branch(:)};
                ctr = ctr+1;
            end
            meanMEP = mean(cellfun(@(x) x.MEP, branches));
        end
        function branch = recursionMEP(obj,MEP,branch,idx,rm,ra)
            %     Add the index of the parent node;
            tmpIdx = find(obj.tree.dA(idx,:) );
            tmpLoc = size(branch,2)+1;
            branch(tmpLoc).id=tmpIdx;
            branch(tmpLoc).MEP = MEP(obj.tree.D(tmpIdx),rm,ra);

            %     If the parent is continuoum point, continue the recursion:
            if branch(end).id ~= 1 
                branch = obj.recursionMEP(MEP,branch,branch(end).id,rm,ra);
            else
                %     If the parent is bifurcation, STOP
                return;
            end
        end
        
    end
    
    methods (Static)
        function tree = fixRnames(subtree)
            insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));
            rnums = cellfun(@(x) str2num(x),subtree.rnames);
            for i=1:4
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
        function diam = regionDiameter(reg)
            % region is a cell of structures containing segments
            for i=1:length(reg)
                for t=1:length(reg{i})-1
                    dtmp(t) = reg{i}(t).d; % is norm() the same?
                end
                diam{i} = dtmp;
            end
        end
    end
    
end

