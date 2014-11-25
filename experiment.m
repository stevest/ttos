% Add needed toolboxes to the path:
% NOTE: minor bugs changed in TREESTOOLBOX, so use the provided one.
% For comments/bugs stamatiad.st at gmail dot com

addpath('~/Documents/MATLAB/Dendrites/');
addpath(genpath('~/Documents/MATLAB/Dendrites/TREES1.15'));
runpath = '~/Documents/frontiers_in_cellular_neuroscience_rev/experiment/network/' ;

% Load morphologies names from Smith's lab
morphDir = 'SmithMorphologies';
names = dir(morphDir);
names = {names(~[names.isdir]).name};
names = names(cellfun(@(x) strcmp(x(end-3:end),'.swc'), names));

% % We are using A LOT of recursions (as much as the length of points of
% % swc files), so set limit accordingly to avoid warnings:
% set(0,'RecursionLimit',5000);

% Load all morphologies as ncell objects:
for i=1:length(names)
    tic;
    PC(i) =  ncell([],1,[morphDir,'/',names{i}]);
    t(i) = toc
end
mean(t)

% % Save morphologies: (gia ta paidia):
% for i=1:length(names)
%     tmpCell = PC(i);
%     tmpCell.tree = sort_tree(tmpCell.tree);
%     % Save as .hoc file using custom function:
%     neuron_tree_mod(tmpCell.tree,['~/Desktop/blah','/',tmpCell.tree.name,'.hoc'],1,'');
% end

% Get soma volume of each ncell:
for i=1:length(PC)
    tic;
    volume_soma(i) = sum( cellfun(@(x) (x),ncell.regionVolume(PC(i).soma)) );
    t(i) = toc
end
mean(t)

% Poorly implemented median:
[~,idx ]= sort(volume_soma);
tmp = volume_soma(idx);
tmp(round(length(tmp)/2));
median_volume_soma_idx = find(tmp(round(length(tmp)/2))==volume_soma);

median_soma = PC(median_volume_soma_idx).getSubtree('soma');

% find the simplest/more complex apical tree (based in bifurcation order):
for i=1:length(PC)
    BO(i) = max(BO_tree(PC(i).getSubtree('apic')));
end
[~,simple_apical_idx] = min(BO);
[~,complex_apical_idx] = max(BO);

simple_apical = PC(simple_apical_idx).getSubtree('apic');

NULL_axon = struct('dA',[],'X',[],'Y',[],'Z',[],'D',[],'R',[],'rnames',{},'name','');

%% Create all the possible combinations of morphologies: (a lot of cells!)
modif_dir = 'modified';
system(['mkdir ',modif_dir]);


for i=1:length(PC)
    i
        % Swap soma:
        tmpCell = PC(i).swapSoma(median_soma);
        
        % Delete axon compartments:
        tmpCell.tree = delete_tree(tmpCell.tree,...
            find(cellfun(@str2num,{tmpCell.tree.rnames{tmpCell.tree.R}}) == 2) );
        % Swap basal:
        tmpCell = tmpCell.swapSubtree(PC(i).getSubtree('dend'), 'dend');
        % Swap apical:
        tmpCell = tmpCell.swapSubtree(simple_apical, 'apic');
        % Save as .hoc file using custom function
        tmpCell.tree = sort_tree(tmpCell.tree);
        neuron_tree_mod(tmpCell.tree,['~/Desktop/blah','/',tmpCell.tree.name,'.hoc'],1,'');
end





%% Find a morphology with a reconstructed axon:
% for i=1:length(PC)
%     if( ~isempty(PC(i).axon) )
%         median_axon = PC(i).getSubtree('axon');
%         median_axon_idx = i;
%         break;
%     end
% end

% Override Kiki's axon:
median_axon = PC(59).getSubtree('axon');
median_axon_idx = 59;

%% find the simplest/more complex apical tree (based in bifurcation order):
for i=1:length(PC)
    BO(i) = max(BO_tree(PC(i).tree));
end
[~,simple_apical_idx] = min(BO);
[~,complex_apical_idx] = max(BO);
%% Create all the possible combinations of morphologies: (a lot of cells!)
all_dir = 'all_possible';
system(['mkdir ',all_dir]);


for i=1:length(PC)
    for j=1:length(PC)
        [i,j]
        % Swap soma:
        tmpCell = PC(i).swapSoma(median_soma);
        % Swap axon:
        tmpCell = tmpCell.swapSubtree(median_axon,'axon');
        % Swap basal:
        tmpCell = tmpCell.swapSubtree(PC(i).getSubtree('dend'), 'dend');
        % Swap apical:
        tmpCell = tmpCell.swapSubtree(PC(j).getSubtree('apic'), 'apic');
        % change name to keep track of whats is going on:
        tmpCell.tree.name = sprintf('comb_pc_b%d_a%d',i,j) ;
        % Save as .hoc file using custom function
        tmpCell.tree = sort_tree(tmpCell.tree);
%         neuron_tree_mod(tmpCell.tree,[all_dir,'/',tmpCell.tree.name,'.hoc'],1,'');
        neuron_tree_mod(tmpCell.tree,['~/Desktop/PFC_Alexandra/morphology/all_possible','/',tmpCell.tree.name,'.hoc'],1,'');
        all_comb_PC(i,j) = tmpCell;
    end
end
%% Save simple apical morphologies as .hoc files (NEURON)
% simple_dir = 'smith_simple_apical';
% system(['mkdir ',simple_dir]);
% simple_apical = PC(simple_apical_idx).getSubtree('apic');
% 
% for i=1:length(PC)
%     % Swap apical:
%     tmpCell = PC(i).swapSubtree(simple_apical,'apic');
%     % Swap axon:
%     tmpCell = tmpCell.swapSubtree(median_axon,'axon');
%     % Swap soma:
%     tmpCell = tmpCell.swapSoma(median_soma);
%     % Save as .hoc file using custom function
% %     neuron_tree_mod(tmpCell.tree,[simple_dir,'/',tmpCell.tree.name,'.hoc'],1,'');
%     simple_PC(i) = tmpCell;
% end

clear t;
% complex_apical_idx = 59; % Override: Kiki's cell
simple_apical_idx = 22; % Override: Apical that gives 50% IB percentage..
simple_dir = 'smith_simple_apical_kiki_rev';
system(['mkdir ',simple_dir]);
simple_apical = PC(simple_apical_idx).getSubtree('apic'); %COMPLEX APIC ID=47,46-3.CNG.swc

for i=1:length(PC)
    tic;
    % Swap apical:
    tmpCell = PC(i).swapSubtree(simple_apical,'apic');
    toc
    % Swap axon:
    tmpCell = tmpCell.swapSubtree(median_axon,'axon');
    toc
    % Swap soma:
    tmpCell = tmpCell.swapSoma(median_soma);
    toc
    % Save as .hoc file using custom function
    neuron_tree_mod(tmpCell.tree,[simple_dir,'/',tmpCell.tree.name,'.hoc'],1,'');
    simple_PC(i) = tmpCell;
    t(i)=toc
end
mean(t)


%% Save complex apical morphologies as .hoc files (NEURON)
clear t;
% complex_apical_idx = 59; % Override: Kiki's cell
complex_apical_idx = 110; % Override: Apical that gives 50% IB percentage..
complex_dir = 'smith_complex_apical_kiki_rev';
system(['mkdir ',complex_dir]);
complex_apical = PC(complex_apical_idx).getSubtree('apic'); %COMPLEX APIC ID=47,46-3.CNG.swc

for i=1:length(PC)
    tic;
    % Swap apical:
    tmpCell = PC(i).swapSubtree(complex_apical,'apic');
    toc
    % Swap axon:
    tmpCell = tmpCell.swapSubtree(median_axon,'axon');
    toc
    % Swap soma:
    tmpCell = tmpCell.swapSoma(median_soma);
    toc
    % Save as .hoc file using custom function
    neuron_tree_mod(tmpCell.tree,[complex_dir,'/',tmpCell.tree.name,'.hoc'],1,'');
    complex_PC(i) = tmpCell;
    t(i)=toc
end
mean(t)
%% Save same basal morphologies as .hoc files (NEURON)
sameb_dir = 'smith_same_basal';
system(['mkdir ',sameb_dir]);
same_basal_idx = 59; % Override: Kiki's cell
smith_same_basal = PC(same_basal_idx).getSubtree('dend'); %SAMEBASAL ID=18, 35-2.CNG.swc

for i=1:length(PC)
    % Swap apical:
    tmpCell = PC(i).swapSubtree(smith_same_basal,'dend');
    % Swap axon:
    tmpCell = tmpCell.swapSubtree(median_axon,'axon');
    % Swap soma:
    tmpCell = tmpCell.swapSoma(median_soma);
    % Save as .hoc file using custom function
    neuron_tree_mod(tmpCell.tree,[sameb_dir,'/',tmpCell.tree.name,'.hoc'],1,'');
    complex_PC_SB(i) = tmpCell;
end


%% Save same soma morphologies as .hoc files (NEURON)
sames_dir = 'smith_same_soma';
system(['mkdir ',sames_dir]);

for i=1:length(PC)
    % Swap axon:
    tmpCell = PC(i).swapSubtree(median_axon,'axon');
    % Swap soma:
    tmpCell = tmpCell.swapSoma(median_soma);
    % Save as .hoc file using custom function
    neuron_tree_mod(tmpCell.tree,[sames_dir,'/',tmpCell.tree.name,'.hoc'],1,'');
    complex_PC_SS(i) = tmpCell;
end

%% Make data for Alexandra. same apical
clear t;
volum = @(l,d) l*pi*(d/2)^2 ;
morphological_data = {} ;
% list_cells = {} ;
contin = 1;
for i=1:length(PC)
    i
    tic;
%     list_cells{i,1} = complex_PC(i).tree.name;
%     list_cells{i,2} = length(complex_PC(i).dends);
    
    morphological_data(contin:contin+length(complex_PC(i).dends)-1,1) =...
        repmat({complex_PC(i).tree.name},[length(complex_PC(i).dends),1]);
    
%     morphological_data(contin:contin+length(complex_PC(i).dends)-1,3) =...
%         mat2cell(cellfun(@(x) mean(x),ncell.regionDiameter(complex_PC(i).dends))',ones(length(complex_PC(i).dends),1),[1]);

%     morphological_data(contin:contin+length(complex_PC(i).dends)-1,3) =...
%         mat2cell(ncell.regionDiameter(complex_PC(i).dends,1)',ones(length(complex_PC(i).dends),1),[1])

    tmptmp = load([runpath,sprintf('Attributes/Attributes_%d.txt',i-1)]) ;
    morphological_data(contin:contin+length(complex_PC(i).dends)-1,3) =...
        mat2cell(tmptmp(:,1),ones(length(complex_PC(i).dends),1),[1]);
    
%     morphological_data(contin:contin+length(complex_PC(i).dends)-1,4) =...
%         mat2cell(cellfun(@(x) sum(x),ncell.regionLength(complex_PC(i).dends))',ones(length(complex_PC(i).dends),1),[1]);
    
    morphological_data(contin:contin+length(complex_PC(i).dends)-1,4) =...
        mat2cell(tmptmp(:,2),ones(length(complex_PC(i).dends),1),[1]);
    
%     morphological_data(contin:contin+length(complex_PC(i).dends)-1,5) =...
%         mat2cell(cellfun(@(x) sum(x),ncell.regionVolume(complex_PC(i).dends))',ones(length(complex_PC(i).dends),1),[1]);
    
    morphological_data(contin:contin+length(complex_PC(i).dends)-1,5) =...
        mat2cell(arrayfun(@(l,d) volum(l,d),tmptmp(:,2) ,tmptmp(:,1) ),ones(length(complex_PC(i).dends),1),[1]);
    
    % replace on cell # the diameter value:
%     D_ID = (cellfun(@(x) arrayfun(@(s) s.id, x(1:end-1))' ,complex_PC(i).dends,'uniformoutput',false))'  ;
D_=[];
    D_ = cell2mat([(cellfun(@(x) arrayfun(@(s) s.id, x(1:end-1)) ,complex_PC(i).dends,'uniformoutput',false))' , ...
(cellfun(@(x) arrayfun(@(s) s.d, x(1:end-1)) ,complex_PC(i).dends,'uniformoutput',false))'  ]) ;
    
%     D_Size=(cellfun(@(x) length(x(1:end-1)) ,complex_PC(i).dends,'uniformoutput',false))'  ;
%     NewD = 
    complex_PC(i).tree.D(D_(:,1)) = D_(:,2) ;
    
    morphological_data(contin:contin+length(complex_PC(i).dends)-1,6) =...
        repmat({complex_PC(i).dendMEP(3,30000,15000,100,10,5)},[length(complex_PC(i).dends),1]); % for basal only
    t(i)=toc
    contin = contin + length(complex_PC(i).dends);
end

%% UNIFORM in basal
clear t;
morphological_data_uni = {} ;
% list_cells = {} ;
contin = 1;
for i=1:length(PC)
    i
    tic;
%     list_cells{i,1} = complex_PC(i).tree.name;
%     list_cells{i,2} = length(complex_PC(i).dends);
    
    morphological_data_uni(contin:contin+length(complex_PC(i).dends)-1,1) =...
        repmat({complex_PC(i).tree.name},[length(complex_PC(i).dends),1]);
    
%     morphological_data_uni(contin:contin+length(complex_PC(i).dends)-1,3) =...
%         mat2cell(cellfun(@(x) mean(x),ncell.regionDiameter(complex_PC(i).dends))',ones(length(complex_PC(i).dends),1),[1]);
%     
%     morphological_data_uni(contin:contin+length(complex_PC(i).dends)-1,4) =...
%         mat2cell(cellfun(@(x) sum(x),ncell.regionLength(complex_PC(i).dends))',ones(length(complex_PC(i).dends),1),[1]);
%     
%     morphological_data_uni(contin:contin+length(complex_PC(i).dends)-1,5) =...
%         mat2cell(cellfun(@(x) sum(x),ncell.regionVolume(complex_PC(i).dends))',ones(length(complex_PC(i).dends),1),[1]);
    
    morphological_data_uni(contin:contin+length(complex_PC(i).dends)-1,6) =...
        repmat({complex_PC(i).dendMEP(3,30000,15000,100,10,5,'u')},[length(complex_PC(i).dends),1]); % for basal only
    t(i) = toc
    contin = contin + length(complex_PC(i).dends);
end

%% Same basal dendrite, different apicals. 
% Make data for Alexandra.
clear t;
morphological_data_sb = {} ;
% list_cells_sb = {} ;
contin = 1;
for i=1:length(PC)
    i
    tic;
%     list_cells_sb{i,1} = complex_PC_SB(i).tree.name;
%     list_cells_sb{i,2} = length(complex_PC_SB(i).apics);
    
    morphological_data_sb(contin:contin+length(complex_PC_SB(i).apics)-1,1) =...
        repmat({complex_PC_SB(i).tree.name},[length(complex_PC_SB(i).apics),1]);
    
%     morphological_data_sb(contin:contin+length(complex_PC_SB(i).apics)-1,3) =...
%         mat2cell(cellfun(@(x) mean(x),ncell.regionDiameter(complex_PC_SB(i).apics))',ones(length(complex_PC_SB(i).apics),1),[1]);

    tmptmp = load([runpath,sprintf('Attributes/Attributes_%d.txt',i-1)]) ;
    morphological_data_sb(contin:contin+length(complex_PC_SB(i).apics)-1,3) =...
        mat2cell(tmptmp(:,1),ones(length(complex_PC_SB(i).apics),1),[1]);
    
%     morphological_data_sb(contin:contin+length(complex_PC_SB(i).apics)-1,4) =...
%         mat2cell(cellfun(@(x) sum(x),ncell.regionLength(complex_PC_SB(i).apics))',ones(length(complex_PC_SB(i).apics),1),[1]);
     morphological_data_sb(contin:contin+length(complex_PC_SB(i).apics)-1,4) =...
        mat2cell(tmptmp(:,2),ones(length(complex_PC_SB(i).apics),1),[1]);
    
%     morphological_data_sb(contin:contin+length(complex_PC_SB(i).apics)-1,5) =...
%         mat2cell(cellfun(@(x) sum(x),ncell.regionVolume(complex_PC_SB(i).apics))',ones(length(complex_PC_SB(i).apics),1),[1]);

    morphological_data_sb(contin:contin+length(complex_PC_SB(i).apics)-1,5) =...
        mat2cell(arrayfun(@(l,d) volum(l,d),tmptmp(:,2) ,tmptmp(:,1) ),ones(length(complex_PC_SB(i).apics),1),[1]);
    
    D_=[];
    D_ = cell2mat([(cellfun(@(x) arrayfun(@(s) s.id, x(1:end-1)) ,complex_PC_SB(i).apics,'uniformoutput',false))' , ...
        (cellfun(@(x) arrayfun(@(s) s.d, x(1:end-1)) ,complex_PC_SB(i).apics,'uniformoutput',false))'  ]) ;
    complex_PC_SB(i).tree.D(D_(:,1)) = D_(:,2) ;
    
    morphological_data_sb(contin:contin+length(complex_PC_SB(i).apics)-1,6) =...
        repmat({complex_PC_SB(i).dendMEP(4,30000,15000,100, 300, 50)},[length(complex_PC_SB(i).apics),1]); % for basal only
    t(i)=toc
    contin = contin + length(complex_PC_SB(i).apics);
end

%% Same soma, different apicals&basals. 
% Make data for Alexandra.
morphological_data_ss = {} ;
% list_cells_ss = {} ;
contin = 1;
for i=1:length(PC)
    i
%     list_cells_ss{i,1} = complex_PC_SS(i).tree.name;
%     list_cells_ss{i,2} = length(complex_PC_SS(i).apics);
    
    morphological_data_ss(contin:contin+length(complex_PC_SS(i).apics)-1,1) =...
        repmat({complex_PC_SS(i).tree.name},[length(complex_PC_SS(i).apics),1]);
    
    morphological_data_ss(contin:contin+length(complex_PC_SS(i).apics)-1,3) =...
        mat2cell(cellfun(@(x) mean(x),ncell.regionDiameter(complex_PC_SS(i).apics))',ones(length(complex_PC_SS(i).apics),1),[1]);
    
    morphological_data_ss(contin:contin+length(complex_PC_SS(i).apics)-1,4) =...
        mat2cell(cellfun(@(x) sum(x),ncell.regionLength(complex_PC_SS(i).apics))',ones(length(complex_PC_SS(i).apics),1),[1]);
    
    morphological_data_ss(contin:contin+length(complex_PC_SS(i).apics)-1,5) =...
        mat2cell(cellfun(@(x) sum(x),ncell.regionVolume(complex_PC_SS(i).apics))',ones(length(complex_PC_SS(i).apics),1),[1]);
    
    morphological_data_ss(contin:contin+length(complex_PC_SS(i).apics)-1,6) =...
        repmat({complex_PC_SS(i).dendMEP(3,30000,15000,100, 300,50)},[length(complex_PC_SS(i).apics),1]); % for basal only
    
    contin = contin + length(complex_PC_SS(i).apics);
end

%% Convert morphologies to HOC format:
cd('Cvapp');
for i=1:length(PC)
    unix(['./run.sh ' ,'../',simple_dir,'/',PC(i).tree.name,'.swc']);
end

for i=1:length(PC)
    unix(['./run.sh ' ,'../',complex_dir,'/',PC(i).tree.name,'.swc']);
end

%% Append morphological lists:

for i=1:length(PC)
    fid=fopen([simple_dir,'/',PC(i).tree.name,'.hoc'],'a');
    
    fprintf(fid,'\n');
    fprintf(fid,'objref all, somatic, axonal, basal, apical\n');
    fprintf(fid,'proc subsets() { local i\n');
    fprintf(fid,'objref all, somatic, axonal, basal, apical\n');
    fprintf(fid,'all = new SectionList()\n');
    fprintf(fid,'soma all.append()\n');
    fprintf(fid,sprintf('for i=0, %d axon[i] all.append()\n',length(PC(i).axon)));
    fprintf(fid,sprintf('for i=0, %d dend[i] all.append()\n',length(PC(i).dends)));
    fprintf(fid,sprintf('for i=0, %d apic[i] all.append()\n',length(PC(i).apics)));
    
    fprintf(fid,'somatic = new SectionList()\n');
    fprintf(fid,'soma somatic.append()\n');
    
    fprintf(fid,'axonal = new SectionList()\n');
    fprintf(fid,sprintf('for i=0, %d axon[i] axonal.append()\n',length(PC(i).axon)));
    
    fprintf(fid,'basal = new SectionList()\n');
    fprintf(fid,sprintf('for i=0, %d dend[i] basal.append()\n',length(PC(i).dends)));
    
    fprintf(fid,'apical = new SectionList()\n');
    fprintf(fid,sprintf('for i=0, %d apic[i] apical.append()\n',length(PC(i).apics)));
    
    fprintf(fid,'}\n');
    fclose(fid);
end


for i=1:length(PC)
    fid=fopen([complex_dir,'/',PC(i).tree.name,'.hoc'],'a');
    
    fprintf(fid,'\n');
    fprintf(fid,'objref all, somatic, axonal, basal, apical\n');
    fprintf(fid,'proc subsets() { local i\n');
    fprintf(fid,'objref all, somatic, axonal, basal, apical\n');
    fprintf(fid,'all = new SectionList()\n');
    fprintf(fid,'soma all.append()\n');
    fprintf(fid,sprintf('for i=0, %d axon[i] all.append()\n',length(PC(i).axon)));
    fprintf(fid,sprintf('for i=0, %d dend[i] all.append()\n',length(PC(i).dends)));
    fprintf(fid,sprintf('for i=0, %d apic[i] all.append()\n',length(PC(i).apics)));
    
    fprintf(fid,'somatic = new SectionList()\n');
    fprintf(fid,'soma somatic.append()\n');
    
    fprintf(fid,'axonal = new SectionList()\n');
    fprintf(fid,sprintf('for i=0, %d axon[i] axonal.append()\n',length(PC(i).axon)));
    
    fprintf(fid,'basal = new SectionList()\n');
    fprintf(fid,sprintf('for i=0, %d dend[i] basal.append()\n',length(PC(i).dends)));
    
    fprintf(fid,'apical = new SectionList()\n');
    fprintf(fid,sprintf('for i=0, %d apic[i] apical.append()\n',length(PC(i).apics)));
    
    fprintf(fid,'}\n');
    fclose(fid);
end
%%

% % Replace median volume soma to all the cells:
% for i=1:length(PC)
%     if( i ~= median_volume_soma_idx)
%         PC(i) = PC(i).swapSoma(median_soma);
%         if(~PC(i).isBroken)
%             warning('Error swapping soma in morphology #%d',i);
%         end
%     end
%     fprintf('Done morphology #%d\n',i);
% end


% % Replace the above axon to all morphologies:
% for i=1:length(PC)
%     if( i ~= median_axon_idx)
%         PC(i) = PC(i).swapSubtree(median_axon,'axon');
%         if(~PC(i).isBroken)
%             warning('Error swapping soma in morphology #%d',i);
%         end
%     end
%     fprintf('Done morphology #%d\n',i);
% end