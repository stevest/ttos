% Add needed toolboxes to the path:
% NOTE: minor bugs changed in TREESTOOLBOX, so use the provided one.
addpath('~/Documents/MATLAB/Dendrites/');
addpath(genpath('~/Documents/MATLAB/Dendrites/TREES1.15'));

% Load morphologies names from Smith's lab
names = dir('SmithMorphologies');
names = {names(~[names.isdir]).name};

% We are using A LOT of recursions (as much as the length of points of
% swc files), so set limit accordingly to avoid warnings:
set(0,'RecursionLimit',1000);

% Load all morphologies as ncell objects:
for i=1:length(names)
    PC(i) =  ncell([],1,['SmithMorphologies/',names{i}]);
end

% Get soma volume of each ncell:
for i=1:length(PC)
    volume_soma(i) = mean(cellfun(@(x) mean(x),ncell.regionVolume(PC(i).soma)));
end

% Poorly implemented median:
[~,idx ]= sort(volume_soma);
tmp = volume_soma(idx);
tmp(round(length(tmp)/2));
median_volume_soma_idx = find(tmp(round(length(tmp)/2))==volume_soma);

median_soma = PC(median_volume_soma_idx).getSubtree('soma');

% find the simplest/more complex apical tree (based in bifurcation order):
for i=1:length(PC)
    BO(i) = max(BO_tree(PC(i).tree));
end
[~,simple_apical_idx] = min(BO);
[~,complex_apical_idx] = max(BO);

% Find a morphology with a reconstructed axon:
for i=1:length(PC)
    if( ~isempty(PC(i).axon) )
        median_axon = PC(i).getSubtree('axon');
        median_axon_idx = i;
        break;
    end
end


% Save simple apical morphologies as .hoc files (NEURON)
simple_dir = 'smith_simple_apical';
system(['mkdir ',simple_dir]);
simple_apical = PC(simple_apical_idx).getSubtree('apic');

for i=1:length(PC)
    % Swap apical:
    tmpCell = PC(i).swapSubtree(simple_apical,'apic');
    % Swap axon:
    tmpCell = tmpCell.swapSubtree(median_axon,'axon');
    % Swap soma:
    tmpCell = tmpCell.swapSoma(median_soma);
    % Save as .hoc file using custom function
    neuron_tree_mod(tmpCell.tree,[simple_dir,'/',tmpCell.tree.name,'.hoc'],1,'');
%     swc_tree(tmpCell.tree,[simple_dir,'/',tmpCell.tree.name,'.swc']);
end

% Save complex apical morphologies as .hoc files (NEURON)
complex_dir = 'smith_complex_apical';
system(['mkdir ',complex_dir]);
complex_apical = PC(complex_apical_idx).getSubtree('apic');

for i=1:length(PC)
    % Swap apical:
    tmpCell = PC(i).swapSubtree(complex_apical,'apic');
    % Swap axon:
    tmpCell = tmpCell.swapSubtree(median_axon,'axon');
    % Swap soma:
    tmpCell = tmpCell.swapSoma(median_soma);
    % Save as .hoc file using custom function
    neuron_tree_mod(tmpCell.tree,[complex_dir,'/',tmpCell.tree.name,'.hoc'],1,'');
%     swc_tree(tmpCell.tree,[complex_dir,'/',tmpCell.tree.name,'.swc']);
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