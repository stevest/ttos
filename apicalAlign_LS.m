
% function AXIS = apicalAlign_LS(idx)
function[rettree] = apicalAlign_LS(idx, trees)

% global trees

%apical subset:
apc = find(4 == cellfun(@str2num,trees{idx}.rnames(trees{idx}.R)) );
wgt = trees{idx}.D(apc);
wgt = wgt/sum(wgt);


% weighted 3d regression:
r0=mean([trees{idx}.X(apc).*wgt,trees{idx}.Y(apc).*wgt,trees{idx}.Z(apc).*wgt]);
xyz=bsxfun(@minus,[trees{idx}.X(apc).*wgt,trees{idx}.Y(apc).*wgt,trees{idx}.Z(apc).*wgt],r0);
[~,~,V]=svd(xyz,0);

%find angle and axis to align the X axis
angle = rad2deg(GetAbsAngleBetweenVectors(V(1:3,1), [0,-1,0] ) );% in rads
axis = cross(V(1:3,1)/norm(V(1:3,1)),[0,-1,0] );
axis = axis/norm(axis);
% Rotate to align PCA first axis with Z axis:
f = angle;
ux = axis(1);
uy = axis(2);
uz = axis(3);
R = [cosd(f)+(ux^2)*(1-cosd(f)) , ux*uy*(1-cosd(f))-uz*sind(f) , ux*uz*(1-cosd(f))+uy*sind(f);...
    uy*ux*(1-cosd(f))+uz*sind(f) , cosd(f) + (uy^2)*(1-cosd(f)) , uy*uz*(1-cosd(f))-ux*sind(f);...
    uz*ux*(1-cosd(f))-uy*sind(f) , uz*uy*(1-cosd(f))+ux*sind(f) , cosd(f)+(uz^2)*(1-cosd(f)) ];
T = [trees{idx}.X(1);trees{idx}.Y(1);trees{idx}.Z(1)];
tempVar = V(1:3,1)/norm(V(1:3,1))
AXIS = [tempVar;T];

% plot_tree(trees{idx}, cellfun(@str2num,trees{idx}.rnames(trees{idx}.R))' );
% hold on;
 
% for i=1:length(trees{idx}.X)
%     tempPoint = R * ([trees{idx}.X(i);trees{idx}.Y(i);trees{idx}.Z(i)]-T);
%     trees{idx}.X(i) = tempPoint(1);
%     trees{idx}.Y(i) = tempPoint(2);
%     trees{idx}.Z(i) = tempPoint(3);
% end

onlyApical = find(4 == cellfun(@str2num,trees{idx}.rnames(trees{idx}.R)) );
for i=onlyApical
    tempPoint = R * ([trees{idx}.X(i);trees{idx}.Y(i);trees{idx}.Z(i)]-T);
    trees{idx}.X(i) = tempPoint(1);
    trees{idx}.Y(i) = tempPoint(2);
    trees{idx}.Z(i) = tempPoint(3);
end

% % Calculate the direction of the apical and reverse the cell if necessary!
% if ( mean(trees{idx}.Z(apc)) < 0 )
%     % Rotate 180 in x axis
%     R = [1,0,0;
%         0 , cosd(180)  , -sind(180);...
%         0 , sind(180) , cosd(180) ];
%     for i=1:length(trees{idx}.X)
%         tempPoint = R * [trees{idx}.X(i);trees{idx}.Y(i);trees{idx}.Z(i)];
%         trees{idx}.X(i) = tempPoint(1);
%         trees{idx}.Y(i) = tempPoint(2);
%         trees{idx}.Z(i) = tempPoint(3);
%     end
%     AXIS = [-tempVar;T];
% end

% % % rotate in Y for !@#$@!#!@#$ NEURON:
% % % Rotate 180 in x axis
% %     Ry = [cosd(90),0,sind(90);...
% %         0 , 1  ,0;...
% %         -sind(90) , 0 , cosd(90) ];
% %     for i=1:length(trees{idx}.X)
% %         tempPoint = Ry * [trees{idx}.X(i);trees{idx}.Y(i);trees{idx}.Z(i)];
% %         trees{idx}.X(i) = tempPoint(1);
% %         trees{idx}.Y(i) = tempPoint(2);
% %         trees{idx}.Z(i) = tempPoint(3);
% %     end

 

%  plot3( [AXIS(4), AXIS(1)*100], [AXIS(5), AXIS(2)*100], [AXIS(6), AXIS(3)*100] );
% pause();
% clf;

% return;
rettree = trees{idx};
end