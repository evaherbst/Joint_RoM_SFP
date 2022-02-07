[Q,F,R,angles] = convertKneeData('SS08_FRZN_RH_KNEE',1:6,...
    '../ScanPoints/SS08_RH_knee_femur_plate_POINTS_EVA.txt',...
    '../ScanPoints/SS08_RH_knee_tibfib_plate_POINTS_EVA.txt');
%% This code (written by Enrico Eberhard) has been amended by Eva Herbst to change the axis definitions,
%to enable comparison with rotoscoping data in Maya, where Z = FE, Y = ABAD, and X = LAR (because the Z rotation is
%calculated first, it should be = to FE, the axis with the most movement expected).


%% break the series of rotation matrices into series of X,Y,Z basis vectors



[X,Y,Z] = deal(zeros(size(R,3),3));
for t = 1:size(R,3)
    X(t,:) = R(1,:,t);
    Y(t,:) = R(2,:,t);
    Z(t,:) = R(3,:,t);
    
end



%% The data are now in the right format for the spherical boundaries:
figure(1); clf(1); hold on;
scatter3(X(:,1),X(:,2),X(:,3),'r.'); %changed to blue
scatter3(Y(:,1),Y(:,2),Y(:,3),'g.'); % 
scatter3(Z(:,1),Z(:,2),Z(:,3),'b.'); %changed to red
hold off; axis equal;

%% Fit spherical polygons to the point clusters

% now that the frame axes have been traced onto a sphere, a boundary must
% be fit around each cluster. See the help for boundSphericalPoints for
% more information. In summary, it takes the original point set, an
% optional rotation matrix that helps avoid singularities in the boundary
% mapping, and a shrink factor between 0 and 1. 

X_ = boundSphericalPoints(X',[],0.3);
Y_ = boundSphericalPoints(Y',rotMatY(3*pi/4)*rotMatX(pi/2),0.3);
Z_ = boundSphericalPoints(Z',rotMatY(pi/4)*rotMatX(pi/2),0.3);


%% Plot the spherical frame projection

pts = [X',Y',Z']

fig = makeROMFigure_exvivo_knee(X_,Y_,Z_,pts);

hold on
scatter3(X(:,1),X(:,2),X(:,3),'r.'); %'SizeData',70
scatter3(Y(:,1),Y(:,2),Y(:,3),'g.');
scatter3(Z(:,1),Z(:,2),Z(:,3),'b.');




%% open in vivo figure from rotoscoping to plot ex vivo on same graph
invivo_exvivo_SS06 = openfig("PATCHES_SS06_knee_ventral_with_in_vivo.fig")
%%
%SS06 = openfig("Fig_SS06_knee_1thru6_including_translations_side.fig")
%%

%%

obj = SFPGUI(Q,4);

%%

set(gcf, 'Renderer', 'painters');
%%



%%

%%
plot(obj,invivo_exvivo_SS06,'XYZ',0.99,1.00,1.01)
view(119.547633895177601,60)
%%
