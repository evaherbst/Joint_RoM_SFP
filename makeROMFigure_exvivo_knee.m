function [fig,pl] = makeROMFigure_exvivo_knee(X,Y,Z,pts,fig)
% MAKEROMFIGURE
%   Generates a spherical frame projection
%
%   makeROMFigure(X,Y,Z) takes 3xN points X, Y and Z that each respectively
%   trace a spherical polygon, and plots them in red, green and blue on a
%   background unit sphere.
%
%   makeROMFigure(...,pts) additionally plots any 3xM points (for instance,
%   the raw point clusters that were used to generate the original
%   boundaries) as small black pts
%
%   makeROMFigure(...,pts,fig) or (...,[],fig) specifies a figure handle
%   for plotting.
%
%   fig = makeROMFigure(...) returns the figure handle
%
%   [fig, pl] = makeRomFigure(...) also returns the handle to the plotted
%   boundaries (for adjusting colors/line thickness)
%
%   See also: boundSphericalPoints

[v,f] = spheretri(600);
[sph.X,sph.Y,sph.Z] = sphere;

if nargin < 4
    pts = [];
end

if nargin < 5
    fig = figure();
end

figure(fig); clf(fig); hold on;

pl.X = plot3(X(1,:),X(2,:),X(3,:));
pl.Y = plot3(Y(1,:),Y(2,:),Y(3,:));
pl.Z = plot3(Z(1,:),Z(2,:),Z(3,:));

if nargin > 3 && ~isempty(pts)
	%pl.sc = scatter3(pts(1,:),pts(2,:),pts(3,:));
   % pl.sc.Marker = '.'; pl.sc.SizeData = 10; pl.sc.CData = [0 0 0];
 


end

sph.p = patch('Vertices',v,'Faces',f,'FaceColor','k','EdgeColor','k');
sph.srf = surf(1*sph.X,1*sph.Y,1*sph.Z,'FaceColor',[0.9 0.9 0.9],'EdgeColor','k');



%note: this axis recoloring accounts for the polygon rotation 
%all rotated by (-pi/2) in ExampleROMScriptRefactored
%This visual rotation is used for the knee, where you want Y up (at the
%pole)

plot3([0 1], [0 0], [0 0], 'r', 'LineWidth', 2); 

plot3([0 0], [0 1], [0 0], 'g', 'LineWidth', 2); 

plot3([0 0], [0 0], [0 1], 'b', 'LineWidth', 2); 

hold off;


sph.p.FaceAlpha = 0.8; sph.p.EdgeAlpha = 0.1; sph.p.Visible = 'off';
sph.srf.FaceAlpha = 0.4; sph.srf.EdgeAlpha = 0.1;

%boundary polygon colors and line widths

pl.X.LineWidth = 1.5; pl.X.Color = [0.8 0.2 0.2]; %orig size 2.5
pl.Y.LineWidth = 1.5; pl.Y.Color = [0.2 0.8 0.2];
pl.Z.LineWidth = 1.5; pl.Z.Color = [0.2 0.2 0.8]; 
 

axis equal; axis([-1 1 -1 1 -1 1]);view(0,180);%view(119.547633895177601,60) for side view, %view(0,180) for ventral for ventral view
xlabel('X'); ylabel('Y'); zlabel('Z'); %changelabels
camup([0 1 0])


end