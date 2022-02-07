function fit = plateTransform(points)
%plateTransform
%
%   fit = plateTransform(points) generates a plate coordinate frame
%       from a structure of XYZ points on specific layers. First
%       get points using the function parsePointsFromRhino.
%
%   fit is a struct containing best fit lines, planes and cylinders
%       for the corresponding layers in points, as well as an Origin
%       point and normalized orthogonal X, Y, Z directional axes.
%
%   fit also contains a transformation matrix TM that transforms points
%       from scan frame into the local plate frame, TMmarker that
%       transforms from scan frame into marker frame, and TMproxworld
%       that transforms from scan frame into world frame (for the static
%       proximal plate).
%
%   fit.Points contains the original input points, fit.PointsMarker is the
%       input points transformed by TMmarker, and fit.PointsProx is the
%       input points transformed by TMproxworld
%
%   See also parsePointsFromRhino.
%
%   2018 Enrico Eberhard





% check what data is available to fit
layers = fieldnames(points);


fit = struct;


%fit some lines
for ind = find(contains(layers,'Axis'))'
    
    xyz = points.(layers{ind});
    
    fit.(layers{ind}) = fitLine(xyz);
    
end



%fit some planes
for ind = find(contains(layers,'Plane'))'
    
    xyz = points.(layers{ind});
    
    fit.(layers{ind}) = fitPlane(xyz);
    
end

%get/update axes from plane intersections
fit = axesFromPlanes(fit);




%fit some cylinders
for ind = find(contains(layers,'Cyl'))'
    
    xyz = points.(layers{ind});
    
    [C, W, ~] = fitCylinder(xyz);
    
    fit.(layers{ind}).u = C;
    fit.(layers{ind}).V = W;
    
end
    

%find or update axes from cylinder fit
fit = axesFromCyls(fit);

%find the corner 
fit.Origin = cornerPoint(fit);


%orthoganalise axes
[X, Y, Z] = orthogonalize(fit.XAxis.V, fit.YAxis.V, fit.ZAxis.V);

fit.X = X;
fit.Y = Y;
fit.Z = Z;


R = ...
[dot([1 0 0], X), dot([1 0 0], Y), dot([1 0 0], Z);
 dot([0 1 0], X), dot([0 1 0], Y), dot([0 1 0], Z);
 dot([0 0 1], X), dot([0 0 1], Y), dot([0 0 1], Z)];

%Scan frame to plate frame
fit.TM = [[R zeros(3,1)]; [-fit.Origin*R 1]];


%(distal)
%Scan frame to marker frame
TMmarker = ...
[        0    1.0000         0         0;
    1.0000         0         0         0;
         0         0   -1.0000         0;
  -30.0000   -1.0500  -13.1500    1.0000];

fit.TMmarker = fit.TM*TMmarker;



%(proximal)
%Scan frame to load cell frame
TMloadcell = ...
[        0    1.0000         0         0;
   -1.0000         0         0         0;
         0         0    1.0000         0;
   30.0000   -1.0500   10.3000    1.0000];

fit.TMloadcell = fit.TM*TMloadcell;


%Scan plate frame to world frame
TMproxworld = ...
[        0         0    1.0000         0;
   -1.0000         0         0         0;
         0   -1.0000         0         0;
   71.0600  -20.9600  103.2500    1.0000];

calib = [3.7050   -2.4708    2.6850];
TMproxworld(4,1:3) = TMproxworld(4,1:3) - calib;

fit.TMproxworld = fit.TM*TMproxworld;



fit.Points = points.All;

pt = [points.All ones(length(points.All),1)]*fit.TM;
fit.PointsPlate = pt(:,1:3);

pt = [points.All ones(length(points.All),1)]*fit.TMmarker;
fit.PointsMarker = pt(:,1:3);

pt = [points.All ones(length(points.All),1)]*fit.TMproxworld;
fit.PointsProx = pt(:,1:3);


end





function L = fitLine(xyz)

    %center points around mean
    u=mean(xyz);
    xyz0=bsxfun(@minus,xyz,u);

    %decompose
    [~,~,V]=svd(xyz0,0);

    L.u = u; L.V = V(:,1)';

end


function P = fitPlane(xyz)

    u = mean(xyz);

    %move points so mean is at (0,0,0)
    R = bsxfun(@minus,xyz,u);

    %linear algebra bit: make the point matrix
    %symmetric 3x3 by multiplying by the transpose
    %then compute principal directions of eigenvectors
    [V,~] = eig(R'*R);

    %first eigenvector is normal to the plane
    P.u = u; P.V = V(:,1)';


end




function fit = axesFromPlanes(fit)
    
    lays = fieldnames(fit);
    planes = {lays{contains(lays,'Plane')}};
    
    %first update each axis direction that is normal to the plane
    
    
    
    
    switch(numel(planes))
        case 2
            u = [fit.(planes{1}).u; fit.(planes{2}).u];
            V = [fit.(planes{1}).V; fit.(planes{2}).V];
            
        case 3
            u = [fit.(planes{1}).u; fit.(planes{2}).u; fit.(planes{3}).u];
            V = [fit.(planes{1}).V; fit.(planes{2}).V; fit.(planes{3}).V];
            
        otherwise
            warning('Not enough planes to find any axes!')
            return
    end
    
    pairs = nchoosek(1:numel(planes),2);
            
    for pair = pairs'

        [ui, Vi] = intersectionOfPlanes(u(pair,:),V(pair,:));


        if sum(contains(planes(pair), 'X')) == 2
            ax = 'XAxis';
        elseif sum(contains(planes(pair), 'Y')) == 2
            ax = 'YAxis';
        elseif sum(contains(planes(pair), 'Z')) == 2
            ax = 'ZAxis';
        else
            error('No plane pairs found matching expected names');
        end


        %check if this fit already exists;

        if isfield(fit, ax)
        %   if so, average it (in same direction)
            if dot(fit.(ax).V, Vi) < 0
                Vi = -Vi;
            end
            fit.(ax).V = normalize((fit.(ax).V + Vi)/2);
            fit.(ax).u = (fit.(ax).u + ui)/2;
        else
            fit.(ax).V = Vi;
            fit.(ax).u = ui;
        end

    end
    
    
    
    
end



function fit = axesFromCyls(fit)

    fits = fieldnames(fit);
    cyls = sort({fits{contains(fits,'Cyl')}});
    
    if ~numel(cyls) %no cylinders!
        return
    end
    
    %direction of cylinder is X axis
    %so, first average all cylinder W
    
    V_ = zeros(1,3);
    for cyl = cyls
        
        if dot(V_, fit.(cyl{1}).V) < 0
            V_ = -V_;
        end
        V_ = V_ + fit.(cyl{1}).V;
        
    end
    
    if isfield(fit, 'XAxis')
        if dot(fit.XAxis.V, V_) < 0
            V_ = -V_;
        end
        fit.XAxis.V = normalize((fit.XAxis.V + V_)/(1 + numel(cyls)));
    else
        fit.XAxis.V = normalize(V_);
    end
    
    
    
    %if there are multiple cyls, Z axis direction runs through
    %centers and must be perpendicular to X axis
    if numel(cyls) > 1
        u = zeros(0,3);
        for cyl = cyls
            u = [u; fit.(cyl{1}).u];
        end
        L = fitLine(u);

        %project fit line onto YZPlane
        Y_ = cross(fit.XAxis.V,L.V);
        Z_ = cross(fit.XAxis.V,Y_);

        %correct direction of Z from order of cyls
        if dot(Z_, u(end,:) - u(1,:)) < 0
            Z_ = -Z_;
        end

        %make or update ZAxis
        if isfield(fit, 'ZAxis')
            if dot(fit.ZAxis.V, Z_) < 0
                fit.ZAxis.V = -fit.ZAxis.V;
            end
            fit.ZAxis.V = normalize(...
                    (fit.ZAxis.V + Z_*(numel(cyls)-1))/numel(cyls));
        else
            fit.ZAxis.V = normalize(Z_);
        end
        
        
        Y_ = cross(fit.ZAxis.V, fit.XAxis.V);
        
        %make or update YAxis
        if isfield(fit, 'YAxis')
            if dot(fit.YAxis.V, Y_) < 0
                fit.YAxis.V = -fit.YAxis.V;
            end
            fit.YAxis.V = normalize(...
                    (fit.YAxis.V + Y_*(numel(cyls)-1))/numel(cyls));
        else
            fit.YAxis.V = normalize(Y_);
        end
    end
    
    
    
    %point on x axis is -2.5 Z, 2.5 Y from Cyl1
    %need to know: 1) correct direction of axes
    %              2) cyl centroid projected on YZPlane
    
    %positive XAxis is from plate interior towards YZPlane
    if isfield(fit, 'YZPlane')
        if dot(fit.XAxis.V,normalize(fit.YZPlane.u - fit.(cyls{1}).u)) < 0
            fit.XAxis.V = -fit.XAxis.V;
        end
    end
    
    %similarly update Z direction
    if isfield(fit, 'XYPlane')
        if ~isfield(fit,'ZAxis')
            fit.ZAxis.V = fit.XYPlane.V;
        end
        if dot(fit.ZAxis.V,normalize(fit.(cyls{1}).u - fit.XYPlane.u)) < 0
            fit.ZAxis.V = -fit.ZAxis.V;
        end
    end
    
    %similarly update Y direction
    if isfield(fit, 'ZXPlane')
        if ~isfield(fit,'YAxis')
            fit.YAxis.V = fit.ZXPlane.V;
        end
        if dot(fit.YAxis.V,normalize(fit.ZXPlane.u - fit.(cyls{1}).u)) < 0
            fit.YAxis.V = -fit.YAxis.V;
        end
    end
    
    %there should now be enough axes in the correct direction to 
    %determine the rest
    
    fits = fieldnames(fit);
    axes = sort({fits{contains(fits,'Axis')}});
    
    if numel(axes) == 2
        if contains(axes{1}, 'YAxis')
            fit.XAxis.V = cross(fit.YAxis,V, fit.ZAxis.V);
        elseif contains(axes{2}, 'ZAxis')
            fit.YAxis.V = cross(fit.ZAxis,V, fit.XAxis.V);
        else
            fit.ZAxis.V = cross(fit.XAxis,V, fit.YAxis.V);
        end
    elseif numel(axes) < 2
        error('Not enough point data to construct all axes');
    end
        
    
    if ~isfield(fit.YAxis, 'u') && ...
       ~isfield(fit.ZAxis, 'u') && ...
       ~isfield(fit, 'YZPlane')
   
        error(['Not enough information in X axis direction. '...
               'Need at least 1 of YAxis, ZAxis or YZPlane']);
    end
    
    %center points can be used to find corner point
    % and points on axes
    
    Xu = 0;
    for cyl = cyls
        
        u = fit.(cyl{1}).u;
        
        cyln = str2double(cyl{1}(4)) - 1;
        
        xu = u + 2.5*fit.YAxis.V - (cyln*5 + 2.5)*fit.ZAxis.V;
        
        Xu = Xu + xu;
        
    end
    
    if isfield(fit.XAxis, 'u')
        fit.XAxis.u = (fit.XAxis.u + Xu) / (numel(cyls) + 1);
    else
        fit.XAxis.u = Xu / numel(cyls);
    end
    
    

    %find intersection point of XAxis and YZPlane,
    %i.e. closest point on XAxis to either Y or Z Axis
    
    if isfield(fit.YAxis, 'u') && ~isfield(fit.ZAxis, 'u')
        fit.ZAxis.u = intersectionOfLines([fit.XAxis.u; fit.YAxis.u],...
                                        [fit.XAxis.V; fit.YAxis.V]);
              
    elseif ~isfield(fit.YAxis, 'u')
        fit.YAxis.u = intersectionOfLines([fit.XAxis.u; fit.ZAxis.u],...
                                        [fit.XAxis.V; fit.ZAxis.V]);
                                    
    end
        
end


function p = cornerPoint(fit)

    fits = fieldnames(fit);
    axes = {fits{contains(fits,'Axis')}};

    p = [];
    if numel(axes) < 2
        error('Not enough axes to determine corner point!')
    else

        u = []; V = [];
        for axis = axes
            u = [u; fit.(axis{1}).u];
            V = [V; fit.(axis{1}).V];
        end
        
        p = intersectionOfLines(u,V);
            
    end
end


