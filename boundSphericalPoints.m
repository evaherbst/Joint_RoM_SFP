function B = boundSphericalPoints(P,R,shrink)
% BOUNDSPHERICALPOINTS
%   Returns boundaries of points on a spherical projection
%
%   B = boundSphericalPoints(P) takes 3xN points P on the surface of a
%       unit sphere and returns the 3xM points B that represent a 
%       spherical boundary.
%
%   First, the spherical points are projected onto a 2D surface using
%   latitude / longitude. Then, the boundary function is used to find
%   points that are on the edge of the cluster. The 2D boundary is 
%   mapped back to the spherical surface, and so the returned points B
%   trace a boundary the spherical points.
%
%   B = boundSphericalPoints(P,R) takes a 3x3 rotation matrix to help with
%       boundary fitting. Latitude and longitude mapping has singularities
%       at the poles (latitude = ±90) and opposite to the prime meridian
%       (longitude = 180); any point clusters near or crossing the
%       singularities will poorly mapped boundaries. R is used to rotate
%       the points onto a better side of the sphere (the mean of the points
%       should be around 0 lat and lon) before boundary fitting, and R'
%       then rotates the boundary back to the original position.
%
%
%   B = boundSphericalPoints(P,R,shrink) or (P,[],shrink) takes an
%       additional input shrink, between 0 and 1, that determines the
%       'tightness' of the boundary mapping. At 0, the boundary is the
%       convex hull of the points, and at 1, the boundary crosses all
%       possible points without self-intersection.
%       Lower shrink values will give 'smoother' shapes but overestimate
%       the area of the cluster. Higher shrink values give jagged edges
%       but better fit the data.
%
%   The best spherical boundary for any point set will be the result of
%   trying various rotations R and shrink values. If the point set is too
%   large to avoid singularities even with rotation, consider splitting
%   along an obvious line (e.g. the equator) and reconstructing the
%   boundaries of the two halves.
%
%   See also: boundary
%
%   2018 Enrico Eberhard

if nargin < 2 || isempty(R)
    R = eye(3);
end
if nargin < 3
    shrink = 0.5;
end


%rotate to central position
P_ = R*P;


%convert to lat,lon
lat = asin(P_(3,:));
lon = atan2(P_(2,:), P_(1,:));

%find 2D boundary indices
j = boundary(lon',lat',shrink);

%reconvert boundary points to spherical points and rotate
B = R'*latlon2cart(lat(j),lon(j));