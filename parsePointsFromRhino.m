function points = parsePointsFromRhino(rhinofile)
%parsePointsFromRhino
%
%   points = parsePointsFromRhino(rhinoFile) takes a text file generated
%       by the Rhino "what" command (known as object property details),
%       and returns a points struct, where each field is a layer name with
%       Nx3 XYZ points.
%
%   2018 Enrico Eberhard


%%Points are unnamed and in one of the following layers:
%
%   Axes:
%       XAxis
%       YAxis
%       ZAxis
%
%   Planes:
%       XYPlane
%       YZPlane
%       ZXPlane
%
%   Cyls:
%       Cyl1
%       Cyl2
%       Cyl3

% Bone Prox
% Bone Dist

% Example entry:

% point  
%   ID: 3af5750b-b4ef-46e9-ae0c-00949a962d53 (212)
%   Layer name: XAxis
%   Render Material: 
%     source = from layer
%     index = -1
%   Geometry:
%     Valid point.
%     Point at (-2.677,0.407,-5.082).




text = fileread(rhinofile);


%%What to look for if SITES
%find "Layer name: " and take all non-whitespace characters following
%skip the minimum number of characters until 
% "Point at (" XVAL "," YVAL "," ZVAL ")."

layerExp = 'point.*?Layer name: (?<name>.*?)\s';
xExp = 'Point at \x{28}(?<x>.*?),';
yExp = '(?<y>.*?),';
zExp = '(?<z>.*?)\x{29}';

%concatenate the expressions
exp = [layerExp,'.*?',xExp,yExp,zExp];


%find the point layers and coordinates using named tokens
pointData = regexp(text,exp,'names');

%put names and x y z coords into nx4 cell array
pointCellArray = squeeze(struct2cell(pointData));

%separate layer and coordinate data
layers = pointCellArray(1,:);
xyz = str2double(pointCellArray(2:4,:));


%parse into point clouds per layer
points = struct;

for n = 1:length(layers)
    
    if ~isfield(points, layers{n})
        points.(layers{n}) = xyz(:,n)';
    else
        points.(layers{n}) = [points.(layers{n}); xyz(:,n)'];
    end
    
end



layers = fieldnames(points);

points.All = zeros(0,3);
for layer = layers'
    points.All = [points.All; points.(layer{1})];
end


end