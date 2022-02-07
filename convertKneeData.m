function [Q,R] = convertKneeData(specimen,subsets,rhinoFileProx,rhinoFileDist)

[p,q,f] = parseROMTrialData(specimen, subsets);


% trim points that translate too much
[c,r] = sphereFit(p);

%spherical error
re = sqrt(sum((p-c).^2,2)) - r;

%threshold
i = (re.^2)< 100; %set high threshold to include all translations. 
%Can also set lower cutoff to exclude some translations. Not spherical
%error is used here, but could also write new function to measure
%translations from joint center for future studies

fprintf('Removing %i highly translated data points from %i original\n',...
        length(p) - sum(i), length(p));

q = q(i,:);
f = f(i,:);

femur_points = parsePointsFromRhino(rhinoFileProx);
femur_fit = plateTransform(femur_points);

tibfib_points = parsePointsFromRhino(rhinoFileDist);
tibfib_fit = plateTransform(tibfib_points);

femur_bone = knee_ACS(femur_points, 'femur'); 
tibfib_bone = knee_ACS(tibfib_points, 'tibfib'); 

femur_in_world = (femur_bone \ femur_fit.TM) * (femur_fit.TM \ femur_fit.TMproxworld);
tibfib_in_marker = (tibfib_bone \ tibfib_fit.TM) * (tibfib_fit.TM \ tibfib_fit.TMmarker);

% get rotation matrices of the marker
R = zeros(3,3,length(q));
Q = q;
for t = 1:length(q)
    markerTM = quatToTM(q(t,:));
    tibfib_in_world = tibfib_in_marker * markerTM;

    R(:,:,t) = rotFromAxes(femur_in_world(1:3,1:3), tibfib_in_world(1:3,1:3));
    Q(t,:) = quatFromRot(R(:,:,t));
    
    
     

 
end



end