function [TM, CS, origin, PointonX, PointonY, PointonZ] = knee_ACS(points, bone)

%Eva Herbst, July 2019 (built on Enrico Eberhard's knee SSCS code for ex vivo
%experiments, but with some axes calculated differently, and with directions
%and names changed)

%This function calculates ACSs for the knee joint (distal femur ACS and
%proximal tibfib ACS).
%These axes are stored in the CS structure (so for X axis unit vector, see
%CS.X)

%Here, Z = FE, Y = ABAD, and X = LAR (following XROMM/rotoscoping
%conventions, based on calculation order in Maya)

%For distal femur, long axis (X) is recalculated to make axes orthogonal, for
%proximal tibia/fibula, FE axis (Y) is recalculated to make axes
%orthogonal.

dist = points.BoneDist;
prox = points.BoneProx;
points.Knee_flex_side = points.YSide  %add or comment out, depending on if pts have YSide or Knee_flex_side..
knee_flex_side = points.Knee_flex_side; 


        
switch(bone)
    case 'femur'
        
        Z = fitLine(dist); 
        
        u = mean(prox);
        
        X.V = normalize(u - Z.u)
        
        Y.V = cross(Z.V, X.V);
        
    case 'tibfib' %add new code that recalculates FE (Z) instead of LAR (X)
        
        Z = fitLine(prox);
        
        u = mean(dist);
        
        X.V = normalize(u - Z.u);
        
        Y.V = cross(X.V, Z.V); 
        
        
        
    otherwise
        warning('Please specify bone as ''femur'' or ''tibfib''.')
        
end


Y.u = Z.u;
X.u = Z.u;
origin = Z.u



%check alignment of Z/Y vectors




switch(bone)

case 'femur' 
    
   if dot(normalize(knee_flex_side - Z.u), Y.V) < 0
   
    Y.V = -Y.V;
    Z.V = - Z.V
    

      
    end
  
    
    if dot(normalize(knee_flex_side - Z.u), Y.V) > 0
       
     Z.V = -Z.V;
    end
     X.V = cross(Y.V, Z.V);
    


case 'tibfib' %update FE (Z)    
    
   if dot(normalize(knee_flex_side - Z.u), Y.V) < 0
   
    Y.V = -Y.V;
    Z.V = - Z.V
    

    end

     
Z.V = cross(X.V, Y.V) 

end
    


CS.X = normalize(X.V);
CS.Y = normalize(Y.V);
CS.Z = normalize(Z.V);

PointonX = origin + CS.X
PointonY = origin + CS.Y
PointonZ = origin + CS.Z

R = ...
[dot([1 0 0], CS.X), dot([1 0 0], CS.Y), dot([1 0 0], CS.Z);
 dot([0 1 0], CS.X), dot([0 1 0], CS.Y), dot([0 1 0], CS.Z);
 dot([0 0 1], CS.X), dot([0 0 1], CS.Y), dot([0 0 1], CS.Z)];

TM = [[R zeros(3,1)]; [-X.u*R 1]];


end


function L = fitLine(xyz)

    %center points around mean
    u=mean(xyz);
    xyz0=bsxfun(@minus,xyz,u);

    %decompose
    [~,~,V]=svd(xyz0,0);

    L.u = u; L.V = V(:,1)';

end
