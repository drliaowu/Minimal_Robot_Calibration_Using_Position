function [ theta ] = rotationTheta( g )
% [ theta ] = rotationTheta( g ) : gets the rotation angle of a homogenous
% transformation
tr=(trace(g(1:3,1:3))-1)/2;
if tr>=1
    tr=1;
elseif tr<=-1
    tr=-1;
end
theta=acos(tr);

end

