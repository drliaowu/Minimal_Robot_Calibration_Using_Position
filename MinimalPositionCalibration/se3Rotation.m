function [ T ] = se3Rotation( w,v,theta )
% [ T ] = se3Rotation( w,v,theta ) generates a homogeneous tranformation
% from a motion consisting both rotation and translation

R=rotationMatrix(w,theta);
p=(theta*eye(3)+(1-cos(theta))*crossMatrix(w)+(theta-sin(theta))*crossMatrix(w)*crossMatrix(w))*v;
T=[R,p;
    0,0,0,1];


end

