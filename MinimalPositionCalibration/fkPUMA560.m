function [ T ] = fkPUMA560( xi,xist,theta,n )
%[ T ] = fkPUMA560( xi,xist,theta,n ) : forward kinematics of PUMA560 type
%robot
%
%Input:
%   xi=6*n, joint twists;
%   xist=6, initial twist
%   theta=n*1, Joint viariables;
%   n, number of joints
%
%Output:
%   T=4*4, homogeneous transformation
%
T=eye(4);
for i=1:n
    T=T*se3Exp(xi(:,i)*theta(i));
end
T=T*se3Exp(xist);
end

