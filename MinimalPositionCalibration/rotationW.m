function [ w ] = rotationW( g,theta )
%[ w ] = rotationW( g,theta ) gives rotation axis of g
if theta==0
    w=[0;0;0];
else
    w=1/(2*sin(theta))*[g(3,2)-g(2,3);g(1,3)-g(3,1);g(2,1)-g(1,2)];
end

