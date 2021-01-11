function [ g ] = rotframe( xi )
%[ g ] = rotframe( xi ): revolute joint
theta=norm(xi(1:3));
omega=xi(1:3)/theta;
v=xi(4:6)/theta;
p=cross(omega,v);
if 0==norm(v)
    g=rotframe1(xi);
else
    g=[v/norm(v),p/norm(p),omega,p;
        0,0,0,1];
end
end

