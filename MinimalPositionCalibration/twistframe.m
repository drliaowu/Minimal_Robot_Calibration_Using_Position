function [ g ] = twistframe( xi )
%[ g ] = twistframe( xi ) gives the transformation (g) from base frame to
%link frame
theta=norm(xi(1:3));
if 0==theta
    g=transframe(xi);
else
    g=rotframe(xi);
end
end

