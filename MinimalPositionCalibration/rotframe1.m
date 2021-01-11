function [ g ] = rotframe1( xi )
%[ g ] = rotframe1( xi ) : joint axis passes through origin
omega=xi(1:3)/norm(xi(1:3));
axis=cross([0;0;1],omega);
phi=acos([0,0,1]*omega);
switch phi
    case 0
        g=eye(4);
    case pi
        g=[1,0,0,0;
            0,-1,0,0;
            0,0,-1,0;
            0,0,0,1];
    otherwise
        g=se3Rotation(axis/norm(axis),[0;0;0],phi);
end

end

