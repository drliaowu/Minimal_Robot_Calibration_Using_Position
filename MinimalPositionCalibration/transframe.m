function [ g ] = transframe( xi )
% [ g ] = transframe( xi ) : prismatic joint
v=xi(4:6)/norm(xi(4:6));
axis=cross([0;0;1],v);
phi=acos([0,0,1]*v);
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

