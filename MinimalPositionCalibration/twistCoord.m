function [ kesi ] = twistCoord( w, p )
%[kesi]=twistCoord(w,p) generates the twist coordinate, xi, from axis direction, w, and a point along axis, p. 

if norm(w)~=0
    kesi=[w;cross(p,w)];
else
    kesi=[0;0;0;p];
end

end

