function [ w ] = vlogR( R )
%[ w ] = vlogR( R ) : logarithm mapping of rotation
tr=(trace(R)-1)/2;
if tr>=1
    tr=1;
elseif tr<=-1
    tr=-1;
end
fai=acos(tr);

if fai==0
    w=[0;0;0];
else
    w=fai/(2*sin(fai))*[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
end

end

