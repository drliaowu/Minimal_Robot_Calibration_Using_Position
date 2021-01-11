function [ kesi ] = vlog( g )
%[ kesi ] = vlog( g ) : logarithm mapping from SE(3) to R6
fai=rotationTheta(g);
w=fai*rotationW(g,fai);
if fai==0
    p=g(1:3,4);
else
    p=(eye(3)-0.5*crossMatrix(w)+(2*sin(fai)-fai*(1+cos(fai)))/(2*fai*fai*sin(fai))*crossMatrix(w)*crossMatrix(w))*g(1:3,4);
end
kesi=[w;p];

end

