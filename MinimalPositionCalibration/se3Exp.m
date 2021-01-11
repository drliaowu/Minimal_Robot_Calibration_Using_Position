function [ T ] = se3Exp( kesi )
%[ T ] = se3Exp( kesi ) generates homogeneous transformation (T) from a twist (kesi).

n1=norm(kesi(1:3));
n2=norm(kesi(4:6));

if n1==0&&n2==0
    T=eye(4);
elseif n1==0
    T=se3Translation(kesi(4:6)/n2,n2);
else
    T=se3Rotation(kesi(1:3)/n1,kesi(4:6)/n1,n1);
end

end

