function [ ad ] = adM( m )
%[ ad ] = adM( m ) : adjoint transformation

r=m(1:3,1:3);
p=m(1:3,4);
ad=zeros(6,6);
ad(1:3,1:3)=r;
ad(4:6,1:3)=crossMatrix(p)*r;
ad(4:6,4:6)=r;
end

