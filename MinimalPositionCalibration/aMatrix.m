function [ aM ] = aMatrix( kesi, q )
%[ aM ] = aMatrix( kesi, q ) : A matrix in (13) and (14)

w=kesi(1:3);
v=kesi(4:6);
bW=zeros(6,6);
bW(1:3,1:3)=crossMatrix(w);
bW(4:6,1:3)=crossMatrix(v);
bW(4:6,4:6)=crossMatrix(w);
n=norm(w);
t=n*q;
if n==0
    aM=q*eye(6);
else
    aM=q*eye(6)+((4-t*sin(t)-4*cos(t))/2/n^2)*bW+((4*t-5*sin(t)+t*cos(t))/2/n^3)*bW*bW+((2-t*sin(t)-2*cos(t))/2/n^4)*bW*bW*bW+((2*t-3*sin(t)+t*cos(t))/2/n^5)*bW*bW*bW*bW;
end

end

