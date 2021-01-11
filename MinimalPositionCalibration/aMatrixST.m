function [ aM ] = aMatrixST( kesi )
%[ aM ] = aMatrixST( kesi ), matrix A_st in Xiangdong Yang, Liao Wu, Jinquan Li, Ken Chen, A minimal kinematic model for serial robot calibration using POE formula, Robotics and Computer-Integrated Manufacturing, Volume 30, Issue 3, June 2014, Pages 326-334

w=kesi(1:3);
v=kesi(4:6);
bW=zeros(6,6);
bW(1:3,1:3)=crossMatrix(w);
bW(4:6,1:3)=crossMatrix(v);
bW(4:6,4:6)=crossMatrix(w);
t=norm(w);
if t==0
    aM=eye(6);
else
    aM=eye(6)+((4-t*sin(t)-4*cos(t))/2/t^2)*bW+((4*t-5*sin(t)+t*cos(t))/2/t^3)*bW*bW+((2-t*sin(t)-2*cos(t))/2/t^4)*bW*bW*bW+((2*t-3*sin(t)+t*cos(t))/2/t^5)*bW*bW*bW*bW;
end
end

