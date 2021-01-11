function [ xi,xist ] = Puma560TraditionalOnlyP( xi0,xist0,vtheta,P01,n1,P02,n2,P03,n3,Pa,M )
%[xi,xist] =
%Puma560TraditionalOnlyP(xi0,xist0,vtheta,P01,n1,P02,n2,P03,n3,Pa,M)
%: Traditional calibration model for puma560 robot using only position
%
%Input:
%   xi0: nominal joint twists; 
%   xist0: nominal initial twist; 
%   vtheta: joint positions;
%   P01: nominal position of first measured point 
%   n1: number of measurements with point 1
%   P02: nominal position of second measured point 
%   n2: number of measurements with point 2
%   P03: nominal position of first measured point 3 
%   n3: number of measurements with point 3
%   Pa=3*(n1+n2+n3): Actual measured position
%   M: iteration steps
%
%Output:
%   xi: calibrated joint twists; 
%   xist: calibrated initial twist; 
%
%Please refer to 
%Xiangdong Yang, Liao Wu, Jinquan Li, Ken Chen, A minimal kinematic model for serial robot calibration using POE formula, Robotics and Computer-Integrated Manufacturing, Volume 30, Issue 3, June 2014, Pages 326-334
%and
%Liao Wu, Xiangdong Yang, Ken Chen, Hongliang Ren. A minimal POE-based model for robotic kinematic calibration with only position measurements. IEEE Transactions on Automation Science and Engineering. 2015, 12(2): 758-763.

xi=xi0;
xist=xist0;
N=size(vtheta,1);
Pn=zeros(3,N);%nominal position of points

for m=1:M
% normdp=inf;
% while (normdp>1e-6)
    for i=1:n1
        temp=fkPUMA560(xi,xist,vtheta(i,:),6)*[P01;1];
        Pn(:,i)=temp(1:3);
    end
    for i=n1+1:n1+n2
        temp=fkPUMA560(xi,xist,vtheta(i,:),6)*[P02;1];
        Pn(:,i)=temp(1:3);
    end
    for i=n1+n2+1:n1+n2+n3
        temp=fkPUMA560(xi,xist,vtheta(i,:),6)*[P03;1];
        Pn(:,i)=temp(1:3);
    end
    
    simY=reshape(Pa-Pn,3*N,1);

    %Jacobian matrix
    for k=0:N-1
        IP=[-crossMatrix(Pn(:,k+1)),eye(3)];
        simJ(1+3*k:3+3*k,1:6)=IP*aMatrix(xi(:,1),vtheta(k+1,1));
        simJ(1+3*k:3+3*k,7:12)=IP*adM(se3Exp(vtheta(k+1,1)*xi(:,1)))*aMatrix(xi(:,2),vtheta(k+1,2));
        simJ(1+3*k:3+3*k,13:18)=IP*adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2)))*aMatrix(xi(:,3),vtheta(k+1,3));
        simJ(1+3*k:3+3*k,19:24)=IP*adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3)))*aMatrix(xi(:,4),vtheta(k+1,4));
        simJ(1+3*k:3+3*k,25:30)=IP*adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3))*se3Exp(vtheta(k+1,4)*xi(:,4)))*aMatrix(xi(:,5),vtheta(k+1,5));
        simJ(1+3*k:3+3*k,31:36)=IP*adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3))*se3Exp(vtheta(k+1,4)*xi(:,4))*se3Exp(vtheta(k+1,5)*xi(:,5)))*aMatrix(xi(:,6),vtheta(k+1,6));
        simJ(1+3*k:3+3*k,37:42)=IP*adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3))*se3Exp(vtheta(k+1,4)*xi(:,4))*se3Exp(vtheta(k+1,5)*xi(:,5))*se3Exp(vtheta(k+1,6)*xi(:,6)))*aMatrixST(xist);
    end
    
    dp=simJ\simY;

    %update
    xi(:,1)=xi(:,1)+dp(1:6);
    xi(1:3,1)=xi(1:3,1)/norm(xi(1:3,1));
    xi(4:6,1)=xi(4:6,1)-xi(1:3,1)'*xi(4:6,1)/(xi(1:3,1)'*xi(1:3,1))*xi(1:3,1);

    xi(:,2)=xi(:,2)+dp(7:12);
    xi(1:3,2)=xi(1:3,2)/norm(xi(1:3,2));
    xi(4:6,2)=xi(4:6,2)-xi(1:3,2)'*xi(4:6,2)/(xi(1:3,2)'*xi(1:3,2))*xi(1:3,2);

    xi(:,3)=xi(:,3)+dp(13:18);
    xi(1:3,3)=xi(1:3,3)/norm(xi(1:3,3));
    xi(4:6,3)=xi(4:6,3)-xi(1:3,3)'*xi(4:6,3)/(xi(1:3,3)'*xi(1:3,3))*xi(1:3,3);

    xi(:,4)=xi(:,4)+dp(19:24);
    xi(1:3,4)=xi(1:3,4)/norm(xi(1:3,4));
    xi(4:6,4)=xi(4:6,4)-xi(1:3,4)'*xi(4:6,4)/(xi(1:3,4)'*xi(1:3,4))*xi(1:3,4);

    xi(:,5)=xi(:,5)+dp(25:30);
    xi(1:3,5)=xi(1:3,5)/norm(xi(1:3,5));
    xi(4:6,5)=xi(4:6,5)-xi(1:3,5)'*xi(4:6,5)/(xi(1:3,5)'*xi(1:3,5))*xi(1:3,5);

    xi(:,6)=xi(:,6)+dp(31:36);
    xi(1:3,6)=xi(1:3,6)/norm(xi(1:3,6));
    xi(4:6,6)=xi(4:6,6)-xi(1:3,6)'*xi(4:6,6)/(xi(1:3,6)'*xi(1:3,6))*xi(1:3,6);

    xist=xist+dp(37:42);

end

end

