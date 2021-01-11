function [ xi,xist ] = Puma560Traditional( xi0,xist0,vtheta,gm,M)
%[xi,xist] = Puma560Traditional(xi0,xist0,vtheta,gm,M) : Traditional calibration model for puma560 robot
%
%Input:
%   xi0: nominal joint twists; 
%   xist0: nominal initial twist; 
%   vtheta: joint positions;
%   gm: Actual measured end-effector poses
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

%memory allocation
xi=xi0;
xist=xist0;
N=size(vtheta,1);

gn=zeros(4,4,N);%nominal end-effector poses
dg=zeros(4,4,N);%pose error
vLog=zeros(6,N);%pose error in twist form

for m=1:M
    for i=1:N
        gn(:,:,i)=fkPUMA560(xi,xist,vtheta(i,:),6);
        dg(:,:,i)=gm(:,:,i)/gn(:,:,i);
        vLog(:,i)=vlog(dg(:,:,i));
    end
    simY=zeros(6*N,1);
    for i=1:N
        simY(6*i-5:6*i,1)=vLog(:,i);
    end
    
    %Jacobian matrix
    for k=0:N-1
        simJ(1+6*k:6+6*k,1:6)=aMatrix(xi(:,1),vtheta(k+1,1));
        simJ(1+6*k:6+6*k,7:12)=adM(se3Exp(vtheta(k+1,1)*xi(:,1)))*aMatrix(xi(:,2),vtheta(k+1,2));
        simJ(1+6*k:6+6*k,13:18)=adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2)))*aMatrix(xi(:,3),vtheta(k+1,3));
        simJ(1+6*k:6+6*k,19:24)=adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3)))*aMatrix(xi(:,4),vtheta(k+1,4));
        simJ(1+6*k:6+6*k,25:30)=adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3))*se3Exp(vtheta(k+1,4)*xi(:,4)))*aMatrix(xi(:,5),vtheta(k+1,5));
        simJ(1+6*k:6+6*k,31:36)=adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3))*se3Exp(vtheta(k+1,4)*xi(:,4))*se3Exp(vtheta(k+1,5)*xi(:,5)))*aMatrix(xi(:,6),vtheta(k+1,6));
        simJ(1+6*k:6+6*k,37:42)=adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3))*se3Exp(vtheta(k+1,4)*xi(:,4))*se3Exp(vtheta(k+1,5)*xi(:,5))*se3Exp(vtheta(k+1,6)*xi(:,6)))*aMatrixST(xist);
    end
    
    dp=pinv(simJ)*simY;

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

