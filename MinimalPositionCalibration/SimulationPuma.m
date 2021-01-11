%%
% Please refer to Section III, SIMULATIONS ON A PUMA 560 ROBOT 
% Liao Wu, Xiangdong Yang, Ken Chen, Hongliang Ren. A minimal POE-based model for robotic kinematic calibration with only position measurements. IEEE Transactions on Automation Science and Engineering. 2015, 12(2): 758-763.
clc;
clear;

%% Prepare data
s=0.001; %unit scale for distance, s=0.001 if use m as unit, s=1 if use mm as unit

%PUMA 560 robot model, nominal parameters
xi0=[   0,  0,  0,      0,      0,      0;
        0,  -1, -1,     0,      -1,     0;
        1,  0,  0,      -1,     0,      -1;
        0,  0,  0,      -50,    -20,    -50;
        0,  0,  0,      250,    0,      250;
        0,  0,  -100,   0,      -250,   0;];
xi0(4:6,:)=xi0(4:6,:)*s;
%initial twist
xist0=[0;0;0;250;50;-20];
xist0(4:6,:)=xist0(4:6,:)*s;

%actual parameters
xi01=[0.0399999800000150;-0.0199999900000075;0.998999500500375;0.0200000000000000;0.0400000000000000;0;];
xi02=[0;-1;0;-0.0200000000000000;0;0.0500000000000000;];
xi03=[0.178005251232368;-0.984029029284552;-0.00100002950130544;-0.0841845888907446;0.0874136824072620;-100.999920311298;];
xi04=[0.0619994730067192;0.0129998895014089;-0.997991517108157;-50.9999969248523;249.000000644789;0.0751505000414993;];
xi05=[0.000999959501660306;-0.999999500040373;0;-20.6000000008239;-0.0205991760337826;-249;];
xi06=[0.0949994775043106;0.0309998295014066;-0.994994527545148;-51.2730269967030;248.910906980023;2.85959854441601;];
xi00=[xi01,xi02,xi03,xi04,xi05,xi06];
xi00(4:6,:)=xi00(4:6,:)*s;
%initial twist
xist00=[0.02;-0.01;0.01;249;51;-20.6];
xist00(4:6,:)=xist00(4:6,:)*s;

P01=[100;0;0;]*s;%nominal position of point 1
P02=[0;100;0;]*s;%nominal position of point 2
P03=[0;0;100;]*s;%nominal position of point 3
PX=[P01,P02,P03];

P04=[-100;-100;-100;]*s;%point for validation

%% calibration
mean_error_before=zeros(10,1);%mean error before calibration
mean_error_afterMOP=zeros(10,1);%mean error after calibration with MOP
mean_error_afterTOP=zeros(10,1);
mean_error_afterMFP=zeros(10,1);
mean_error_afterTFP=zeros(10,1);

max_error_before=zeros(10,1);
max_error_afterMOP=zeros(10,1);
max_error_afterTOP=zeros(10,1);
max_error_afterMFP=zeros(10,1);
max_error_afterTFP=zeros(10,1);

for N=10:10:100 %number of measurements
    disp(N)
    CaliJointAngles=rand(N,6)*2*pi-pi; %N groups of random joint angles

    gm=zeros(4,4,N);%measured end-effector poses
    gn=zeros(4,4,N);%nominal end-effector poses
    
    Pa1=zeros(3,N);%measured position of point 1
    Pa2=zeros(3,N);%measured position of point 2
    Pa3=zeros(3,N);%measured position of point 3
    
    %simulate measurement data, add random error
    for i=1:N
        Pa1(:,i)=[eye(3),zeros(3,1)]*fkPUMA560(xi00,xist00,CaliJointAngles(i,:),6)*[P01;1]+rand(3,1)*0.2*s-0.1*s;
        Pa2(:,i)=[eye(3),zeros(3,1)]*fkPUMA560(xi00,xist00,CaliJointAngles(i,:),6)*[P02;1]+rand(3,1)*0.2*s-0.1*s;
        Pa3(:,i)=[eye(3),zeros(3,1)]*fkPUMA560(xi00,xist00,CaliJointAngles(i,:),6)*[P03;1]+rand(3,1)*0.2*s-0.1*s;
        PY=[Pa1(:,i),Pa2(:,i),Pa3(:,i)];
        
        [R,t,~,~,~]=Registration(PX,PY,eye(3),zeros(3,1),1);
        gm(1:3,1:3,i)=R;
        gm(1:3,4,i)=t;
        
        gn(:,:,i)=fkPUMA560(xi00,xist00,CaliJointAngles(i,:),6);
    end
    
    [xiTFP,xistTFP]=Puma560Traditional(xi0,xist0,CaliJointAngles,gm,10); %TFP, Calibration with Traditional model and Full Pose data
    [xiMFP,xistMFP]=Puma560Minimal(xi0,xist0,CaliJointAngles,gm,10); %MFP, Calibration with Minimal model and Full Pose data
    [xiTOP,xistTOP]=Puma560TraditionalOnlyP(xi0,xist0,repmat(CaliJointAngles,3,1),P01,N,P02,N,P03,N,[Pa1,Pa2,Pa3],10); %TOP, Calibration with Traditional model and Only Position data
    [xiMOP,xistMOP]=Puma560MinimalOnlyP(xi0,xist0,repmat(CaliJointAngles,3,1),P01,N,P02,N,P03,N,[Pa1,Pa2,Pa3],10); %MOP, Calibration with Minimal model and Only Position data

% evaluation
    Nt=40; %number of test data
    TestJointAngles=[rand(Nt,1)*2*pi,rand(Nt,1)*2*pi,rand(Nt,1)*2*pi,rand(Nt,1)*2*pi,rand(Nt,1)*2*pi,rand(Nt,1)*2*pi]; %Nt groups of test configurations
    
    error_before=zeros(1,Nt);error_afterMOP=zeros(1,Nt);error_afterTOP=zeros(1,Nt);error_afterMFP=zeros(1,Nt);error_afterTFP=zeros(1,Nt);
    
    for i=1:Nt
        error_before(i)=norm((fkPUMA560(xi0,xist0,TestJointAngles(i,:),6)-fkPUMA560(xi00,xist00,TestJointAngles(i,:),6))*[P04;1])/s;
        error_afterMOP(i)=norm((fkPUMA560(xiMOP,xistMOP,TestJointAngles(i,:),6)-fkPUMA560(xi00,xist00,TestJointAngles(i,:),6))*[P04;1])/s;
        error_afterTOP(i)=norm((fkPUMA560(xiTOP,xistTOP,TestJointAngles(i,:),6)-fkPUMA560(xi00,xist00,TestJointAngles(i,:),6))*[P04;1])/s;
        error_afterMFP(i)=norm((fkPUMA560(xiMFP,xistMFP,TestJointAngles(i,:),6)-fkPUMA560(xi00,xist00,TestJointAngles(i,:),6))*[P04;1])/s;
        error_afterTFP(i)=norm((fkPUMA560(xiTFP,xistTFP,TestJointAngles(i,:),6)-fkPUMA560(xi00,xist00,TestJointAngles(i,:),6))*[P04;1])/s;
    end
    mean_error_before(N/10,1)=mean(error_before);
    mean_error_afterMOP(N/10,1)=mean(error_afterMOP);
    mean_error_afterTOP(N/10,1)=mean(error_afterTOP);
    mean_error_afterMFP(N/10,1)=mean(error_afterMFP);
    mean_error_afterTFP(N/10,1)=mean(error_afterTFP);

    max_error_before(N/10,1)=max(error_before);
    max_error_afterMOP(N/10,1)=max(error_afterMOP);
    max_error_afterTOP(N/10,1)=max(error_afterTOP);
    max_error_afterMFP(N/10,1)=max(error_afterMFP);
    max_error_afterTFP(N/10,1)=max(error_afterTFP);

end

%% result
figure(1)
hold off;plot(10:10:100,mean_error_afterMOP,'r-','LineWidth',1.5)
hold on; plot(10:10:100,mean_error_afterTOP,'g-','LineWidth',1.5)
hold on; plot(10:10:100,mean_error_afterMFP,'m-','LineWidth',1.5)
hold on; plot(10:10:100,mean_error_afterTFP,'b-','LineWidth',1.5)
hold on; plot(10:10:100,max_error_afterMOP,'r:','LineWidth',1.5)
hold on; plot(10:10:100,max_error_afterTOP,'g:','LineWidth',1.5)
hold on; plot(10:10:100,max_error_afterMFP,'m:','LineWidth',1.5)
hold on; plot(10:10:100,max_error_afterTFP,'b:','LineWidth',1.5)

xlabel('Number of Mesurement Configurations')
ylabel('Position Error of the Test Point (mm)')
legend('Mean(MOP)','Mean(TOP)','Mean(MFP)','Mean(TFP)','Max(MOP)','Max(TOP)','Max(MFP)','Max(TFP)')