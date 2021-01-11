%%
% Please refer to Section IV, EXPERIMENTS ON AN ABB IRB 120 ROBOT
% Liao Wu, Xiangdong Yang, Ken Chen, Hongliang Ren. A minimal POE-based model for robotic kinematic calibration with only position measurements. IEEE Transactions on Automation Science and Engineering. 2015, 12(2): 758-763.
clc;
clear;

%% measurement data processing
s=0.001;%unit scale for distance, s=0.001 if use m as unit, s=1 if use mm as unit

%BaseCirclePoints, measure a fixed point on the end-effector when
%rotating the first axis, to find the Z direction of base frame wrt
%measurement frame
load('BaseCirclePoints.mat');
BaseCirclePoints=BaseCirclePoints*s;

%CaliPosition=[point1_1;point2_1;point3_1;point1_2;point2_2;point3_2;...;point1_36;point2_36;point3_36]
%36 measured configurations; in each configuration 3 points are measured wrt measurement frame 
load('CaliPosition.mat');
CaliPosition=CaliPosition*s;

%CaliFullPose, 37 measured configurations, the first 36 measured
%configurations are random, the last measured configuration is at the
%zero-reference configuration. In each configuration, 8 points on a plane
%are measured to construct the XY plane and 4 points evenly distributed around
%the centre are measured to construct the the Z axis 
load('CaliFullPose.mat');
CaliFullPose=CaliFullPose*[s*eye(3),zeros(3);zeros(3),eye(3)];

%36 groups of joint angles for calibration
load('CaliJointAngles.mat');
CaliJointAngles=CaliJointAngles*pi/180;

%first 10 measurements of a point on end-effector for test
load('TestPosition.mat');
TestPosition=TestPosition*s;

%another 30 measurements of a point on end-effector for test
load('TestPosition2.mat');
TestPosition2=TestPosition2*s;

%first 10 groups of joint angles for test
load('TestJointAngles.mat');
TestJointAngles=TestJointAngles*pi/180;

TestPoint_Meas=[-82.538708	502.156798	652.749355]*s;%Test Point 1

%another 30 groups of joint angles for test
load('TestJointAngles2.mat');
TestJointAngles2=TestJointAngles2*pi/180;

TestPoint_Meas2=[-82.389469	502.253113	654.004118]*s;%Test Point 2

%base frame measured wrt measurement frame
BaseDirZ_Meas=[-0.000289	-0.002459	0.999997];
BasePZ_Meas=[-126.930293	871.797701	11.929126]*s;
BaseDirX_Meas=[0.140895     -0.990024	0.00132];

%find transformation from measurement frame to base frame
DirOZ=BaseDirZ_Meas/norm(BaseDirZ_Meas);
DirOY=cross(DirOZ,BaseDirX_Meas);
DirOY=DirOY/norm(DirOY);
DirOX=cross(DirOY,DirOZ);
[~,Circle]=FindTwist(BaseCirclePoints(1:6,:));
O=Circle+DirOZ*(BasePZ_Meas'-Circle)*DirOZ';
R=[DirOX',DirOY',DirOZ'];
Tmb=[R,O;0,0,0,1];%transformation from measurement frame to base frame

%nominal positions of three points measured for calibration wrt measurement
%frame
CaliPoint_Meas=[-82.401255	502.191363	652.794976	;
-69.225085	489.81311	653.858786	;
-60.642945	505.693111	653.060532	;]*s;%

%transform data from measurement frame to base frame
%for methods using only position, base frame = end-effector frame
CaliPosition_Base=(CaliPosition-repmat(O',108,1))*R;
TestPosition_Base=(TestPosition-repmat(O',10,1))*R;
TestPosition_Base2=(TestPosition2-repmat(O',30,1))*R;

TestPoint_Base=(TestPoint_Meas-O')*R;
TestPoint_Base2=(TestPoint_Meas2-O')*R;
CaliPoint_Base1=(CaliPoint_Meas(1,:)-O')*R;
CaliPoint_Base2=(CaliPoint_Meas(2,:)-O')*R;
CaliPoint_Base3=(CaliPoint_Meas(3,:)-O')*R;

%process full pose measurements
%measured end-effector frame for calibration
EndFrame=zeros(4,4,36);
for i=1:36
    OP=CaliFullPose((i-1)*12+8,1:3);
    DZ=CaliFullPose((i-1)*12+8,4:6);
    DZ=DZ/norm(DZ);
    DX=CaliFullPose((i-1)*12+9,1:3)-OP;
    DY=cross(DZ,DX);
    DY=DY/norm(DY);
    DX=cross(DY,DZ);
    EndFrame(1:4,1:4,i)=Tmb\[DX',DY',DZ',OP';0,0,0,1];
end

%measured initial end-effector frame
EndFrame0=zeros(4,4);
OP=CaliFullPose((37-1)*12+8,1:3);
DZ=CaliFullPose((37-1)*12+8,4:6);
DZ=DZ/norm(DZ);
DX=CaliFullPose((37-1)*12+9,1:3)-OP;
DY=cross(DZ,DX);
DY=DY/norm(DY);
DX=cross(DY,DZ);
EndFrame0(1:4,1:4)=Tmb\[DX',DY',DZ',OP';0,0,0,1];

%transform data from base frame to end-effector frame
TestPoint_End=EndFrame0(1:3,1:3)\(TestPoint_Base'-EndFrame0(1:3,4));
TestPoint_End2=EndFrame0(1:3,1:3)\(TestPoint_Base2'-EndFrame0(1:3,4));
xist0_full=vlog(EndFrame0);

%nominal robot model
w1=[0;0;1];p1=[0;0;0]*s;
w2=[0;1;0];p2=[0;0;290]*s;
w3=[0;1;0];p3=[0;0;560]*s;
w4=[1;0;0];p4=[0;0;630]*s;
w5=[0;1;0];p5=[302;0;630]*s;
w6=[1;0;0];p6=[0;0;630]*s;
xi0=[twistCoord(w1,p1),twistCoord(w2,p2),twistCoord(w3,p3),twistCoord(w4,p4),twistCoord(w5,p5),twistCoord(w6,p6)];
xist0=[0;0;0;0*s;0*s;0*s];

%% calibration

%calibration with only position
[xiTOP,xistTOP]=Puma560TraditionalOnlyP(xi0,xist0,[CaliJointAngles;CaliJointAngles;CaliJointAngles],CaliPoint_Base1',36,CaliPoint_Base2',36,CaliPoint_Base3',36,transpose([CaliPosition_Base(1:3:108,:);CaliPosition_Base(2:3:108,:);CaliPosition_Base(3:3:108,:);]),10);
[xiMOP,xistMOP]=Puma560MinimalOnlyP(xi0,xist0,[CaliJointAngles;CaliJointAngles;CaliJointAngles],CaliPoint_Base1',36,CaliPoint_Base2',36,CaliPoint_Base3',36,transpose([CaliPosition_Base(1:3:108,:);CaliPosition_Base(2:3:108,:);CaliPosition_Base(3:3:108,:);]),10);

%calibration with full pose
[xiTFP,xistTFP]=Puma560Traditional(xi0,xist0_full,CaliJointAngles,EndFrame,10);
[xiMFP,xistMFP]=Puma560Minimal(xi0,xist0_full,CaliJointAngles,EndFrame,10);

%% evaluation
P_before=zeros(4,10); P_afterMOP=zeros(4,10); P_afterTOP=zeros(4,10); P_afterMFP=zeros(4,10); P_afterTFP=zeros(4,10);
P_before2=zeros(4,30); P_afterMOP2=zeros(4,30); P_afterTOP2=zeros(4,30); P_afterMFP2=zeros(4,30); P_afterTFP2=zeros(4,30);

error_before=zeros(1,10); error_afterMOP=zeros(1,10); error_afterTOP=zeros(1,10); error_afterMFP=zeros(1,10); error_afterTFP=zeros(1,10);
error_before2=zeros(1,30); error_afterMOP2=zeros(1,30); error_afterTOP2=zeros(1,30); error_afterMFP2=zeros(1,30); error_afterTFP2=zeros(1,30);

for i=1:10
    P_before(:,i)=fkPUMA560(xi0,xist0,TestJointAngles(i,:),6)*[TestPoint_Base';1];
    error_before(i)=norm(P_before(1:3,i)-transpose(TestPosition_Base(i,:)));
    P_afterTFP(:,i)=fkPUMA560(xiTFP,xistTFP,TestJointAngles(i,:),6)*[TestPoint_End;1];
    error_afterTFP(i)=norm(P_afterTFP(1:3,i)-transpose(TestPosition_Base(i,:)));
    P_afterMFP(:,i)=fkPUMA560(xiMFP,xistMFP,TestJointAngles(i,:),6)*[TestPoint_End;1];
    error_afterMFP(i)=norm(P_afterMFP(1:3,i)-transpose(TestPosition_Base(i,:)));
    P_afterTOP(:,i)=fkPUMA560(xiTOP,xistTOP,TestJointAngles(i,:),6)*[TestPoint_Base';1];
    error_afterTOP(i)=norm(P_afterTOP(1:3,i)-transpose(TestPosition_Base(i,:)));
    P_afterMOP(:,i)=fkPUMA560(xiMOP,xistMOP,TestJointAngles(i,:),6)*[TestPoint_Base';1];
    error_afterMOP(i)=norm(P_afterMOP(1:3,i)-transpose(TestPosition_Base(i,:)));
end

for i=1:30
    P_before2(:,i)=fkPUMA560(xi0,xist0,TestJointAngles2(i,:),6)*[TestPoint_Base2';1];
    error_before2(i)=norm(P_before2(1:3,i)-transpose(TestPosition_Base2(i,:)));
    P_afterTFP2(:,i)=fkPUMA560(xiTFP,xistTFP,TestJointAngles2(i,:),6)*[TestPoint_End2;1];
    error_afterTFP2(i)=norm(P_afterTFP2(1:3,i)-transpose(TestPosition_Base2(i,:)));
    P_afterMFP2(:,i)=fkPUMA560(xiMFP,xistMFP,TestJointAngles2(i,:),6)*[TestPoint_End2;1];
    error_afterMFP2(i)=norm(P_afterMFP2(1:3,i)-transpose(TestPosition_Base2(i,:)));
    P_afterTOP2(:,i)=fkPUMA560(xiTOP,xistTOP,TestJointAngles2(i,:),6)*[TestPoint_Base2';1];
    error_afterTOP2(i)=norm(P_afterTOP2(1:3,i)-transpose(TestPosition_Base2(i,:)));
    P_afterMOP2(:,i)=fkPUMA560(xiMOP,xistMOP,TestJointAngles2(i,:),6)*[TestPoint_Base2';1];
    error_afterMOP2(i)=norm(P_afterMOP2(1:3,i)-transpose(TestPosition_Base2(i,:)));
end

error_beforet=[error_before,error_before2]/s;%before
error_afterMFPt=[error_afterMFP,error_afterMFP2]/s;%MFP
error_afterTFPt=[error_afterTFP,error_afterTFP2]/s;%TFP
error_afterTOPt=[error_afterTOP,error_afterTOP2]/s;%MOP
error_afterMOPt=[error_afterMOP,error_afterMOP2]/s;%TOP

%% result
figure(1)
bar([mean(error_beforet),max(error_beforet);mean(error_afterMOPt),max(error_afterMOPt);mean(error_afterTOPt),max(error_afterTOPt);mean(error_afterMFPt),max(error_afterMFPt);mean(error_afterTFPt),max(error_afterTFPt)])
set(gca,'XTickLabel',{'Before Calibration','MOP','TOP', 'MFP', 'TFP'})
ylabel('PositionError of the Test Point (mm)')
legend('Mean Error','Max Error')