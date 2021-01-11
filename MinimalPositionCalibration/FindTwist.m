function [ xi,origin] = FindTwist( P )
%[ xi,origin] = FindTwist( P ) gives twist from circle points
N=size(P,1); %number of points
%fitplane
b=-1*ones(N,1);
x_inv=(P'*P)\(P'*b);
x_inv=x_inv/norm(x_inv);
if(sum(x_inv)<0)
    x_inv=-x_inv;
end
w=x_inv;%xi=[w,v];

z_new=x_inv;
x_new=cross(z_new,[1;1;1]);
x_new=x_new/norm(x_new);
y_new=cross(z_new,x_new);
C_new=[x_new,y_new,z_new];

P_reg=P*C_new;
Circle=CircleFitByKasa(P_reg(:,1:2));
origin_reg=[Circle(1);Circle(2);mean(P_reg(:,3))];
origin=C_new*origin_reg;
v=cross(origin,w);
xi=[w;v];
end

