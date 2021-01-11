function [ T ] = se3Translation( v,theta )
%[ T ] = se3Translation( v,theta ) generates a pure translation

T=[eye(3),v*theta;
    0,0,0,1];
end

