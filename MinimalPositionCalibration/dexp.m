function [ J ] = dexp( kesi, theta )
%[ J ] = dexp( kesi, theta ) : Jacobian matrix for a joint (see (13) and (14))
J=[aMatrix(kesi,theta),kesi];

end

