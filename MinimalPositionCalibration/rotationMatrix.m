function [ R ] = rotationMatrix( w,theta )
% [ R ] = rotationMatrix( w,theta ) generates a rotation matrix (R) from unit axis
% direction (w) and rotation angle (theta)

R=eye(3)+crossMatrix(w)*sin(theta)+crossMatrix(w)*crossMatrix(w)*(1-cos(theta));

end

