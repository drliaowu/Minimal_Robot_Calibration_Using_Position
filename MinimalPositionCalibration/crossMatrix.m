function [ R ] = crossMatrix( a )
% [ R ] = crossMatrix( a ) gives skew matrix of a vector (a).
R=[0,-a(3),a(2);
   a(3),0,-a(1);
  -a(2),a(1),0];

end

