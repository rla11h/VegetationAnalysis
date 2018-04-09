function K = StiffnessProfile(theta,k,H)
% Function generates a line to define a stiffness profile of the stem based
% on the two stiffness values calculated during each hold phase. Solving
% for the equation of the line, it returns the expected stiffness for the
% current amount of bending

% Inputs are two stiffness values, k, two matching angles, H, and the angle
% at which to evaluate, theta

      m = (k(2)-k(1))/(H(2)-H(1));
      b = k(2)-m*H(2);
      
      K = m*theta+b;
end
