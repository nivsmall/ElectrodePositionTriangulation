function [lb,ub,A,B,Aeq,Beq] = initBoundaries(coordCoeff, ACoeff, alphaCoeff)
% initialize boundaries for optimization
global Coordinates alpha_mu len coord_coeff  A_coeff  alpha_coeff;

Coord_lb = Coordinates(:, 1:2) - 3500.*ones(len,2);
Coord_ub = Coordinates(:, 1:2) + 3500.*ones(len,2);
A_Vect_lb = zeros(len, 1);
A_Vect_ub = 20000*ones(len, 1);
alpha_lb = 20*alpha_mu;
alpha_ub = 0.005*alpha_mu;
A =         [];
B =         [];
Aeq =       [];
Beq =       [];
if coord_coeff && A_coeff && alpha_coeff
    lb = [Coord_lb A_Vect_lb; alpha_lb 0 0];
    ub = [Coord_ub A_Vect_ub; alpha_ub 0 0];
elseif coord_coeff && A_coeff && ~alpha_coeff
    lb = [Coord_lb A_Vect_lb];
    ub = [Coord_ub A_Vect_ub];
elseif coord_coeff && ~A_coeff && ~alpha_coeff
    lb =  Coord_lb;
    ub =  Coord_ub;
elseif alpha_coeff && ~coord_coeff && ~A_coeff
    lb = alpha_lb;
    ub = alpha_ub;
elseif alpha_coeff && ~coord_coeff && A_coeff
    lb = [A_Vect_lb; alpha_lb];
    ub = [A_Vect_ub; alpha_ub];
elseif ~coord_coeff && A_coeff && ~alpha_coeff
    lb = A_Vect_lb;
    ub = A_Vect_ub;
end

end

