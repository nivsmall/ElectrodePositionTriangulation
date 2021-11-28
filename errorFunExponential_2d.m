function [fun] = errorFunExponential_2d(x)
% the computed loss/objective function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   ln(Vij)-ln(Aj)+ln(alpha)*rij                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Pot_Mat len alpha_mu A_Vect v0 coord_coeff A_coeff alpha_coeff;

prefun = zeros(len);
% exclude noisy measurements:
%Pot_Mat = Pot_Mat.*~needs_clean + needs_clean.*(1+v0).*10^12;

if coord_coeff && A_coeff && alpha_coeff
    r = dist(x(1:len,1:2)');
    A_Vector = A_Vect;
    alpha = x(len+1);
elseif coord_coeff && A_coeff && ~alpha_coeff
    r = dist(x(1:len,1:2)');
    A_Vector = x(1:len, 3);
    alpha = alpha_mu;
elseif coord_coeff && ~A_coeff && ~alpha_coeff
    r = dist(x(1:len,1:2)');
    A_Vector = A_Vect;
    alpha = alpha_mu;
elseif alpha_coeff && ~coord_coeff && ~A_coeff
end

for idx = 1:len
    for jdx = 1:len
        if idx == jdx
            continue
        end
        prefun(idx,jdx) = log(abs(Pot_Mat(idx,jdx))) - log(A_Vector(jdx)) - r(idx,jdx)*log(-alpha);
    end
end
fun = sum(sum(prefun.^2));
end

