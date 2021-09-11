function fun=errorFun_2d(x)
% the computed loss/objective function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   ln(Vij-Vo)-ln(Aj)+ln(rij)*alpha                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Pot_Mat len alpha_mu A_Vect v0 coord_coeff A_coeff alpha_coeff distFilterMask Coordinates CSDTrain CSD;

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
    r = dist(Coordinates');
    A_Vector = A_Vect;
    alpha = x;
elseif alpha_coeff && ~coord_coeff && A_coeff
    r = dist(Coordinates');
    A_Vector = x(1:len);
    alpha = x(len+1);
elseif  ~alpha_coeff && ~coord_coeff && A_coeff
    r = dist(Coordinates');
    A_Vector = x(1:len);
    alpha = alpha_mu;
end


for idx = 1:len
    for jdx = 1:len
        if idx == jdx
            continue
        end
        if distFilterMask(idx, jdx)>0
            continue
        end
        if CSDTrain
            %prefun(idx,jdx) = CSD(idx, jdx) - Pot_Mat(idx, jdx)*A_Vector(jdx)*r(idx, jdx);
            prefun(idx, jdx) = CSD(idx, jdx) - A_Vector(jdx)*r(idx, jdx)*r(idx, jdx)^alpha;
        else
            prefun(idx,jdx) = log(abs(Pot_Mat(idx,jdx)*10^-12)) - log(A_Vector(jdx)) - alpha*log(r(idx,jdx));
        end
    end
end
fun = sum(sum(prefun.^2));
end
