function fun=errorFun_2d(x)
% the computed loss/objective function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   ln(Vij-Vo)-ln(Aj)+ln(rij)*alpha                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Pot_Mat len x0_init coord_coeff A_coeff alpha_coeff distFilterMask  CSDTrain CSD alpha_mu A_Vect Coordinates;

prefun = zeros(len);
% exclude noisy measurements:
%Pot_Mat = Pot_Mat.*~needs_clean + needs_clean.*(1+v0).*10^12;

if coord_coeff && A_coeff && alpha_coeff
    r = dist(x(1:len,1:2)');
    A_Vector = x(1:len,3);
    alpha = x(len+1);
elseif coord_coeff && A_coeff && ~alpha_coeff
    r = dist(x(1:len,1:2)');
    A_Vector = x(1:len, 3);
    alpha = x0_init(len+1);
elseif coord_coeff && ~A_coeff && ~alpha_coeff
    r = dist(x(1:len,1:2)');
    A_Vector = x0_init(1:len,3);
    alpha = x0_init(len+1);
elseif alpha_coeff && ~coord_coeff && ~A_coeff
    r = dist(x0_init(1:len,1:2)');
    A_Vector = x0_init(1:len,3);
    alpha = x;
elseif alpha_coeff && ~coord_coeff && A_coeff
    r = dist(x0_init(1:len,1:2)');
    A_Vector = x(1:len);
    alpha = x(len+1);
elseif  ~alpha_coeff && ~coord_coeff && A_coeff
    r = dist(x0_init(1:len,1:2)');
    A_Vector = x(1:len);
    alpha = x0_init(len+1);
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
        else %+3.5
            prefun(idx,jdx) = log(abs(Pot_Mat(idx,jdx)*10^-(12+3.5))) - log(A_Vector(jdx)) - alpha*log(r(idx,jdx));
        end
    end
end
fun = sum(sum(prefun.^2));
end
