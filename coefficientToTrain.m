function [coordCoeff, ACoeff,alphaCoeff, STR] = coefficientToTrain(coeffString)
% choosing coefficient to train in current optimization run

coordCoeff=false;
ACoeff=false;
alphaCoeff=false;

if coeffString == "coord"
    coordCoeff=true;
    STR = "->coordinates";
elseif coeffString == "A"
    ACoeff=true;
    STR = "->coupling";
elseif coeffString == "alpha"
    alphaCoeff=true;
    STR = "->decay";
end

end

