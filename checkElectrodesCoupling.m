function [PearsonCorrelation] = checkElectrodesCoupling(PotMat, A_Vect)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
side1 = reshape(mean(PotMat, 1), length(PotMat), 1);
side2 = mean(PotMat, 2);
final = (side1+side2)./2;
final = final./mean(final);
corr = corrcoef(A_Vect, final);
PearsonCorrelation = corr(1, 2);
end

