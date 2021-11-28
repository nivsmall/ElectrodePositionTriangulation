function [] = addNoiseManually(r)
%   add noise to simulated data according to MEA data <--> regression_analysis(Coordinates, Pot_MatMea, len)
global distances res Pot_Mat;
figure
plot(distances,res)
title("Residual vs electrode distance plot")
xlabel('Distance [m]')
ylabel('MRE [V]')
    for d = 1:length(distances)
        mask = r == distances(d);
        Pot_Mat = Pot_Mat - mask.*res(d)*10^12;
    end
end

