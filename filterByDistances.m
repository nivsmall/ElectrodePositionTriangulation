function [mask] = filterByDistances(r, min, max)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global Pot_Mat distances;
mask = zeros(size(Pot_Mat));
for d=1:length(distances)
    if distances(d)<min || distances(d)>max
        mask0 = r == distances(d);
        mask = mask+mask0;
    end
end
end

