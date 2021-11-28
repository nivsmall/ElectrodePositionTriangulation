function [dist] = euc_dist3D(coord_mat)
%   distance from the electrode at (0,0)
dist = sqrt(sum(coord_mat.^2, 2));
end

