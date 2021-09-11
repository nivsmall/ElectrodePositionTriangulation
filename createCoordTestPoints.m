function [testPointsMat] = createCoordTestPoints(truePosition,boxSize, n)
%   testPointsMat is a 3 dimensional matrix
%   1st dim represents X's (elements in coloumns identical)
%   2nd dim represents Y's (elements in rows identical)
%   3rd dim will represent the loss at (X, Y) at adjacent indices
xTests = linspace(truePosition(1)-boxSize/2, truePosition(1)+boxSize/2, n);
xTests = ones(n,1)*xTests;
yTests = linspace(truePosition(2)-boxSize/2, truePosition(2)+boxSize/2, n);
yTests = yTests'*ones(1, n);
testPointsMat = zeros(length(yTests), length(xTests), 3);
testPointsMat(:,:,1) = xTests;
testPointsMat(:, :, 2) = yTests;
end

