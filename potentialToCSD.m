function [CSD] = potentialToCSD(PotMat, len, Coordinates)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
spatialPotMat = zeros(10, 6);
plotMat = nan(10, 6, 3);
CSD = nan(len);
for jdx = 1:len
    for idx = 1:len
        if idx == jdx
        %    spatialPotMat(idx) = inf;
            Y = Coordinates(idx, 1)/500 + 1;
            X = Coordinates(idx, 2)/500 + 1;
        end
        y = Coordinates(idx, 1)/500 + 1;
        x = Coordinates(idx, 2)/500 + 1;
        spatialPotMat(y, x) = PotMat(jdx, idx);
        plotMat(y,x,1) = Coordinates(idx, 1);
        plotMat(y,x, 2) = Coordinates(idx, 2);
        plotMat(y, x, 3) = spatialPotMat(y, x);
    end
    
    %plotMat(6,1,3) = (spatialPotMat(5,1)+spatialPotMat(7,1))/2;
    %spatialPotMat(6,1) = plotMat(6,1,3);
    plotMat(6,1,1) = (plotMat(5, 1, 1)+plotMat(7, 1, 1))/2;
    plotMat(6,1,2) = (plotMat(5, 1, 2)+plotMat(7, 1, 2))/2;
    %plotMat(Y,X,3) = ( spatialPotMat(Y-1,X)+spatialPotMat(Y,X-1)+spatialPotMat(Y,X+1)+spatialPotMat(Y+1,X) ) / 4;
%     plotMat(Y, X,3) = inf;
%     figure;
%     surfacePlot = surf(plotMat(:,:, 1), plotMat(:,:,2), plotMat(:,:,3));
%     daspect([500 500 10^8]);
%     if Y<10
%         datatip(surfacePlot, 'DataIndex', (X-1)*10+Y+1);
%     end
%     if X<6
%         datatip(surfacePlot, 'DataIndex', (X)*10+Y);
    
    actuatorStr = '('+string(plotMat(Y,X,1))+'_'+string(plotMat(Y,X,2))+')'; 
    %saveas(surfacePlot, 'C:\Users\nivsm\Documents\FinalProject\images and plots\PotentialPlots\PotentialMeasuremetns_ActuationBy'+actuatorStr+'.png');
    laplacian_ = del2(spatialPotMat, 500, 500);
    %figure;
    %surfacePlot = surf(plotMat(:,:, 1), plotMat(:,:,2), laplacian_);
    %saveas(surfacePlot, 'C:\Users\nivsm\Documents\FinalProject\images and plots\CSDPlots\CSDMeasuremetns_ActuationBy'+actuatorStr+'.png');


for idx = 1:len
        y =  Coordinates(idx, 1)/500 + 1;
        x = Coordinates(idx, 2)/500 + 1;
        CSD(idx, jdx) = laplacian_(y, x);
    end
end

end

