function [] = plotLossLandscape(titleStr)
%UNTITLED3 Summary of this function goes here
%%
    global Coordinates chooser alpha_mu oneAtATime;
    dAlpha=1;
    alpha_mu = alpha_mu*dAlpha;
    minVal = 10^20;
    %%
    testBoxSize = 2000;
    pointsCount = 101;   % number of test points will be pointsCount^2
    if oneAtATime
        testPointsMat = createCoordTestPoints(Coordinates(chooser, :), testBoxSize, pointsCount);
    else
        testPointsMat = createCoordTestPoints([0, 0], testBoxSize, pointsCount);
    end
    for ydx = 1:length(testPointsMat(:,1,1))
        for xdx = 1:length(testPointsMat(1,:,1))
            x = testPointsMat(ydx, xdx, 1);
            y = testPointsMat(ydx, xdx, 2);
            if oneAtATime
                testPointsMat(ydx,xdx,3) = ObjectiveFunctionOneAtaTime([x, y]);
            else
                ElecToView = 11;
                testCoord = Coordinates(:,1:2);
                testCoord(ElecToView, 1:2) = testCoord(ElecToView, 1:2)+[x, y];
                localError = errorFun_2d(testCoord);
                testPointsMat(ydx,xdx,3) = localError;
                if localError < minVal
                    minVal = localError;
                    globalMinimum = testCoord(ElecToView, 1:2);
                    Y=ydx;
                    X=xdx;
                end
            end
        end
    end
    figure;
    s = surf(testPointsMat(:,:,1), testPointsMat(:,:,2), testPointsMat(:,:,3));
    id=sub2ind(size(testPointsMat(:,:,3)),Y,X);
    datatip(s, 'DataIndex', id);
    title(titleStr);
    legend('Global Minimum');
    alpha_mu = alpha_mu/dAlpha;
end

