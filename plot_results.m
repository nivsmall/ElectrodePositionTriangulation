function [] = plot_results(trueLabels,initialization, optResults ,titleStr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global len
X_results = optResults(1:len, 1);
Y_results = optResults(1:len, 2);

X_TLabel = trueLabels(1:len, 1);
Y_TLabel = trueLabels(1:len, 2);

X_init = initialization(1:len, 1);
Y_init = initialization(1:len, 2);

figure
plot(X_init, Y_init, 'sk');
hold on
plot(X_TLabel, Y_TLabel, 'ob');
plot(X_results, Y_results,'*r');
legend('Optimization Initialization', 'True Label', 'Optimization Results');
title(titleStr);
hold off;
end

