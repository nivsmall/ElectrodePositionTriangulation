function [] = plot_single_train_results(trueLabels,initialization, optResults, Coord_lb, Coord_ub)

    X_TLabel = trueLabels(:, 1);
    Y_TLabel = trueLabels(:, 2);

    X_init = initialization(1);
    Y_init = initialization(2);
 
    X_results = optResults(1);
    Y_results = optResults(2);

    
    figure
    
    plot(X_TLabel, Y_TLabel, '.b');
    hold on
    
    plot(X_init, Y_init, 'og');
       
    plot([Coord_lb(1) Coord_ub(1) Coord_ub(1) Coord_lb(1)], [Coord_lb(2) Coord_ub(2) Coord_lb(2) Coord_ub(2)]);
    
    plot(X_results, Y_results, '*r');
    
    legend('True Label', 'Optimization Initialization', 'Lower&Upper Bounds', 'Optimization Results');
    hold on
end

