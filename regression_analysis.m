function [alpha, distances, residuals] = regression_analysis(Coordinates,Pot_Mat, len, show)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global needs_clean A_Vect v0;
neads_clean = zeros(len);   % binary - out of acccepted range data cleaning
r = dist(Coordinates.');
v=Pot_Mat;
r(logical(eye(length(r)))) = 0;
v(logical(eye(length(Pot_Mat)))) = inf;

% model function: Vij = Aj*rij^(alpha) + Const
modelfun = @(c,x)(c(1).*(x.^c(2)))+c(3);
modelfit = cell(len,1);
reg_Coef = zeros(len,3);
R_sq = zeros(len,1);
for idx = 1:len
    x = r(:,idx).*10^-6;
    y = v(:,idx).*10^-12;
    %[x, y, cleaned(idx)]=clean_outliers(x, y);
    %Pot_Mat(:,idx) = y*10^12;
    %[x, y, needs_clean(:,idx)]=clean_out_of_range_potentials([10^-4,10^-3], x,y);
    test = 1;
    while test <= length(y)
        if abs(y(test))==inf
            y(test) = [];
            x(test) = [];
            continue
        end
        test = test+1;
    end
    %x(idx) = [];
    %y(idx) = [];
                        %initial values:[Aj,alpha,Const]
    modelfit{idx} = fitnlm(x,y,modelfun,[1e-4 -2 -10^-5]);
    reg_Coef(idx,:) = modelfit{idx}.Coefficients.Estimate';
    R_sq(idx) = modelfit{idx}.Rsquared.Adjusted;
    
    yhat=modelfit{idx}.Fitted;
    if show
        plot(x,yhat)
        hold on
        plot (x, y, '.r')
        hold off
    end
 end
residuals = zeros(len-1);
for idx = 1:len
    tab = table2array(modelfit{idx}.Variables);
    [~,xorder]=sort(tab(:,1));
    restab = table2array(modelfit{idx}.Residuals(:,1));
    residuals(idx,:) = restab(xorder);
    
end
if show
    figure
    plot(mean(abs(residuals)))
    title("residuals plot (ordered from high potential measurements to low)")
end


res_sum = [];
sum_size = [];
domain = [];

% since a distance between electrodes can repeat, lets average the fit for
% those cases:
for idx = 1:len
   range = table2array(modelfit{idx}.Variables(:,1)); 
   res = table2array(modelfit{idx}.Residuals(:,1));
   for jdx = 1:len-1
       if isempty(find(domain==range(jdx)))
           dom_size = length(domain);
           domain(dom_size+1) = range(jdx);
           sum_size(dom_size+1) = 1;
           res_sum(dom_size+1) = abs(res(jdx));
       else
           sum_size(find(domain==range(jdx))) = sum_size(find(domain==range(jdx))) +1;
           res_sum(find(domain==range(jdx))) = res_sum(find(domain==range(jdx))) + res(jdx);
       end
   end
end
res_avg = res_sum./sum_size;
[~,domain_order] = sort(domain);
%figure
%plot(domain(domain_order),res_avg(domain_order))
%title("Residual vs electrode distance plot")

distances = domain(domain_order).*10^6;
residuals = res_avg(domain_order);

%%plotting a set of regressions
if show
    for reg_num = 33:34
        figure
        tab = table2array(modelfit{reg_num}.Variables);
        [x,xorder]=sort(tab(:,1));
        y = tab(xorder,2);
        yhat = modelfit{reg_num}.Fitted(xorder);
        plot(x,y,'.r');
        hold on
        t = plot(x,yhat);
        title(['Regression of potential decay from electrode ' num2str(reg_num)])
        xlabel('Distance [m]')
        ylabel('Potential [V]')
        %saveas(t, ['C:\Users\nivsm\Documents\FinalProject\images and plots\New folder\not clean\Regression of potential decay from electrode ' num2str(reg_num) '.png'])
    end
end
% [alpha0, alpha0_pci]: 
alpha = mle(reg_Coef(:,2));
v0 = mle(reg_Coef(:,3));
A_Vect = reg_Coef(:, 1);
end

