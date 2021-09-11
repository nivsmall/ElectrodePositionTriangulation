function [x, y, cleaned] = clean_outliers(x, y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[x_sorted,xorder]=sort(x);
y_sorted = y(xorder);
cleaned = 0;
values = unique(x_sorted);
cnt=0;
for dist = 1:length(values)
    repeats=0;
    for idx = 1:length(x_sorted)
        if x_sorted(cnt+idx) == values(dist)
            repeats = repeats+1;
            continue;
        else
            tempV = y_sorted(cnt+idx-repeats:idx-1);
            if length(tempV)>2
               [~, TF] = rmoutliers(tempV);
               cleaned = cleaned + sum(TF);
               tempV = tempV.*~TF;
               meanVal = sum(tempV)/(length(tempV)-sum(TF));
               tempV = meanVal.*TF+tempV;
               y_sorted(cnt+idx-repeats:idx-1) = tempV;
               cnt = cnt + idx-1;
            end
        end
        break;
    end
    
end
cleaned
unsorted = 1:length(x);
original_order(xorder)=unsorted;
x = x_sorted(original_order);
y = y_sorted(original_order);
end

