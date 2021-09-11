function [x,y, cleaned_indices] = clean_out_of_range_potentials(range,x, y)
cleaned_indices=zeros(length(y),1);
for idx = 1:length(y)
   if range(1)<y(idx) && y(idx)<range(2)  
   else
       disp('value out of range: ( X, Y ):');
       fprintf('\t\t ( %d, %d ) \n',x(idx), y(idx));
       y(idx)=[NaN];
       cleaned_indices(idx,1)=1;
   end
end
end

