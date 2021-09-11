clear
global Pot_Mat len A_Vect;
num_trials = 10;
load_data = true;
read_data = true;

%% Read the data:
%%%% Data Flow prior MATLAB: MCRack -> .cmd file -> MCDataManager -> .h5 file %%%%
if load_data
    samplerate = 1/10000;       % 10[kHz]
    if read_data
        load('continuousStream51to60.mat');
        load('continuousStream1to50.mat');
        load('triggers_51to60.mat');
        load('triggers_1to50.mat');
        load('ChannelDataTimeStamps_51to60.mat');
        load('ChannelDataTimeStamps_1to50.mat');
    else
        % configure the segment of data to read:
        cfg_1to50=[];
        cfg_1to50.window=[55496  15008090].*samplerate; % data slicing w.r.t time [s] [5549.7, 1999894.1] / 2031874 / 5026385 / 10037200]
        cfg_51to60=[];
        cfg_51to60.window=[42200 3016750].*samplerate;  % 1st trigger @ [4220.1, 401491.9] 421506 / 102040.8 / 2218212]
        data_path_1to50 = 'C:/Users/nivsm/Documents/FinalProject/MEA_experiments/MultichannelDataManager/50_first_excites_approx20trials0001.h5';
        data_path_51to60 = 'C:/Users/nivsm/Documents/FinalProject/MEA_experiments/MultichannelDataManager/51to60_excites_approx20trials.h5';
        data_1to50 = McsHDF5.McsData(data_path_1to50);
        data_51to60 = McsHDF5.McsData(data_path_51to60);
        % triggers:                 1 X #excitations
        triggers_51to60 = data_51to60.Recording{1, 1}.EventStream{1, 1}.readPartialEventData(cfg_51to60).Events{1, 1};
        %save('triggers_51to60.mat', 'triggers_51to60');
        triggers_1to50 = data_1to50.Recording{1, 1}.EventStream{1, 1}.readPartialEventData(cfg_1to50).Events{1, 1};
        %save('triggers_1to50.mat', 'triggers_1to50');
        % continuous data stream:   #channels X #samples
        M_51to60 = data_51to60.Recording{1, 1}.AnalogStream{1, 2}.readPartialChannelData(cfg_51to60);
        M_51to60Cont =M_51to60.ChannelData;
        ChannelDataTimeStamps_51to60 = M_51to60.ChannelDataTimeStamps;
        save('ChannelDataTimeStamps_51to60.mat', 'ChannelDataTimeStamps_51to60');
        %save('continuousStream51to60.mat', 'M_51to60Cont')
        M_1to50 = data_1to50.Recording{1, 1}.AnalogStream{1, 2}.readPartialChannelData(cfg_1to50);
        M_1to50Cont =M_1to50.ChannelData;
        ChannelDataTimeStamps_1to50 = M_1to50.ChannelDataTimeStamps;
        save('ChannelDataTimeStamps_1to50.mat', 'ChannelDataTimeStamps_1to50');
        %save('continuousStream1to50.mat', 'M_1to50Cont')
    end

    % extract pulses, rearrange data:
    trig_indices_1to50 = trigger2triggerIndices(triggers_1to50, ChannelDataTimeStamps_1to50);
    trig_indices_51to60 = trigger2triggerIndices(triggers_51to60, ChannelDataTimeStamps_51to60);
    M_1to50 = continuousStream2triggered(M_1to50Cont, trig_indices_1to50, 50, num_trials);
    M_51to60 = continuousStream2triggered(M_51to60Cont, trig_indices_51to60, 10, num_trials);


    % Potential Matrix - coloumn # is the exciting electrode, row # is the
    % measuring electrode
    Pot_Mat = abs(horzcat(M_1to50, M_51to60));
    save('Pot_Mat_experiment_24_2_2021_X2.mat', 'Pot_Mat')
    return
    % TODO: verify setup of Pot_Mat is correct
    %for idx = 1:num_trials
    %    Pot_Mat = Pot_Mat + mM(:,60*(idx-1)+1:60*idx);
    %end
else
    Pot_Mat = load('Pot_Mat_experiment_24_2_2021.mat');
    Pot_Mat = Pot_Mat.Pot_Mat;
end

%% remove bad electrodes - in this case 15 is REF/GND thus it is not used
% use the mea60layout-> mea60Pin2index to enter electrode by index

bad_electrodes = [15];

for idx = sort(bad_electrodes,'descend')
    Pot_Mat(idx,:) = [];
    Pot_Mat(:,idx) = [];
end



%% defining hyperparameters
dX=         500; %distance between electrodes
alpha0=     1.5; %exponential decay constant
noiseStd=   0.1; 
len =       length(Pot_Mat); %number of electrodes

%% Regression analysis of potentials and distances to see if data conforms to model
load('mea60layout.mat');
%load('Solution_Coordinates.mat')
Coordinantes = mea60_coordinantes;
% TODO: eliminate Z axis in Coordinantes, ->near_sol accordingly
for idx = sort(bad_electrodes,'descend')
    Coordinantes(idx,:) =[];
end
r = dist(Coordinantes.');
v=Pot_Mat;
r(logical(eye(length(r)))) = 0;
v(logical(eye(length(Pot_Mat)))) = inf;

% model function: Vij = Aj*rij^(-alpha) + Const
modelfun = @(c,x)(c(1).*(x.^-c(2)))+c(3);
modelfit = cell(len,1);
reg_Coef = zeros(len,3);
R_sq = zeros(len,1);
for idx = 1:len
    x = r(:,idx).*10^-6;
    x(idx) = [];
    y = v(:,idx).*10^-12;
    y(idx) = [];
                        %initial values:[Aj,alpha,Const]
    modelfit{idx} = fitnlm(x,y,modelfun,[10^-08 1 10^-5]);
    reg_Coef(idx,:) = modelfit{idx}.Coefficients.Estimate';
    R_sq(idx) = modelfit{idx}.Rsquared.Adjusted;
    
    plot(x,modelfit{idx}.Fitted)
    
 end
residuals = zeros(len-1);
for idx = 1:len
    tab = table2array(modelfit{idx}.Variables);
    [~,xorder]=sort(tab(:,1));
    restab = table2array(modelfit{idx}.Residuals(:,1));
    residuals(idx,:) = restab(xorder);
    
end
figure
plot(mean(abs(residuals)))
title("residuals plot")


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
           res_sum(dom_size+1) = res(jdx);
       else
           sum_size(find(domain==range(jdx))) = sum_size(find(domain==range(jdx))) +1;
           res_sum(find(domain==range(jdx))) = res_sum(find(domain==range(jdx))) + abs(res(jdx));
       end
   end
end
res_avg = res_sum./sum_size;
[~,domain_order] = sort(domain);
figure
plot(domain(domain_order),res_avg(domain_order))
title("Residual vs electrode distance plot")


%%plotting a single regression
reg_num = 58;
figure
tab = table2array(modelfit{reg_num}.Variables);
[x,xorder]=sort(tab(:,1));
y = tab(xorder,2);
yhat = modelfit{reg_num}.Fitted(xorder);
plot(x,y,'.r')
hold on
plot(x,yhat)
title(['Regression of potential decay from electrode ' num2str(reg_num)])
xlabel('Distance [m]')
ylabel('Potential [V]')
%%
% We want to define constraints to the optimization problem:
% Coordinantes: physical contraints in terms of physical distances
% Elec_constants: ???
% alpha: must be larger than zero (otherwise exponent is not deminishing, as physically expected) ???
% We assume Aj and Ai should be similar, so our contraint is on A's standard deviation stdMax = 2

A =         [];
B =         [];
Aeq =       [];
beq =       [];
%           [Coordinantes                    ,     A_vect      ;alpha,   Vo      ,   ,  ]
lb =        [     zeros(len,3)               ,1e-5*ones(len,1) ; 0 , -1000^10^-6 , 0 , 0];
ub =        [0.050.*ones(len,2) zeros(len, 1),1e4*ones(len,1)  ; 4 , 1000^10^-6 , 0 , 0];
nonlcon =   [@(A) A_VectConstraint];
%reference for fmincon optimoptions
%https://www.mathworks.com/help/releases/R2020a/optim/ug/fmincon.html?browser=F1help#busog7r-options
options1 =  optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-20,'StepTolerance',1e-20,'ConstraintTolerance',1e-16,'MaxIterations',1e6,'MaxFunctionEvaluations',1e5,'Algorithm','sqp', 'PlotFcn', 'optimplotfval');
%options2 =  optimoptions('fmincon','Display','final-detailed','OptimalityTolerance',1e-20,'StepTolerance',1e-20,'ConstraintTolerance',1e-16,'MaxIterations',100000000,'MaxFunctionEvaluations',100000000);
%options3 =  optimoptions('fminunc','Display','final-detailed','OptimalityTolerance',1e-20,'MaxIterations',1000000,'MaxFunctionEvaluations',1000000);
%options4 =  optimset('Display','final-detailed','MaxIter',1000000,'MaxFunEvals',100000,'TolFun',1e-10);


%%
% We will define a mesh grid on the z=0 plane as the initial conditions of
% the optimization. If the number of electrodes is not a perfect square we
% will create a mesh grid that with a larger number of electrodes and then
% % eliminate excess electrodes at random.
% Mesh_Size = nearest_upper_square(len);
% [X,Y]=meshgrid(0:dX:(dX*sqrt(Mesh_Size)-1),0:dX:(dX*sqrt(Mesh_Size)-1));
% Mesh=[X(:) Y(:) zeros(Mesh_Size,1) ones(Mesh_Size,1)]; %the mesh guess for the A constants is the mean of the constants calculated in the regression
% while length(Mesh) > len
%     Mesh(round(randinterval(1,1,1,1,120)),:) = [];
% end
% Mesh(:,4) = mean(reg_Coef(:,1));
% Mesh(end+1,:) =         [mean(reg_Coef(:,2)) 0 0 0]; 
% Pot_Mat = abs(Pot_Mat - mean(reg_Coef(:,3))).*~eye(size(Pot_Mat));
% A_struct = load('As_solved_by_optimization');
% A_vect = A_struct.Newt.Vector;
v0 = Pot_Mat;
v0(logical(eye(len))) = inf;
v0 = min(min(v0))*10^-12-eps;
A_Vect = reg_Coef(:,1);
alpha = mean(reg_Coef(:,2));
near_sol_ = [Coordinantes.*10^-6 A_Vect ;alpha v0 0 0];
near_sol = near_sol_ + near_sol_.*abs([noiseStd*randn(len,3)  noiseStd*(randn(len,1)) ; 0 0 0 0]);
%% importing impedances from optimization and testing value at initial conditions with approximated impedances
%Con_As = load('Con_As.mat');
%GA_As = load('GA_As.mat');
%Newt_As = load('Newt_As.mat');
%PSO_As = load('PSO_As.mat');
%SQP_As = load('SQP_As.mat');

%Con_As.init_val = errorFun_3d([Coordinantes Con_As.Con.Vector; alpha v0 0 0])
%GA_As.init_val = errorFun_3d([Coordinantes GA_As.GA.Vector'; alpha v0 0 0])
%Newt_As.init_val = errorFun_3d([Coordinantes Newt_As.Newt.Vector; alpha v0 0 0])
%PSO_As.init_val = errorFun_3d([Coordinantes PSO_As.PSO.Vector'; alpha v0 0 0])
%SQP_As.init_val = errorFun_3d([Coordinantes SQP_As.SQP.Vector; alpha v0 0 0])
%regression_As.init_val = errorFun_3d([Coordinantes A_Vect; alpha v0 0 0])
x0 = near_sol;
Sol_Vector = fmincon(@(x) errorFun_3d(x),x0,A,B,Aeq,beq,lb,ub,nonlcon,options1);

return




%%
n=100;
nvars = 4*(len+1);
x0 =                zeros(len+1,4,n);
% Con.Vector =        zeros(len+1,4,n);
% Con.Time =          zeros(n,1);
% Con.Value =         zeros(n,1);
% Con.Dist =          zeros(n,1); 
% Con.Var_As =        zeros(n,1);
% Con.Mean_As =       zeros(n,1);

SQP.Vector =   zeros(len+1,4,n);
SQP.Time =     zeros(n,1);
SQP.Value =         zeros(n,1);
SQP.Dist =          zeros(n,1); 
SQP.Var_As =        zeros(n,1);
SQP.Mean_As =       zeros(n,1);

% Newt.Vector =   zeros(len+1,4,n);
% Newt.Time =     zeros(n,1);
% Newt.Value =         zeros(n,1);
% Newt.Dist =          zeros(n,1); 
% Newt.Var_As =        zeros(n,1);
% Newt.Mean_As =       zeros(n,1);
% 
% PSO.Vector =   zeros(len+1,4,n);
% PSO.Time =     zeros(n,1);
% PSO.Value =         zeros(n,1);
% PSO.Dist =          zeros(n,1); 
% PSO.Var_As =        zeros(n,1);
% PSO.Mean_As =       zeros(n,1);
% 
% GA.Vector =   zeros(len+1,4,n);
% GA.Time =     zeros(n,1);
% GA.Value =         zeros(n,1);
% GA.Dist =          zeros(n,1); 
% GA.Var_As =        zeros(n,1);
% GA.Mean_As =       zeros(n,1);


for idx = 1:n

    x0(:,:,idx) = near_sol;
%     x0(:,:,idx) = near_sol;
    %Con(InteriorPoint)
%     tic
%     Con.Vector(:,:,idx) = fmincon(@(x) errorFun_3d(x),x0(:,:,idx),A,B,Aeq,beq,lb,ub,nonlcon,options2);
%     Con.Time(idx) = toc;
%     Con.Value(idx) = errorFun_3d(Con.Vector(:,:,idx));
%     Con.Dist(idx) = sum(sum(abs(dist(Con.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
%     Con.Var_As(idx) = var(Con.Vector(1:len,4,idx));
%     Con.Mean_As(idx) = mean(Con.Vector(1:len,4,idx));
%     
    %SQP
    tic
    SQP.Vector(:,:,idx) = fmincon(@(x) errorFun_3d(x),x0(:,:,idx),A,B,Aeq,beq,lb,ub,nonlcon,options1);
    SQP.Time(idx) = toc;
    SQP.Value(idx) = errorFun_3d(SQP.Vector(:,:,idx));
    SQP.Dist(idx) = sum(sum(abs(dist(SQP.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
    SQP.Var_As(idx) = var(SQP.Vector(1:len,4,idx));
    SQP.Mean_As(idx) = mean(SQP.Vector(1:len,4,idx));
    
    %Newt
%     tic
%     Newt.Vector(:,:,idx) = fminunc(@(x) errorFun_3d(x),x0(:,:,idx),options3);
%     Newt.Time(idx) = toc;
%     Newt.Value(idx) = errorFun_3d(Newt.Vector(:,:,idx));
%     Newt.Dist(idx) = sum(sum(abs(dist(Newt.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
%     Newt.Var_As(idx) = var(Newt.Vector(1:len,4,idx));
%     Newt.Mean_As(idx) = mean(Newt.Vector(1:len,4,idx));
     
%     %PSO
%     tic
%     PSO.Vector(:,:,idx) = vect2mat(particleswarm(@(x) errorFun_vector(x),nvars,mat2vect(lb),mat2vect(ub)));
%     PSO.Time(idx) = toc;
%     PSO.Value(idx) = errorFun_vector(mat2vect(PSO.Vector(:,:,idx)));
%     PSO.Dist(idx) = sum(sum(abs(dist(PSO.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
%     PSO.Var_As(idx) = var(PSO.Vector(1:len,4,idx));
%     PSO.Mean_As(idx) = mean(PSO.Vector(1:len,4,idx));
%      
%     %GA
%     tic
%     GA.Vector(:,:,idx) = vect2mat(ga(@(x) errorFun_vector(x),nvars,A,B,Aeq,beq,mat2vect(lb),mat2vect(ub)));
%     GA.Time(idx) = toc;
%     GA.Value(idx) = errorFun_vector(mat2vect(GA.Vector(:,:,idx)));
%     GA.Dist(idx) = sum(sum(abs(dist(GA.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
%     GA.Var_As(idx) = var(GA.Vector(1:len,4,idx));
%     GA.Mean_As(idx) = mean(GA.Vector(1:len,4,idx));
%     
     
    
    disp(['Now finishing iteration number:' num2str(idx)]) 
    
    
    
    
end
% [~,Con.best_iter] = min(Con.Value);
% [~,SQP.best_iter] = min(SQP.Value);
% [~,Newtbest_iter] = min(Newt.Value);
% [~,PSO.best_iter] = min(PSO.Value);
% [~,GA.best_iter] = min(GA.Value);

% 
% Best_Ans_Con = Con_Vector_iter(:,:,best_iter_Con);
% Best_Ans_SQP = SQP_Vector_iter(:,:,best_iter_SQP);
% Best_Ans_Newt = Newt_Vector_iter(:,:,best_iter_Newt);
% Best_Ans_PSO = PSO_Vector_iter(:,:,best_iter);
% Best_Ans_GA = GA_Vector_iter(:,:,best_iter_GA);

% 
% plot3(Best_Ans_Con(:,1),Best_Ans_Con(:,2),Best_Ans_Con(:,3),'.')
% plot3(Best_Ans_SQP(:,1),Best_Ans_SQP(:,2),Best_Ans_SQP(:,3),'.')
% plot3(Best_Ans_Newt(:,1),Best_Ans_Newt(:,2),Best_Ans_Newt(:,3),'.')
% plot3(Best_Ans_PSO(:,1),Best_Ans_PSO(:,2),Best_Ans_PSO(:,3),'.')
% plot3(Best_Ans_GA(:,1),Best_Ans_GA(:,2),Best_Ans_GA(:,3),'.')



% plot(1:length(Con_Value),sort(Con_Value))
% plot(1:n,Mean_As_Con)
% plot3(Con_Vector_iter(:,1,m),Con_Vector_iter(:,2,m),Con_Vector_iter(:,3,m),'.')

for idx = 1:10
    Near_Sol.f_x0(idx) = errorFun_3d(Near_Sol.x0(:,:,idx));
    
    
end