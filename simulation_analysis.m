global Pot_Mat len A_Vect ACurr alpha_mu x0_init v0 coord_coeff A_coeff alpha_coeff Coordinates Pot_Mat_input chooser oneAtATime distances res distFilterMask CSDTrain harmonicSignal;

%%choosing coeeficients to train %%
% coord_coeff=true;
% alpha_coeff=false;
% A_coeff=false;

%% choosing experimental setup %%
simulated=true;
oneAtATime=false;
lossLandscape=false;
CSDTrain = false;
fiveByfive=true;
harmonicSignal=false;

%% defining hyperparameters %%
dX=         500; %distance between electrodes
alpha0=     -1.1539; %exponential decay constant
%alpha0=     -1.0769; %exponential decay constant
noiseStd=  150;
noiseStsStr = "(STD = " + string(noiseStd) + ")";
V0 = 0;
len = 60;

%%
%pltTitle = {'protocol: ', noiseStsStr, ''};
pltTitle = {'Corrected Positions of Electrodes'};
lossPltTitle = "Loss Landscape - MEA Data | " + noiseStsStr;

%% loading the data %%
load('mea60layout.mat');
Coordinates = mea60_coordinantes;
A_Vect = ones(len, 1);  % abs(0.25.*randn(len, 1));
r = dist(Coordinates.');
excited =   100*ones(len); %[mv]%100*10^-3;
excitedMeasured_ratio = (A_Vect.*r.^(alpha0)+V0./excited);
Expected_Pot = excited.*excitedMeasured_ratio;
if harmonicSignal
    Pot_Mat_ = load('Pot_Mat_harmonicMeas.mat');
else
    Pot_Mat_ = load('Pot_Mat_experiment_24_2_2021_X2.mat');
end
Pot_MatMea = Pot_Mat_.Pot_Mat;
meas_theory_ratio = Pot_MatMea./excitedMeasured_ratio;

if simulated
    Pot_Mat = excitedMeasured_ratio*10^12;
else
    Pot_Mat = Pot_MatMea;
end


%% remove bad electrodes - in this case 15 is REF/GND thus it is not used %%
% use the mea60layout-> mea60Pin2index to enter electrode by index %
if fiveByfive
    participating_electrodes = [26, 29, 32, 35, 37, 27, 28, 33, 34, 38, 25, 30, 31, 36, 40, 20, 22, 39, 41, 42, 17, 18, 43, 44, 45];%, ...
%          24, 23, 21, 19, 16, ...
%          14, 13, 48, 47, 46, ...
%          12, 11, 9, 52, 50, 49, ...
%          10, 6, 1, 60, 55, 51, ...
%          8, 4, 3, 58, 57, 53, ...
%          7, 5, 2, 59, 56, 54];
    bad_electrodes = nan(1,len-length(participating_electrodes));
    cnt=1;
    for idx = 1:len
        if find(participating_electrodes == idx)
            continue;
        else
            bad_electrodes(cnt)=idx;
            cnt=cnt+1;
        end
    end
else
    if true
        bad_electrodes = [15];
    else
        bad_electrodes = [15, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60];
    end
end

for idx = sort(bad_electrodes,'descend')
    Pot_Mat(idx,:) = [];
    Pot_Mat(:,idx) = [];
    Pot_MatMea(idx,:) = [];
    Pot_MatMea(:,idx) = [];
    Expected_Pot(idx,:) = [];
    Expected_Pot(:,idx) = [];
    excitedMeasured_ratio(idx,:) = [];
    excitedMeasured_ratio(:,idx) = [];
    Coordinates(idx, :) = [];
    A_Vect(idx) = [];
end

len = length(Pot_Mat);

%% Current Source Density (CSD):
if CSDTrain
    CSD = potentialToCSD(Pot_Mat, len, Coordinates);
    %Pot_Mat = CSD;
end

%% importing impedances, alpha from optimization and testing value at initial conditions with approximated impedances %%
[alpha, distances, res] = regression_analysis(Coordinates, Pot_MatMea, len, false);
[alpha, distances, res] = regression_analysis(Coordinates, CSD, len, true);
                     
alpha_mu = alpha0;
%A_Vect = load("A_VectPulse.mat").A_Vect;
v0 = V0;
r = dist(Coordinates');

%correlation = checkElectrodesCoupling(Pot_Mat, A_Vect);
%addNoiseManually(r)
distFilterMask = zeros(len);
%distFilterMask = filterByDistances(r, 520, 10000);
%% boundaries for optimization %%


for i=0:0
% if coord_coeff && A_coeff && alpha_coeff
%     near_sol = [Coordinates(:, 1:2) A_Vect; alpha_mu 0 0];
% elseif coord_coeff && A_coeff && ~alpha_coeff
%     near_sol = [Coordinates(:, 1:2)  A_Vect];
% elseif coord_coeff && ~A_coeff && ~alpha_coeff
%     near_sol = Coordinates(:, 1:2);
% elseif ~coord_coeff && ~A_coeff && alpha_coeff
%     near_sol = alpha_mu;
% elseif ~coord_coeff && A_coeff && alpha_coeff
%     near_sol = [A_Vect; alpha_mu];
% elseif ~coord_coeff && A_coeff && ~alpha_coeff
%     near_sol = A_Vect;
% end
    
% Coord_lb = Coordinates(:, 1:2) - 3500.*ones(len,2);
% Coord_ub = Coordinates(:, 1:2) + 3500.*ones(len,2);
% A_Vect_lb = zeros(len, 1);
% A_Vect_ub = 20000*ones(len, 1);
% alpha_lb = 20*alpha_mu;
% alpha_ub = 0.005*alpha_mu;
% A =         [];
% B =         [];
% Aeq =       [];
% beq =       [];
% if coord_coeff && A_coeff && alpha_coeff
%     lb = [Coord_lb A_Vect_lb; alpha_lb 0 0];
%     ub = [Coord_ub A_Vect_ub; alpha_ub 0 0];
% elseif coord_coeff && A_coeff && ~alpha_coeff
%     lb = [Coord_lb A_Vect_lb];
%     ub = [Coord_ub A_Vect_ub];
% elseif coord_coeff && ~A_coeff && ~alpha_coeff
%     lb =  Coord_lb;
%     ub =  Coord_ub;
% elseif alpha_coeff && ~coord_coeff && ~A_coeff
%     lb = alpha_lb;
%     ub = alpha_ub;
% elseif alpha_coeff && ~coord_coeff && A_coeff
%     lb = [A_Vect_lb; alpha_lb];
%     ub = [A_Vect_ub; alpha_ub];
% elseif ~coord_coeff && A_coeff && ~alpha_coeff
%     lb = A_Vect_lb;
%     ub = A_Vect_ub;
% end
end

nonlcon =   [];
options1 =  optimoptions('fmincon','Display','final-detailed','OptimalityTolerance',1e-6,'StepTolerance',1e-8,'ConstraintTolerance',1e-16,'MaxIterations',10000000,'MaxFunctionEvaluations',1000000,'Algorithm','interior-point','PlotFcn', {'optimplotfval','optimplotx','optimplotfirstorderopt'});
options2 =  optimoptions('fmincon','Display','final-detailed','OptimalityTolerance',1e-20,'StepTolerance',1e-20,'ConstraintTolerance',1e-16,'MaxIterations',100000000,'MaxFunctionEvaluations',100000000);
options3 =  optimoptions('fminunc','Display','final-detailed','OptimalityTolerance',1e-20,'MaxIterations',1000000,'MaxFunctionEvaluations',1000000);
options4 =  optimset('Display','final-detailed','MaxIter',1000000,'MaxFunEvals',100000,'TolFun',1e-10);



%% session run %%
protocol = {'alpha','coord','A', 'coord'};%,'alpha','coord', 'alpha','coord', 'alpha','coord'};
near_sol = [Coordinates(:, 1:2) A_Vect; alpha_mu 0 0];
x0 = near_sol + [noiseStd*randn(len,2) zeros(len, 1); 0 0 0];
x0_init = x0;

n=length(protocol);
for i=0:0
% if coord_coeff && A_coeff && alpha_coeff
%     x0 = zeros(len+1,3,n);
% elseif coord_coeff && A_coeff && ~alpha_coeff
%     x0 = zeros(len,3,n);
% elseif coord_coeff && ~A_coeff && ~alpha_coeff
%     x0 = zeros(len,2,n);
% elseif ~coord_coeff && ~A_coeff && alpha_coeff
%     x0 = zeros(1, n);
% elseif ~coord_coeff && A_coeff && alpha_coeff
%     x0 = zeros(len+1, 1, n);
% elseif ~coord_coeff && A_coeff && ~alpha_coeff
%     x0=zeros(len, 1, n);
% end
%x0 = zeros(len+1,3,n);
end

Session.Vector =   x0;
Session.Time =     zeros(n,1);
Session.Value =    zeros(n,1);
Session.Dist =     zeros(n,1); 
Session.Var_As =   zeros(n,1);
Session.Mean_As =  zeros(n,1);


for idx = 1:n
    
    if oneAtATime
        Coordinates_init_copy = Coordinates;
        Coordinates = x0(1:len,1:2,idx);
        for chooser = 1:len
            
            Pot_Mat_input = Pot_Mat(chooser, :);
            Pot_Mat_input(chooser) = [];
            
            ACurr = A_Vect(chooser);
            CoordCurr = x0(chooser, 1:2, idx);
            
            if A_coeff
                lbCurr = lb(chooser, 1:3);
                ubCurr = ub(chooser, 1:3);
                x0Curr = [CoordCurr ACurr];
            else            
                lbCurr = lb(chooser, 1:2);
                ubCurr = ub(chooser, 1:2);
                x0Curr = CoordCurr;
            end
            if lossLandscape
                plotLossLandscape();
                continue
            end
            x_optimal = fmincon(@ObjectiveFunctionOneAtaTime,x0Curr,A,B,Aeq,Beq,lbCurr,ubCurr,nonlcon,options1);

            plot_single_train_results(Coordinates_init_copy, CoordCurr, x_optimal, lbCurr, ubCurr);

            %Coordinates_init_copy = Coordinates_init;
        end
        
    else
            if lossLandscape
                plotLossLandscape(lossPltTitle);
                continue
            end
            
            [coord_coeff, A_coeff, alpha_coeff, STR] = coefficientToTrain(protocol(idx));
            [lb, ub, A, B, Aeq, Beq] = initBoundaries(coord_coeff, A_coeff, alpha_coeff);
            %pltTitle(1,1) = {pltTitle(1,1) + STR};
            
            %SQP
            tic
            
            if coord_coeff
                Session.Vector(1:len,1:2) = fmincon(@errorFun_2d,x0(1:len,1:2),A,B,Aeq,Beq,lb,ub,nonlcon,options1);
            elseif A_coeff
                Session.Vector(1:len,3) = fmincon(@errorFun_2d,x0(1:len,3),A,B,Aeq,Beq,lb,ub,nonlcon,options1);
            elseif alpha_coeff
                Session.Vector(len+1,1) = fmincon(@errorFun_2d,x0(len+1,1),A,B,Aeq,Beq,lb,ub,nonlcon,options1);
                Session.Vector(len+1)
            Session.Time(idx) = toc;
            end
            %SQP.Value(idx) = errorFun_2d(SQP.Vector(:,:,idx));
            if coord_coeff
                Session.Dist = mean(abs(dist(Session.Vector(1:len,1:2)')-r), 'all');
                Session.Dist
                ME = mean(abs(Session.Vector(1:len,1:2) - Coordinates(1:len, 1:2)), 'all');
                %pltTitle(1, 3) = {'Residuals Mean: ' + string(SQP.Dist)};
                plot_results(Coordinates, x0(1:len,1:2), Session.Vector(1:len,1:2), pltTitle)
            end
            if A_coeff && ~coord_coeff
                A_Vect = Session.Vector(1:len,3);
                %save("A_Vect.mat", "A_Vect");
            end
            
            x0_init = Session.Vector; %x0;
            
            
%             alpha_coeff = false;
%             A_coeff = false;
%             coord_coeff = true;
%             plotLossLandscape(lossPltTitle);
    end
    disp(['Now finishing iteration number:' num2str(idx)])
end
return


