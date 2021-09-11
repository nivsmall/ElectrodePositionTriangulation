global Pot_Mat len A_Vect ACurr alpha_mu v0 coord_coeff A_coeff alpha_coeff Coordinates Pot_Mat_input chooser oneAtATime distances res distFilterMask CSDTrain CSD;

%%choosing coeeficients to train:
coord_coeff=true;
alpha_coeff=false;
A_coeff=true;

simulated=false;
oneAtATime=false;
lossLandscape=false;
CSDTrain = false;

pltTitle = 'protocol: decay->coupling&coordinates(STD=1500) ';
lossPltTitle = "Loss Landscape - MEA Data | decay->coupling&coordinates(STD=1500)";


%% defining hyperparameters
dX=         500; %distance between electrodes
alpha0=     -1.1413; %exponential decay constant
noiseStd=  1500;
V0 = 0;%72.2*10^6;
len = 60;

load('mea60layout.mat');
Coordinates = mea60_coordinantes;
A_Vect = ones(len, 1);  % abs(0.25.*randn(len, 1));
r = dist(Coordinates.');
excited =   100*ones(len); %[mv]%100*10^-3;
excitedMeasured_ratio = (A_Vect.*r.^(alpha0)+V0./excited);
Expected_Pot = excited.*excitedMeasured_ratio;
Pot_Mat_ = load('Pot_Mat_experiment_24_2_2021_X2.mat');
Pot_MatMea = Pot_Mat_.Pot_Mat;
meas_theory_ratio = Pot_MatMea./excitedMeasured_ratio;

if simulated
    Pot_Mat = excitedMeasured_ratio;
else
    Pot_Mat = Pot_MatMea;
end

%% remove bad electrodes - in this case 15 is REF/GND thus it is not used
% use the mea60layout-> mea60Pin2index to enter electrode by index
if false
    bad_electrodes = [15];
else
    bad_electrodes = [15, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60];
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

%% importing impedances, alpha from optimization and testing value at initial conditions with approximated impedances
%[alpha, distances, res] = regression_analysis(Coordinates, Pot_MatMea, len, false);
%[alpha, distances, res] = regression_analysis(Coordinates, CSD, len, true);

alpha_mu = alpha0;
A_Vect = load("A_Vect.mat").A_Vect;
v0 = V0;
r = dist(Coordinates');

correlation = checkElectrodesCoupling(Pot_Mat, A_Vect);
%addNoiseManually(r)
distFilterMask = zeros(len);
%distFilterMask = filterByDistances(r, 520, 10000);
%% boundaries for optimization:
if coord_coeff && A_coeff && alpha_coeff
    near_sol = [Coordinates(:, 1:2) A_Vect; alpha_mu 0 0];
elseif coord_coeff && A_coeff && ~alpha_coeff
    near_sol = [Coordinates(:, 1:2)  A_Vect];
elseif coord_coeff && ~A_coeff && ~alpha_coeff
    near_sol = Coordinates(:, 1:2);
elseif ~coord_coeff && ~A_coeff && alpha_coeff
    near_sol = alpha_mu;
elseif ~coord_coeff && A_coeff && alpha_coeff
    near_sol = [A_Vect; alpha_mu];
elseif ~coord_coeff && A_coeff && ~alpha_coeff
    near_sol = A_Vect;
end


Coord_lb = Coordinates(:, 1:2) - 3500.*ones(len,2);
Coord_ub = Coordinates(:, 1:2) + 3500.*ones(len,2);
A_Vect_lb = zeros(len, 1);
A_Vect_ub = 20000*ones(len, 1);
alpha_lb = 20*alpha_mu;
alpha_ub = 0.005*alpha_mu;
A =         [];
B =         [];
Aeq =       [];
beq =       [];
if coord_coeff && A_coeff && alpha_coeff
    lb = [Coord_lb A_Vect_lb; alpha_lb 0 0];
    ub = [Coord_ub A_Vect_ub; alpha_ub 0 0];
elseif coord_coeff && A_coeff && ~alpha_coeff
    lb = [Coord_lb A_Vect_lb];
    ub = [Coord_ub A_Vect_ub];
elseif coord_coeff && ~A_coeff && ~alpha_coeff
    lb =  Coord_lb;
    ub =  Coord_ub;
elseif alpha_coeff && ~coord_coeff && ~A_coeff
    lb = alpha_lb;
    ub = alpha_ub;
elseif alpha_coeff && ~coord_coeff && A_coeff
    lb = [A_Vect_lb; alpha_lb];
    ub = [A_Vect_ub; alpha_ub];
elseif ~coord_coeff && A_coeff && ~alpha_coeff
    lb = A_Vect_lb;
    ub = A_Vect_ub;
end


nonlcon =   [];
options1 =  optimoptions('fmincon','Display','final-detailed','OptimalityTolerance',1e-7,'StepTolerance',1e-8,'ConstraintTolerance',1e-16,'MaxIterations',10000000,'MaxFunctionEvaluations',1000000,'Algorithm','interior-point','PlotFcn', {'optimplotfval','optimplotx','optimplotfirstorderopt'});
options2 =  optimoptions('fmincon','Display','final-detailed','OptimalityTolerance',1e-20,'StepTolerance',1e-20,'ConstraintTolerance',1e-16,'MaxIterations',100000000,'MaxFunctionEvaluations',100000000);
options3 =  optimoptions('fminunc','Display','final-detailed','OptimalityTolerance',1e-20,'MaxIterations',1000000,'MaxFunctionEvaluations',1000000);
options4 =  optimset('Display','final-detailed','MaxIter',1000000,'MaxFunEvals',100000,'TolFun',1e-10);



%% 
n=5;
if coord_coeff && A_coeff && alpha_coeff
    x0 = zeros(len+1,3,n);
elseif coord_coeff && A_coeff && ~alpha_coeff
    x0 = zeros(len,3,n);
elseif coord_coeff && ~A_coeff && ~alpha_coeff
    x0 = zeros(len,2,n);
elseif ~coord_coeff && ~A_coeff && alpha_coeff
    x0 = zeros(1, n);
elseif ~coord_coeff && A_coeff && alpha_coeff
    x0 = zeros(len+1, 1, n);
elseif ~coord_coeff && A_coeff && ~alpha_coeff
    x0=zeros(len, 1, n);
end


SQP.Vector =   x0;
SQP.Time =     zeros(n,1);
SQP.Value =    zeros(n,1);
SQP.Dist =     zeros(n,1); 
SQP.Var_As =   zeros(n,1);
SQP.Mean_As =  zeros(n,1);

for idx = 1:n
    if coord_coeff && A_coeff && alpha_coeff
        x0(:,:,idx) = near_sol+[noiseStd*randn(len,2) noiseStd/100*randn(len, 1); noiseStd/100*randn(1) 0 0]; %noiseStd/100*randn(len, 1)
    elseif coord_coeff && A_coeff && ~alpha_coeff
        x0(:,:,idx) = near_sol+noiseStd*[randn(len,2) noiseStd/100*randn(len, 1)];
    elseif coord_coeff && ~A_coeff && ~alpha_coeff
        x0(:,:,idx) = near_sol+noiseStd*randn(len,2);
        %Coordinates(:, 1:2) = Coordinates(:, 1:2)+noiseStd.*randn(len, 2);
    elseif alpha_coeff && ~coord_coeff && ~A_coeff
        x0 = near_sol;%+noiseStd;
        Coordinates(:, 1:2) = Coordinates(:, 1:2)+noiseStd.*randn(len, 2);
    elseif alpha_coeff && ~coord_coeff && A_coeff
        x0 = near_sol;%+noiseStd.*randn(len+1,1);
    elseif ~coord_coeff && A_coeff && ~alpha_coeff
        x0 = near_sol;
    end
    


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
            x_optimal = fmincon(@ObjectiveFunctionOneAtaTime,x0Curr,A,B,Aeq,beq,lbCurr,ubCurr,nonlcon,options1);

            plot_single_train_results(Coordinates_init_copy, CoordCurr, x_optimal, lbCurr, ubCurr);

            %Coordinates_init_copy = Coordinates_init;
        end
        
    else
            if lossLandscape
                plotLossLandscape(lossPltTitle);
                continue
            end
            %SQP
            tic
            SQP.Vector(:,:,idx) = fmincon(@errorFun_2d,x0(:,:,idx),A,B,Aeq,beq,lb,ub,nonlcon,options1)
            %SQP.Time(idx) = toc;
            %SQP.Value(idx) = errorFun_2d(SQP.Vector(:,:,idx));
            %SQP.Dist(idx) = sum(sum(abs(dist(SQP.Vector(1:len,1:2,idx)')-r)))/(-len+len^2);
            if coord_coeff
                plot_results(Coordinates, x0(:,:,idx), SQP.Vector(:, :, idx), pltTitle)
            end
            if A_coeff && ~coord_coeff
                A_Vect = SQP.Vector(:,:,idx);
                save("A_Vect.mat", "A_Vect");
            end
%             alpha_coeff = false;
%             A_coeff = false;
%             coord_coeff = true;
%             plotLossLandscape(lossPltTitle);
    end
    disp(['Now finishing iteration number:' num2str(idx)])
end

if coord_coeff
    finalResult = zeros(len, 2);
    for coordIdx = 1:len
        Ymle = mle(reshape(SQP.Vector(coordIdx, 1, :), [n, 1]));
        Xmle = mle(reshape(SQP.Vector(coordIdx, 2, :), [n, 1]));
        finalResult(coordIdx, 1) = Ymle(1);
        finalResult(coordIdx, 2) = Xmle(1);
    end
    pltTitle = strcat(pltTitle, '  -  final');
    plot_results(Coordinates, x0(:,:,idx), finalResult, pltTitle);
end

