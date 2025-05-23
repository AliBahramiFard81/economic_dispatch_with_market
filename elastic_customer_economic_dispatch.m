%% by Ali Bahrami Fard
% alibt1313@gmail.com
% ED using optimization functions
clc;
clear;
close all;
format short g;

%% -------------------- INPUTS --------------------
% columns => 1.A, 2.B, 3.C, 4.Pmin(MW), 5.Pmax(MW), 6.Cost
% F(Pg) = A + B * Pg + C * Pg ^ 2
gen_inputs = [500, 5.3, 0.004, 50, 200, 1;
          400, 5.5, 0.006, 50, 150, 1;
          200, 4.8, 0.009, 30, 100, 1];

demand_inputs = [0, 120, -0.3, 0, 100, 1;
                 0, 100, -0.2, 0, 80, 1;
                 0, 95, -0.15, 0, 120, 1];
%{
inputs = [510, 7.2, 0.00142, 150, 600, 1.1;
          310, 7.85, 0.00194, 100, 400, 1;
          78, 7.97, 0.00482, 50, 200, 1];
inputs = [20, 0, 0, 0, 20, 1;
          50, 20, 0, 0, 50, 1;
          100, 30, 0, 0, 100, 1];

%}

% variables to solve
num_vars = length(gen_inputs(:,1)) + length(demand_inputs(:,1)) + 1;
syms Vars [num_vars, 1];

% Load(MW)
pload = 300;

%% -------------------- CALCULATING A,B,C GENERATION COEFFICIENTS --------------------
% calculating A coefficient
A_Co_Gen = zeros(length(gen_inputs(:,1)),1);
for i=1:length(gen_inputs(:,1))
    A_Co_Gen(i,1) = gen_inputs(i,1) * gen_inputs(i, 6);
end

% calculating B coefficient
B_Co_Gen = zeros(length(gen_inputs(:,2)),1);
for i=1:length(gen_inputs(:,2))
    B_Co_Gen(i,1) = gen_inputs(i,2) * gen_inputs(i, 6);
end

% calculating C coefficient
C_Co_Gen = zeros(length(gen_inputs(:,3)),1);
for i=1:length(gen_inputs(:,3))
    C_Co_Gen(i,1) = gen_inputs(i,3) * gen_inputs(i, 6);
end

%% -------------------- CALCULATING A,B,C DEMAND COEFFICIENTS --------------------
% calculating A coefficient
A_Co_D = zeros(length(demand_inputs(:,1)),1);
for i=1:length(demand_inputs(:,1))
    A_Co_D(i,1) = demand_inputs(i,1) * demand_inputs(i, 6);
end

% calculating B coefficient
B_Co_D = zeros(length(demand_inputs(:,2)),1);
for i=1:length(demand_inputs(:,2))
    B_Co_D(i,1) = demand_inputs(i,2) * demand_inputs(i, 6);
end

% calculating C coefficient
C_Co_D = zeros(length(demand_inputs(:,3)),1);
for i=1:length(demand_inputs(:,3))
    C_Co_D(i,1) = demand_inputs(i,3) * demand_inputs(i, 6);
end
%% -------------------- CONSTRUCTING F --------------------
F_gen = Vars(length(gen_inputs(:,1)));
F_demand = Vars( length(demand_inputs(:,1)));

for i=1:length(gen_inputs(:,1))
    F_gen(i,1) = A_Co_Gen(i,1) + B_Co_Gen(i,1) * Vars(i) + C_Co_Gen(i,1) * Vars(i) ^ 2;
end

for i=1:length(demand_inputs(:,1))
    F_demand(i,1) = A_Co_D(i,1) + B_Co_D(i,1) * Vars(i + length(gen_inputs(:,1))) + C_Co_D(i,1) * Vars(i + length(gen_inputs(:,1))) ^ 2;
end
%% -------------------- SOLVING USING OPTIMIZATION --------------------
% number of thermal units
G = length(F_gen);
D = length(F_demand);

% creating the main variables
Pi = optimvar('Pi', G);
Di = optimvar('Di', D);
Lambda = optimvar('L', 1);

C1_Gen = zeros(G, G);
for i=1:length(C_Co_Gen)
    C1_Gen(i,i) = C_Co_Gen(i);
end

C1_D = zeros(D, D);
for i=1:length(C_Co_D)
    C1_D(i,i) = C_Co_D(i);
end

objective_fun = (A_Co_Gen' * ones(G,1) + B_Co_Gen' * Pi + Pi' * C1_Gen * Pi) - (A_Co_D' * ones(D,1) + B_Co_D' * Di + Di' * C1_D * Di);
sum_F = sum(objective_fun);
problem = optimproblem("Objective", sum_F);

% power min and max constrains
problem.Constraints.gen_min = Pi >= gen_inputs(:,4);
problem.Constraints.gen_max = Pi <= gen_inputs(:,5);

% demnad min and max constrains
problem.Constraints.demand_min = Di >= demand_inputs(:,4);
problem.Constraints.demand_max = Di <= demand_inputs(:,5);

% load generation equality constrain
problem.Constraints.pd = sum(Pi) == sum(Di);

options = optimoptions('quadprog');
[sol,TotalCost,exitflag,output] = solve(problem, 'Options', options);

%% -------------------- CAACULATING LAMBDA --------------------
% selecting the unit that is not at its min or max power generation
% check the KKT rules for more informations
selected_p = 0;
index = 0;
for i=1:length(sol.Pi)
    index = index + 1;
    if sol.Pi(i) ~= gen_inputs(i,4) && sol.Pi(i) ~= gen_inputs(i,5)
       selected_p = sol.Pi(i);
       break;
    end
end

disp(selected_p)
disp(index)
lambda = subs(diff(F_gen(index)), selected_p);

%% -------------------- RESULTS --------------------    
clc;
disp("Unit Demand Results : ")
fprintf("\tUnit    \tDemand(MW)\n")
for i=1:length(sol.Di)
    if i < 10
        if sol.Di(i) ~= 0
            fprintf("\t%d  \t\t\t%.4f\n",[i,sol.Di(i)]);
        end
    else
        if sol.Di(i) ~= 0
            fprintf("\t%d  \t\t%.4f\n",[i,sol.Di(i)]);
        end
    end
end
fprintf("Total Demand    \t%.3f\n\n",sum(sol.Di))

disp("Unit Power Results : ")
fprintf("\tUnit    \tPower(MW)\n")
for i=1:length(sol.Pi)
    if i < 10
        if sol.Pi(i) ~= 0
            fprintf("\t%d  \t\t\t%.4f\n",[i,sol.Pi(i)]);
        end
    else
        if sol.Pi(i) ~= 0
            fprintf("\t%d  \t\t%.4f\n",[i,sol.Pi(i)]);
        end
    end
end
fprintf("Total Generation    \t%.3f\n\n",sum(sol.Pi))

fprintf("\n\n")

fprintf("Total cost    \t%.3f\n",TotalCost)

fprintf("\n\n")

fprintf("Lambda    \t\t%.3f\n",lambda)
