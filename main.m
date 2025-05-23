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
inputs = [500, 5.3, 0.004, 50, 200, 1;
          400, 5.5, 0.006, 50, 150, 1;
          200, 4.8, 0.009, 30, 100, 1];

%{
inputs = [510, 7.2, 0.00142, 150, 600, 1.1;
          310, 7.85, 0.00194, 100, 400, 1;
          78, 7.97, 0.00482, 50, 200, 1];
inputs = [20, 0, 0, 0, 20, 1;
          50, 20, 0, 0, 50, 1;
          100, 30, 0, 0, 100, 1];

%}

% variables to solve
num_vars = length(inputs(:,1)) + 1;
syms Vars [num_vars, 1];

% Load(MW)
pload = 300;

%% -------------------- CALCULATING A,B,C COEFFICIENTS --------------------
% calculating A coefficient
A_Co = zeros(length(inputs(:,1)),1);
for i=1:length(inputs(:,1))
    A_Co(i,1) = inputs(i,1) * inputs(i, 6);
end

% calculating B coefficient
B_Co = zeros(length(inputs(:,2)),1);
for i=1:length(inputs(:,2))
    B_Co(i,1) = inputs(i,2) * inputs(i, 6);
end

% calculating C coefficient
C_Co = zeros(length(inputs(:,3)),1);
for i=1:length(inputs(:,3))
    C_Co(i,1) = inputs(i,3) * inputs(i, 6);
end

%% -------------------- CONSTRUCTING F --------------------
F = Vars(length(inputs(:,1)));
for i=1:length(inputs(:,1))
    F(i,1) = A_Co(i,1) + B_Co(i,1) * Vars(i) + C_Co(i,1) * Vars(i) ^ 2;
end

%% -------------------- SOLVING USING OPTIMIZATION --------------------
% number of thermal units
N = length(F);

% creating the main variables
Pi = optimvar('Pi', N);
Lambda = optimvar('L', 1);

C1 = zeros(N, N);
for i=1:length(C_Co)
    C1(i,i) = C_Co(i);
end

objective_fun = A_Co' * ones(N,1) + B_Co' * Pi + Pi' * C1 * Pi;
sum_F = sum(objective_fun);
problem = optimproblem("Objective", sum_F);

% power min and max constrains
problem.Constraints.min = Pi >= inputs(:,4);
problem.Constraints.max = Pi <= inputs(:,5);

% load generation equality constrain
problem.Constraints.pd = sum(Pi) == pload;

options = optimoptions('quadprog');
[sol,TotalCost,exitflag,output] = solve(problem, 'Options', options);

%% -------------------- CAACULATING LAMBDA --------------------
% selecting the unit that is not at its min or max power generation
% check the KKT rules for more informations
selected_p = 0;
index = 0;
for i=1:length(sol.Pi)
    index = index + 1;
    if sol.Pi(i) ~= inputs(i,4) && sol.Pi(i) ~= inputs(i,5)
       selected_p = sol.Pi(i);
       break;
    end
end

disp(selected_p)
disp(index)
lambda = subs(diff(F(index)), selected_p);

%% -------------------- RESULTS --------------------    
clc;
disp("Unit Power results : ")
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

fprintf("\n\n")

fprintf("Total cost    \t%.3f\n",TotalCost)

fprintf("\n\n")

fprintf("Lambda    \t\t%.3f\n",lambda)