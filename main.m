%% Clear workspace
clear; clc; close all;

%% Physical constants
% Needle diameter [mm]
d = 0.9; 

% Stiffness matrix
B = stiffmatrix(d);

% Insertion length [m]
L = 0.12;

% Insertion case
% Options: 'single_layer_C', 'double_layer_C', 'single_layer_S'
insertion_case = 'single_layer_C';


%% FBG Curvature Data
% Calibration using simulated data
[C, weights] = Csimuldata();

% Simulated data for demonstration
[s_FBG, kappa_FBG] = FBGsimuldata(C);


%% Optimization parameters
[params0, lb, ub] = InsertionCase(insertion_case);
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');

% Cost function handle
CostFunctionHandle = @(optimvars) CostFunction(optimvars, B, s_FBG, kappa_FBG, L, insertion_case);

%% Run optimization
[params_opt, fval] = fmincon(CostFunctionHandle, params0, [], [], [], [], lb, ub, [], options);

% Optimized parameters
[omega_init, kappa_c, s_star] = OptimParam(insertion_case, params_opt);

%% Compute needle shape and orientation along the shaft
% Solve Euler-Poincaré equations with optimized parameters
[kappa_model_s, R_s] = EulerPoincareSolver(params_opt, B, L, insertion_case);

% 3D coordinates
r = coordinates(R_s, kappa_model_s.s);