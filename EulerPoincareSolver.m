function [kappa_model_s, R_s] = EulerPoincareSolver(optimvars, B, L, insertion_case)

% Euler-Poincaré Equations Solver

% Initilization
N = 100; % Number of points
s_vals = linspace(0, L, N)'; % Longitudinal positions
kappa = zeros(N,3); % Curvature values
R_s = zeros(3,3,N); % Rotation matrices at position s
R = eye(3); % Initial rotation matrix assuming initial rotation is identity
R_s(:,:,1) = R;

% Extract parameters
kappa0 = optimvars(1:3); % Inital curvature
kappa_params = optimvars(4:end); % Intrinsic curvatures

% Set up ODE solver options
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

% Solve the Euler-Poincaré equations using an ODE solver
[s_out, kappa_out] = ode45(@(s, k) EulerPoincareODE(s, k, kappa_params, B, L, insertion_case), s_vals, kappa0, options);

% Reshape kappa_out to get curvature at each s
kappa = kappa_out; % reshape(kappa_out', [3, N])'

% Compute rotation matrices along the needle
for i = 2:N

    ds = s_vals(i) - s_vals(i-1);
    omega = kappa(i-1,:)';

    % Update rotation matrix using exponential map
    R = R * expm(skew(omega * ds));
    R_s(:,:,i) = R;

end

% Output curvature and rotation matrices
kappa_model_s.s = s_vals;
kappa_model_s.kappa = kappa;

end % function EulerPoincareSolver