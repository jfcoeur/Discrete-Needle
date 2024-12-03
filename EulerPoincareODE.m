function [dkappa_ds] = EulerPoincareODE(s, kappa, kappa_params, B, L, insertion_case)

% Euler-Poincaré ODE Function

% Intrinsic curvature at position s
kappa_i = IntrinsicCurvature(s, kappa_params, L, insertion_case);

% Compute internal moment M
M = B * (kappa - kappa_i);

% Compute derivative of M
dM_ds = -cross(kappa, M);

% Compute derivative of kappa
dkappa_ds = B \ dM_ds;

end % function EulerPoincareODE