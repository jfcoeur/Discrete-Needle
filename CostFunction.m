function [cost] = CostFunction(optimvars, B, s_FBG, kappa_FBG, L, insertion_case, weights)

% Cost Function for Optimization

% Solve Euler-Poincaré equations to get model curvature
[kappa_model_s, ~] = EulerPoincareSolver(optimvars, B, L, insertion_case);

% Interpolate model curvature at FBG sensor positions
kappa_model_FBG = interp1(kappa_model_s.s, kappa_model_s.kappa, s_FBG, 'linear');

% Compute cost as the sum of squared differences
diff = kappa_FBG - kappa_model_FBG;
diff_wt = weights'.*diff;
cost = sum(vecnorm(diff_wt, 2, 2).^2);

end % function CostFunction
