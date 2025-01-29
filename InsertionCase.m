function [params0, lb, ub] = InsertionCase(insertion_case)

% We set omega_init_z to 0, since we assume they are.

switch insertion_case
    
    case 'single_layer_C'
        
        % Optimization variables: [omega_init_x, omega_init_y, omega_init_z, kappa_c]
        params0 = [0, 0, 0, 0.1];        % Initial guesses
        lb = [-Inf, -Inf, 0, 0];      % Lower bounds
        ub = [Inf, Inf, 0, Inf];       % Upper bounds
    
    case 'double_layer_C'
        
        % Optimization variables: [omega_init_x, omega_init_y, omega_init_z, kappa_c1, kappa_c2, s_star]
        params0 = [0, 0, 0, 0.1, 0.05, L/2];  % Initial guesses
        lb = [-Inf, -Inf, 0, 0, 0, 0];     % Lower bounds
        ub = [Inf, Inf, 0, Inf, Inf, L];    % Upper bounds
    
    case 'single_layer_S'
        
        % Optimization variables: [omega_init_x, omega_init_y, omega_init_z, kappa_c, s_star]
        params0 = [0, 0, 0, 0.1, L/2];   % Initial guesses
        lb = [-Inf, -Inf, 0, 0, 0];   % Lower bounds
        ub = [Inf, Inf, 0, Inf, L];    % Upper bounds
    
    otherwise
        error('Invalid insertion case selected.');

end

end  % function InsertionCase
