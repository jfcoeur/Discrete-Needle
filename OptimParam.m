function [omega_init, kappa_c, s_star] = OptimParam(insertion_case, params_opt)

% Extract optimized parameters based on the insertion case

switch insertion_case
    
    case 'single_layer_C'
        
        omega_init = params_opt(1:3);
        kappa_c = params_opt(4);
        s_star = [];
    
    case 'double_layer_C'
        
        omega_init = params_opt(1:3);
        kappa_c = [params_opt(4), params_opt(5)];
        s_star = params_opt(6);
    
    case 'single_layer_S'
        
        omega_init = params_opt(1:3);
        kappa_c = params_opt(4);
        s_star = params_opt(5);
        
end

end % function OptimParam