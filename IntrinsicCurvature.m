function [kappa_i] = IntrinsicCurvature(s, params, L, insertion_case)

% Intrinsic Curvature Function for Different Cases

switch insertion_case
    
    case 'single_layer_C'
        
        % Parameters: [kappa_c]
        kappa_c = params(1);
        kappa_0 = kappa_c * (1 - s / L)^2;
        kappa_i = [kappa_0; 0; 0];  % Assuming bending in x-direction
    
    case 'double_layer_C'
        
        % Parameters: [kappa_c1, kappa_c2, s_star]
        kappa_c1 = params(1);
        kappa_c2 = params(2);
        s_star = params(3);
        
        if s <= s_star
            kappa_0 = kappa_c1 * ((s_star - s) / L)^2 + kappa_c2 * (1 - s_star / L) * (1 + s_star / L - 2 * s / L);                      
        
        else
            kappa_0 = kappa_c2 * (1 - s / L)^2;
        
        end

        kappa_i = [kappa_0; 0; 0];  % Assuming bending in x-direction
    
    case 'single_layer_S'
        
        % Parameters: [kappa_c, s_star]
        kappa_c = params(1);
        s_star = params(2);
        
        if s < s_star
            kappa_0 = kappa_c * (s_star / L)^(2/3) * (1 - s / L)^2;
        
        elseif s == s_star
            kappa_0 = (kappa_c / 2) * ((s_star / L)^(2/3) - (1 - s_star / L)^(2/3)) * (1 - s_star / L)^2;
        
        else
            kappa_0 = -kappa_c * (1 - s_star / L)^(2/3) * (1 - s / L)^2;
        
        end

        kappa_i = [kappa_0; 0; 0];  % Assuming bending in x-direction
    
    otherwise
        error('Invalid insertion case in IntrinsicCurvature function.');

end

end 