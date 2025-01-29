function dkappa_i_ds = IntrinsicCurvatureDerivative(s, params, L, insertion_case)
% Compute the derivative of the intrinsic curvature kappa_i with respect to s.
%
% This function mirrors the piecewise definitions in IntrinsicCurvature.m
% and returns dkappa_i/ds as a 3x1 vector, assuming bending in the x-direction only.

switch insertion_case
    
    case 'single_layer_C'
        % Parameters: [kappa_c]
        kappa_c = params(1);
        dkappa_0_ds = -2*kappa_c*(1 - s/L)/L;
        dkappa_i_ds = [dkappa_0_ds; 0; 0];
    
    case 'double_layer_C'
        % Parameters: [kappa_c1, kappa_c2, s_star]
        kappa_c1 = params(1);
        kappa_c2 = params(2);
        s_star = params(3);
        
        if s <= s_star            
            dkappa_0_ds = kappa_c1 * (-2*(s_star - s)/L^2) + kappa_c2 * (1 - s_star/L) * (-2/L);                        
        else
            dkappa_0_ds = -2*kappa_c2*(1 - s/L)/L;
        end
        dkappa_i_ds = [dkappa_0_ds; 0; 0];
    
    case 'single_layer_S'
        % Parameters: [kappa_c, s_star]
        kappa_c = params(1);
        s_star = params(2);
        
        % Precompute powers:
        c1 = (s_star/L)^(2/3);
        c2 = (1 - s_star/L)^(2/3);
        
        if s < s_star
            C = kappa_c*c1;
            dkappa_0_ds = -2*C*(1 - s/L)/L;
        
        elseif s > s_star
            D = -kappa_c*c2;
            dkappa_0_ds = -2*D*(1 - s/L)/L;
            
        else
            % At s == s_star, the curvature definition is a single point expression.
            % The derivative can be defined as a limit approaching s_star.
            % For simplicity, choose the derivative from the left or right.
            % Here, we use the left-hand limit:
            C = kappa_c*c1;
            dkappa_0_ds = -2*C*(1 - s/L)/L;  % s == s_star
        end
        
        dkappa_i_ds = [dkappa_0_ds; 0; 0];
        
    otherwise
        error('Invalid insertion case in IntrinsicCurvatureDerivative function.');
end

end
