function [r] = coordinates(R_s, s_vals)

% Compute position along the needle by integrating the tangent
N = length(s_vals);
r = zeros(N, 3);
e1 = [1; 0; 0];  % Initial tangent vector

for i = 2:N
    
    ds = s_vals(i) - s_vals(i-1);
    
    % Tangent vector at position s is R(s) * e1
    t = R_s(:,:,i-1) * e1;
    
    % Update position by integrating the tangent vector
    r(i,:) = r(i-1,:) + t' * ds;

end

end % function coordinates