function [C, weights] = Cmatrix(wvl0, wvl90, k, alpha)

% Inputs: wvl0 & wvl90 (k x fiber x AA), k (N x 1), alpha (1 x 2) (for 0 and 90) 

% Curvature weights
bool = k <= 1;
w = double(bool);
w(w == 0) = 0.05;
W = diag([w; w]);

% Curvatures
kx = reshape(k*cosd(alpha), [], 1);
ky = reshape(k*sind(alpha), [], 1);
k = [kx, ky];

% Calibration matrix
for i = 1:size(wvl0,3)
    AA = [wvl0(:,:,i) ; wvl90(:,:,i)];
    C(:,:,i) = (AA' * W * AA) \ (AA' * W * k);

    % Mean Squared Error (MSE)
    k_est(:,:,i) = AA*C(:,:,i);
    err = k - k_est(:,:,i);
    MSE(i) = mean(err.^2, "all");
end

% Reliability weights for each sensing area
weights = 1./MSE;
weights = weights/sum(weights); % Normalize so sum = 1

% Weighted shape
for j = 1:size(k_est,3)
    k_weight(:,:,j) = weights(j)*k_est(:,:,j);
end
k_weight = sum(k_weight,3);

end % function Cmatrix
