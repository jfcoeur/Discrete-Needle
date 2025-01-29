function [dkappa_ds] = EulerPoincareODE(s, kappa, kappa_params, B, L, insertion_case)

    % Intrinsic curvature and its derivative
    kappa_i = IntrinsicCurvature(s, kappa_params, L, insertion_case);
    dkappa_i_ds = IntrinsicCurvatureDerivative(s, kappa_params, L, insertion_case);

    % Compute internal moment M
    M = B*(kappa - kappa_i);

    % From Euler-Poincaré:
    % B(dkappa/ds - dkappa_i/ds) + cross(kappa, M) = 0
    % => B(dkappa/ds) = cross(kappa, M)*(-1) + B(dkappa_i/ds)
    % => dkappa/ds = dkappa_i/ds - B^{-1}(cross(kappa, M))

    dkappa_ds = dkappa_i_ds - (B \ cross(kappa, M));

end
