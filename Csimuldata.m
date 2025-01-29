function [C, weights] = Csimuldata(N)

% Calibration shapes
kappa = [0; 0.5; 1.6; 2; 2.5; 3.2; 4]; % [m^-1]
alpha0 = 270*rand; % [deg]
alpha = [alpha0, alpha0 + 90]; % [deg]

% Wavelength shifts
for i = 1:length(kappa)

    wvl0 = simuldata(kappa(i)*ones(N,1), alpha(1)*ones(N,1));
    wvl90 = simuldata(kappa(i)*ones(N,1), alpha(2)*ones(N,1));

    % Temperature compensation
    wvl0_Tcomp(i,:,:) = Tcomp(wvl0);
    wvl90_Tcomp(i,:,:) = Tcomp(wvl90);

end

% Calibration matrix
[C, weights] = Cmatrix(wvl0_Tcomp, wvl90_Tcomp, kappa, alpha);

end % function Csimuldata
