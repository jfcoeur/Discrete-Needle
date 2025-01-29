function [wvl_shift] = simuldata(k, alpha)

% Cross-sectional dimensions
tet12 = 120; % [deg]
tet23 = 120; % [deg]
d = 70E-6; % [m]
r = [d*cosd(alpha), d*cosd(tet12 - alpha), d*cosd(tet12 + tet23 - alpha)]; % [m]

% Strain
sig = k.*r;

% Wavelength shift
k_e = 6.668e-6; % strain
wvl_shift = k_e*sig';

end % function simuldata
