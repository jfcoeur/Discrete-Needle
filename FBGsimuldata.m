function [s_FBG, kappa_FBG] = FBGsimuldata(C)

% Experimental shape (values per sensing area (4))
kappa = 4*rand(4,1);
alpha = 360*rand(4,1);

% x-y Curvatures
kx = kappa.*cosd(alpha);
ky = kappa.*sind(alpha);
k = [kx, ky];

% FBG positions
s_FBG = [0.0190; 0.0540; 0.0890; 0.1090]; % Arc length positions of FBG sensors

% Wavelength shift
wvl_shift = simuldata(kappa, alpha);

% Temperature compensation
wvl_shift = Tcomp(wvl_shift);

% FBG experimental curvature
for i = 1:size(wvl_shift,2)
    
    AA = wvl_shift(:,i);
    kappa_FBG(i,:) = AA'*C(:,:,i);

end

kappa_FBG(:,3) = zeros(size(kappa_FBG(:,1)));

end % function FBGsimuldata