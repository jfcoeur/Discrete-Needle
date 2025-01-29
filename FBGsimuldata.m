function [s_FBG, kappa_FBG] = FBGsimuldata(C,N,L)

% Experimental shape (values per sensing area (N))
kappa = RandomWalk(4*rand, 2*rand, N);
alpha = zeros(4,1); % x curvature only

% FBG positions
if N == 4
    s_FBG = [0.0190; 0.0540; 0.0890; 0.1090]; % Arc length positions of FBG sensors
else
    s_FBG = linspace(0, L, N);
end

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

%%%%%%%%
function out = RandomWalk(initial_value, step_size, N)

    out = zeros(N,1);
    out(1) = initial_value;
    for j = 2:N
        out(j) = out(j-1) + step_size*randn;
    end

end
%%%%%%%%

end % function FBGsimuldata
