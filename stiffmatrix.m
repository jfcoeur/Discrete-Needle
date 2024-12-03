function [B] = stiffmatrix(d)

% d = 0.9 (diameter in mm)

% Stpainless Steel 304
E = 200e9*1e-6; % 200 GPa, conversion from N/m^2 to N/mm^2
Pratio = 0.29; % Poisson's ratio

% Shear modulus
G = E/2/(1 + Pratio);

% Bending moment of inertia for a circular cross-section
I = pi*d^4/64;

% Polar moment of inertia for a circular cross-section
J = pi*d^4/32;

% Stiffness
BendStiff = E*I; % Bending
TorStiff = G*J; % Torsional

% Stiffness matrix
B = diag([BendStiff,BendStiff,TorStiff]);

end % function stiffmatrix