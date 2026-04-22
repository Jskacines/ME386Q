clear;clc;

nodal_geometry.type = "regular_grid";
nodal_geometry.nodes = [[0,0];[1,0];[0,1];[1,1];[2,0];[2,1]];
mesh = generate_mesh(nodal_geometry,"triangular",false);

nodal_geometry.connections = [[1,2];[1,3];[2,3];[2,4];[3,4];[2,5];[4,5];[4,6];[5,6]];

N_subdivisions = 2;
truss = generate_truss(nodal_geometry, N_subdivisions, false);

plot_elements(mesh)
plot_elements(truss)

nodal_geometry.thickness = .1; % m
nodal_geometry.E = 10e6; % Pa
nodal_geometry.A = .1; % m 
nodal_geometry.v = .1; % poisson's ratio 

%% Set loads 
% [node_number magnitude_x magnitude_y]
loads = [
    6   0   -1000
    5 200      0
    ];
F = Load_v(loads,mesh); % need to change between mesh/truss depending on what it is
%% Set boundary conditions
fixed_dofs = [1 2 3 4]; % 1 = 1x, 2 = 1y, 4 = 2y
ndof = length(F);
all_dofs = 1:ndof;
free_dofs = setdiff(all_dofs, fixed_dofs);


%% Stiffness matrices
Klin1D = K_matrix(nodal_geometry, "linear_1D",truss);
Ktri = K_matrix(nodal_geometry,"triangular",mesh);
Kquad1D = K_matrix(nodal_geometry, "quadratic_1D", truss);


%% Solve Uf * Kff = Pf & R + F = K*u
% Ku = F 

u = zeros(ndof,1);
u(free_dofs) = Ktri(free_dofs, free_dofs) \ F(free_dofs)
R = Ktri*u - F
