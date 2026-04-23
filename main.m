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
nodal_geometry.I = .4; % m^4
nodal_geometry.v = .1; % poisson's ratio 

%% Set loads 
% [node_number magnitude_x magnitude_y]
loads = [
    6   0   -1000
    5 200      0
    ];
Ftri = Load_v(loads,mesh,"triangular"); % need to change between mesh/truss depending on what it is
F1D = Load_v(loads,truss, "linear_1D");

%% Set boundary conditions
fixed_dofs = [1 2 3 4 5 6]; % 1 = 1x, 2 = 1y,3 = 1rot, 4 = 2x,5 = 2y, 6 = 2rot

% triangular 
ndof_tri = length(Ftri);
all_dofs = 1:ndof_tri;
free_dofs_tri = setdiff(all_dofs, fixed_dofs);

% linear
ndof_lin = length(F1D);
all_dofs = 1:ndof_lin;
free_dofs_lin = setdiff(all_dofs, fixed_dofs);

%% Stiffness matrices
Klin1D = K_matrix(nodal_geometry, "linear_1D",truss);
Ktri = K_matrix(nodal_geometry,"triangular",mesh);
% Kquad1D = K_matrix(nodal_geometry, "quadratic_1D", truss);


%% Solve Uf * Kff = Pf & R + F = K*u for triangular elements
% Ku = F 

utri = zeros(ndof_tri,1);
utri(free_dofs_tri) = Ktri(free_dofs_tri, free_dofs_tri) \ Ftri(free_dofs_tri)
Rtri = Ktri*utri - Ftri

%% Solve Uf * Kff = Pf & R + F = K*u for triangular elements
size(Klin1D)
size(F1D)
ulin1D = zeros(ndof_lin,1);
ulin1D(free_dofs_lin) = Klin1D(free_dofs_lin, free_dofs_lin) \ F1D(free_dofs_lin)
size(ulin1D)
Rlin1D = Klin1D*ulin1D - F1D