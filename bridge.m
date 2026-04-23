clear;clc;

nodal_geometry.nodes = [
    % Bottom chord (deck level) — 9 nodes spanning x = 0..8
    [0,0]; [1,0]; [2,0]; [3,0]; [4,0]; [5,0]; [6,0]; [7,0]; [8,0];
    % Top chord (arch profile) — 7 nodes, low parabolic rise, symmetric
    [1,0.5]; [2,1]; [3,1.5]; [4,1.75]; [5,1.5]; [6,1]; [7,0.5]
];

nodal_geometry.connections = [
    % Bottom chord
    [1,2]; [2,3]; [3,4]; [4,5]; [5,6]; [6,7]; [7,8]; [8,9];
    % Top chord (arch)
    [10,11]; [11,12]; [12,13]; [13,14]; [14,15]; [15,16];
    % Verticals
    [2,10]; [3,11]; [4,12]; [5,13]; [6,14]; [7,15]; [8,16];
    % Diagonals (mirrored about centre)
    [3,10]; [4,11]; [5,12]; [5,14]; [6,15]; [7,16];
    % End portals
    [1,10]; [9,16]
];

disp(size(nodal_geometry.connections))


N_subdivisions = 0;
truss = generate_truss(nodal_geometry, N_subdivisions, false);

plot_elements(truss)

nodal_geometry.thickness = .1; % m
nodal_geometry.E = 10e6; % Pa
nodal_geometry.A = .1; % m 
nodal_geometry.I = .4; % m^4
nodal_geometry.v = .1; % poisson's ratio 

%% Set loads 
% [node_number magnitude_x magnitude_y]

loads_1D = [
    8   0   -1000 
    ];
Flin1D = Load_v(loads_1D,truss, "linear_1D");
Fquad1D = Load_v(loads_1D,truss,"quadratic_1D")


%% Set boundary conditions
fixed_dofs_1D = [1 2 31 32]; % 1 = 1x, 2 = 1y,3 = 1rot, 4 = 2x,5 = 2y, 6 = 2rot
fixed_dofs_quad1D = [1 2 31 32];

% linear
ndof_1D = length(Flin1D);
all_dofs_1D = 1:ndof_1D;
free_dofs_1D = setdiff(all_dofs_1D, fixed_dofs_1D);

% quadratic 
% linear
ndof_quad1D = length(Fquad1D);
all_dofs_quad1D = 1:ndof_quad1D;
free_dofs_quad1D = setdiff(all_dofs_quad1D, fixed_dofs_quad1D);

%% Stiffness matrices
Klin1D = K_matrix(nodal_geometry, "linear_1D",truss);
Kquad1D = K_matrix(nodal_geometry, "quadratic_1D", truss)



%% Solve Uf * Kff = Pf & R + F = K*u for linear 1D elements

ulin1D = zeros(ndof_1D,1);
ulin1D(free_dofs_1D) = Klin1D(free_dofs_1D, free_dofs_1D) \ Flin1D(free_dofs_1D);
Rlin1D = Klin1D*ulin1D - Flin1D;

%% Solve Uf * Kff = Pf & R + F = K*u for qudratic 1D elements

uquad1D = zeros(ndof_quad1D,1);
uquad1D(free_dofs_quad1D) = Kquad1D(free_dofs_quad1D, free_dofs_quad1D) \ Fquad1D(free_dofs_quad1D)
Rquad1D = Kquad1D*uquad1D - Fquad1D

%% plots

truss = compute_deformations(truss,ulin1D);
plot_displacement(truss,20);
plot_displacement_and_force(truss,20,Flin1D,Rlin1D,1e-3);

