clear;clc;

nodal_geometry.type = "regular_grid";
nodal_geometry.nodes = [[0,0];[1,0];[0,1];[1,1];[2,0];[2,1]];
mesh = generate_mesh(nodal_geometry,"triangular",false);

nodal_geometry.connections = [[1,2];[1,3];[2,3];[2,4];[3,4];[2,5];[4,5];[4,6];[5,6]];

N_subdivisions = 0;
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
loads_tri = [
    6   0   -1000 
    5 200      0 
    ];

loads_1D = [
    6   0   -1000 
    ];
Ftri = Load_v(loads_tri,mesh,"triangular"); % need to change between mesh/truss depending on what it is
Flin1D = Load_v(loads_1D,truss, "linear_1D");
Fquad1D = Load_v(loads_1D,truss,"quadratic_1D")


%% Set boundary conditions
fixed_dofs_1D = [1 2 3 4]; % 1 = 1x, 2 = 1y,3 = 1rot, 4 = 2x,5 = 2y, 6 = 2rot
fixed_dofs_quad1D = [1 2 3 4];
fixed_dofs_tri = [1 2 3 4 5 6];

% triangular 
ndof_tri = length(Ftri);
all_dofs_tri = 1:ndof_tri;
free_dofs_tri = setdiff(all_dofs_tri, fixed_dofs_tri);


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
Ktri = K_matrix(nodal_geometry,"triangular",mesh);
Kquad1D = K_matrix(nodal_geometry, "quadratic_1D", truss)


%% Solve Uf * Kff = Pf & R + F = K*u for triangular elements
% Ku = F 

utri = zeros(ndof_tri,1);
utri(free_dofs_tri) = Ktri(free_dofs_tri, free_dofs_tri) \ Ftri(free_dofs_tri);
Rtri = Ktri*utri - Ftri;

%% Solve Uf * Kff = Pf & R + F = K*u for linear 1D elements

ulin1D = zeros(ndof_1D,1);
ulin1D(free_dofs_1D) = Klin1D(free_dofs_1D, free_dofs_1D) \ Flin1D(free_dofs_1D);
Rlin1D = Klin1D*ulin1D - Flin1D;

%% Solve Uf * Kff = Pf & R + F = K*u for qudratic 1D elements

uquad1D = zeros(ndof_quad1D,1);
uquad1D(free_dofs_quad1D) = Kquad1D(free_dofs_quad1D, free_dofs_quad1D) \ Fquad1D(free_dofs_quad1D)
Rquad1D = Kquad1D*uquad1D - Fquad1D

%% plots
mesh = compute_deformations(mesh,utri);
% plot_displacement(mesh,1);
plot_displacement_and_force(mesh,10,Ftri,Rtri,1e-3);

truss = compute_deformations(truss,ulin1D);
plot_displacement(truss,20);
plot_displacement_and_force(truss,20,Flin1D,Rlin1D,1e-3);

