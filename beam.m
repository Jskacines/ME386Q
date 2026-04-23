clear;clc;

nodal_geometry.type = "regular_grid";
nodal_geometry.nodes = [[0,0];[1,0];[0,1];[1,1];[2,0];[2,1];[3,0];[3,1];[4,0];[4,1];[5,0];[5,1];[6,0];[6,1]];
mesh = generate_mesh(nodal_geometry,"triangular",false);



plot_elements(mesh)

nodal_geometry.thickness = .1; % m
nodal_geometry.E = 10e6; % Pa
nodal_geometry.A = .1; % m 
nodal_geometry.I = .4; % m^4
nodal_geometry.v = .1; % poisson's ratio 

%% Set loads 
% [node_number magnitude_x magnitude_y]
loads_tri = [
    13   0   -1000 
    ];

Ftri = Load_v(loads_tri,mesh,"triangular"); % need to change between mesh/truss depending on what it is


%% Set boundary conditions
fixed_dofs_tri = [1 2 3 4 5 6];

% triangular 
ndof_tri = length(Ftri);
all_dofs_tri = 1:ndof_tri;
free_dofs_tri = setdiff(all_dofs_tri, fixed_dofs_tri);

%% Stiffness matrices
Ktri = K_matrix(nodal_geometry,"triangular",mesh);


%% Solve Uf * Kff = Pf & R + F = K*u for triangular elements
% Ku = F 

utri = zeros(ndof_tri,1);
utri(free_dofs_tri) = Ktri(free_dofs_tri, free_dofs_tri) \ Ftri(free_dofs_tri);
Rtri = Ktri*utri - Ftri;

%% plots
mesh = compute_deformations(mesh,utri);
% plot_displacement(mesh,1);
plot_displacement_and_force(mesh,10,Ftri,Rtri,5e-4);

