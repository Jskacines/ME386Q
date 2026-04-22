nodal_geometry.type = "regular_grid";
nodal_geometry.nodes = [[0,0];[1,0];[0,1];[1,1];[2,0];[2,1]];
mesh = generate_mesh(nodal_geometry,"triangular",false);

nodal_geometry.connections = [[1,2];[1,3];[2,3];[2,4];[3,4];[2,5];[4,5];[4,6];[5,6]];
truss = generate_truss(nodal_geometry, false);

plot_elements(mesh)
plot_elements(truss)

nodal_geometry.thickness = .1; % m
nodal_geometry.E = 10e6; % Pa
nodal_geometry.A = .1; % m 
nodal_geometry.v = .1; % poisson's ratio 


element_nodes = mesh.local_nodes;
global_nodes = mesh.global_nodes;


N = size(nodal_geometry.nodes);

K = K_matrix(nodal_geometry, "triangular", mesh);
