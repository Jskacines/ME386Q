nodal_geometry.type = "regular_grid";
nodal_geometry.nodes = [[0,0];[1,0];[0,1];[1,1];[2,0];[2,1]];
mesh = generate_mesh(nodal_geometry,"triangular",false);

nodal_geometry.connections = [[1,2];[1,3];[2,3];[2,4];[3,4];[2,5];[4,5];[4,6];[5,6]];

N_subdivisions = 2;
truss = generate_truss(nodal_geometry, N_subdivisions, false);

plot_elements(mesh)
plot_elements(truss)