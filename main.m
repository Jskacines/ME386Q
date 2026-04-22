nodal_geometry.type = "regular_grid";
nodal_geometry.nodes = [[0,0];[1,1];[1,0];[0,1];[2,0];[2,1]];
mesh = generate_mesh(nodal_geometry,"triangular",true);