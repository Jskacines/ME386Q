function K_matrix = K_matrix(nodal_geometry, element_type,mesh)

    % N is number of nodes 
    
    %% Building global stiffness matrix 

    if element_type == "triangular"
        K_matrix = Ktri_global(mesh,nodal_geometry);

    elseif element_type == "linear_1D"
        K_matrix = Klin1D_global(mesh,nodal_geometry);
    elseif element_type == "quadratic_1D"
        K_matrix = Kquad1D_global(mesh,nodal_geometry);
    % elseif element_type == "quadratic"
    %     K_matrix = Kquad(N,nodal_geometry.E,nodal_geometry.A);
    else 
        fprintf("Unsupported element type")
    end

end

function K_linear = Klin1D(E,A,element_nodes)

    % element_nodes = [x1 y1; x2 y2]

    x1 = element_nodes(1,1);
    y1 = element_nodes(1,2);
    x2 = element_nodes(2,1);
    y2 = element_nodes(2,2);

    L = sqrt((x2-x1)^2 + (y2-y1)^2);

    if L == 0
        error('Zero-length truss element detected.');
    end

    c = (x2-x1)/L;
    s = (y2-y1)/L;

    K_linear = (E*A/L) * ...
        [ c^2   c*s   -c^2   -c*s
          c*s   s^2   -c*s   -s^2
         -c^2  -c*s    c^2    c*s
         -c*s  -s^2    c*s    s^2 ];

end

function K_global = Klin1D_global(mesh,nodal_geometry)
    % nnode = 2; % line has 2 nodes 
    E = nodal_geometry.E;
    A = nodal_geometry.A;

    all_global_nodes = [];
    for e = 1:length(mesh)
        all_global_nodes = [all_global_nodes, mesh(e).global_nodes];
    end
    num_nodes = max(all_global_nodes);

    K_global = zeros(2*num_nodes, 2*num_nodes);

    for e = 1:length(mesh)

        % element geometry
        element_nodes = mesh(e).local_nodes;

        % element stiffness
        K_1D = Klin1D(E,A,element_nodes);

        % node connectivity
        global_conn = mesh(e).global_nodes;
        n1 = global_conn(1);
        n2 = global_conn(2);

        % DOF connectivity: [u1 v1 u2 v2]
        edof = [
            2*n1-1, 2*n1, ...
            2*n2-1, 2*n2
        ];

        % global assembly
        K_global(edof,edof) = K_global(edof,edof) + K_1D;
    end
end

function K_quadratic_1D = Kquad1D(E,A,element_nodes)

    x1 = element_nodes(1,1);
    y1 = element_nodes(1,2);
    x2 = element_nodes(3,1);
    y2 = element_nodes(3,2);

    L = sqrt((x2-x1)^2 + (y2-y1)^2);

    if L == 0
        error('Zero-length element detected.');
    end

    c = (x2-x1)/L;
    s = (y2-y1)/L;

    K_quadratic_1D = (E*A/L) * ...
        [ c^2   c*s   -c^2   -c*s
          c*s   s^2   -c*s   -s^2
         -c^2  -c*s    c^2    c*s
         -c*s  -s^2    c*s    s^2 ];
end

function K_global = Kquad1D_global(mesh,nodal_geometry)

    E = nodal_geometry.E;
    A = nodal_geometry.A;

    all_global_nodes = [];
    for e = 1:length(mesh)
        all_global_nodes = [all_global_nodes, mesh(e).global_nodes];
    end
    num_nodes = max(all_global_nodes);

    K_global = zeros(2*num_nodes, 2*num_nodes);

    for e = 1:length(mesh)

        % Endpoints from the linear mesh
        element_nodes = mesh(e).local_nodes;

        x1 = element_nodes(1,1);
        y1 = element_nodes(1,2);
        x2 = element_nodes(2,1);
        y2 = element_nodes(2,2);

        % Midpoint only for geometric interpolation bookkeeping
        xm = (x1 + x2)/2;
        ym = (y1 + y2)/2;

        element_nodes_quad = [
            x1 y1;
            xm ym;
            x2 y2
        ];

        % Condensed quadratic truss stiffness
        K_e = Kquad1D(E,A,element_nodes_quad);

        global_conn = mesh(e).global_nodes;
        n1 = global_conn(1);
        n2 = global_conn(2);

        edof = [
            2*n1-1, 2*n1, ...
            2*n2-1, 2*n2
        ];

        K_global(edof,edof) = K_global(edof,edof) + K_e;
    end
end

function K_tri = Ktri(E,t,v, element_nodes) % element thickness is t 

    D_mat = D_matrix(E,v);
    A = 1/2 * (element_nodes(1,1) * (element_nodes(2,2)-element_nodes(3,2)) + element_nodes(2,1) * (element_nodes(3,2)-element_nodes(1,2)) + element_nodes(3,1) * (element_nodes(1,2) - element_nodes(2,2)));
    B_mat = B_matrix(A,element_nodes);
    % local stiffness matrix 
    K_tri = t * A * transpose(B_mat) * D_mat *B_mat;

end

function K_global = Ktri_global(mesh,nodal_geometry)
    
    nnode = 3; % triangle has 3 nodes 
    E = nodal_geometry.E;
    t = nodal_geometry.thickness;
    v = nodal_geometry.v;

    num_nodes = size(nodal_geometry.nodes,1);
    K_global = zeros(2*num_nodes, 2*num_nodes);
    
    for e = 1:length(mesh)

        % element geometry
        element_nodes = mesh(e).local_nodes;
    
        % element stiffness
        K_tri = Ktri(E,t,v,element_nodes);
    
        % node connectivity
        nnode = size(element_nodes,1);
        global_conn = zeros(1,nnode);
    
        for i = 1:nnode
            coord = element_nodes(i,:);
            global_conn(i) = mesh(e).local_to_global({coord});
        end
    
        % DOF connectivity
        edof = zeros(1,2*nnode);
        for i = 1:nnode
            n = global_conn(i);
            edof(2*i-1:2*i) = [2*n-1, 2*n];
        end
    
        % global assembly
        K_global(edof,edof) = K_global(edof,edof) + K_tri;
    end
end

function B_matrix = B_matrix(A,element_nodes)
    beta1 = element_nodes(2,2) - element_nodes(3,2); % where x coords are column 1, y coords are column 2 
    beta2 = element_nodes(3,2) - element_nodes(1,2);
    beta3 = element_nodes(1,2) - element_nodes(2,2);
    gamma1 = element_nodes(3,1) - element_nodes(2,1); 
    gamma2 = element_nodes(1,1) - element_nodes(3,1);
    gamma3 = element_nodes(2,1) - element_nodes(1,1);
    B_matrix = 1/(2*A) * [beta1 0 beta2 0 beta3 0 ; 0 gamma1 0 gamma2 0 gamma3 ; gamma1 beta1 gamma2 beta2 gamma3 beta3];
end

function D_matrix = D_matrix(E,v)
    D_matrix = E/(1-v^2) * [ 1 v 0 ; v 1 0 ; 0 0 (1-v)/2];
end
