function K_matrix = K_matrix(nodal_geometry, element_type,mesh)

    % N is number of nodes 
    
    %% Building global stiffness matrix 

    if element_type == "triangular"
        K_matrix = Ktri_global(mesh,nodal_geometry);

    elseif element_type == "linear_1D"
        K_matrix = Klin1D(N,nodal_geometry.E,nodal_geometry.A);
    elseif element_type == "linear2D"
        K_matrix = Klin2D(N,nodal_geometry.E,nodal_geometry.A);
    elseif element_type == "quadratic"
        K_matrix = Kquad(N,nodal_geometry.E,nodal_geometry.A);
    else 
        fprintf("Unsupported element type")
    end

end

function K_linear = Klin1D(N,E,A)

    K_linear(1,1) = 1; 
    K_linear(1,2) = -1; 
    K_linear(N,N) = 1; 
    K_linear(N,N-1) = -1; 
    
    c = 0;
    
    for i = 2:N-1
        for j = 2:N-1
    
            K_linear(i,j-1+c) = E*A/N * -1;
            K_linear(i,j+c) = E*A/N * 2;
            K_linear(i,j+1+c) = E*A/N * -1;
    
            c = c+1;
            break
        end
    end

end

function K_linear = Klin2D(N,E,A)

    K_linear = NaN;

end

% function K_quad = Kquad(N,E,A)
% 
%     K_quad = NaN;
% 
% end

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
