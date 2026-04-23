function F = Load_v(loads,mesh,element_type)

    if element_type == "linear_1D"
        nnode = 3;
    elseif element_type == "quadratic_1D"
        nnode = 3;
    elseif element_type == "triangular"
        nnode = 2;
    end

    all_global_nodes = [];
    for e = 1:length(mesh)
        all_global_nodes = [all_global_nodes, mesh(e).global_nodes];
    end

    num_nodes = max(all_global_nodes);

    F = zeros(nnode*num_nodes,1);

    for i = 1:size(loads,1)
        node = loads(i,1);
        Fx   = loads(i,2);
        Fy   = loads(i,3);
        if element_type == "triangular"
            F(nnode*node-1) = F(nnode*node-1) + Fx;
            F(nnode*node)   = F(nnode*node)   + Fy;
        else
            Mz = loads(i,4);   % nodal moment
            F(3*node-2) = F(3*node-2) + Fx;
            F(3*node-1) = F(3*node-1) + Fy;
            F(3*node)   = F(3*node)   + Mz;
        end
    end

end