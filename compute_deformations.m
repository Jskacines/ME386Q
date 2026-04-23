function elem_struct = compute_deformations(elem_struct, u)
    if elem_struct(1).element_type == "beam"
        N_dof = 3;
    elseif elem_struct(1).element_type == "tri"
        N_dof = 2;
    end
    u_array = reshape(u, [N_dof, numel(u)/N_dof])';
    N_elems = numel(elem_struct);
    for i = 1:N_elems
        elem = elem_struct(i);
        node_deformations = zeros(size(elem.local_nodes,1),2);
        disp(size(elem.local_nodes))
        for j = 1:size(elem.local_nodes,1)
            global_idx = elem.global_nodes(j);
            node_deformations(j,:) = u_array(global_idx, 1:2);
        end
        elem_struct(i).node_deformations = node_deformations;
    end
end