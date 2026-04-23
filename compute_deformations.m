function elem_struct = compute_deformations(elem_struct, u)
    N_elems = numel(elem_struct);
    for i = 1:N_elems
        elem = elem_struct(i);
        node_deformations = zeros(size(elem.local_nodes,1),2);
        for j = 1:size(elem.local_nodes,1)
            global_idx = elem.global_nodes(j);
            dx = u(global_idx*2 - 1);
            dy = u(global_idx*2);
            node_deformations(j,:) = [dx,dy];
        end
        elem_struct(i).node_deformations = node_deformations;
    end
end