function truss = generate_truss(nodal_geometry, do_plots)
    if do_plots
        figure()
        scatter(nodal_geometry.nodes(:,1), nodal_geometry.nodes(:,2))
        axis padded
        axis equal
    end
    nodes = nodal_geometry.nodes;
    N_nodes = size(nodes, 1);
    conn = nodal_geometry.connections;
    N_elems = size(conn,1);
    sorted_nodes = sortrows(nodes, [1 2]);
    keys = num2cell(sorted_nodes,2);
    vals = (1:N_nodes)';
    global_node_dict = dictionary(keys, vals);
    elem_struct = struct('element_type', cell(N_elems,1), 'local_nodes', cell(N_elems,1), 'global_nodes', cell(N_elems,1));
    [elem_struct.element_type] = deal('beam');
    [elem_struct.local_to_global] = deal(global_node_dict);
    [elem_struct.global_to_local] = deal(dictionary(vals, keys));
    for k = 1:N_elems
        idx = conn(k,:);
        elem_struct(k).local_nodes = [nodes(idx(1),:);nodes(idx(2),:)];
        elem_struct(k).global_nodes = global_node_dict({nodes(idx(1),:),nodes(idx(2),:)});
    end
    truss = elem_struct;
end