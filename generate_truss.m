function truss = generate_truss(nodal_geometry, N_subdivisions, do_plots)
    if do_plots
        figure()
        scatter(nodal_geometry.nodes(:,1), nodal_geometry.nodes(:,2))
        axis padded
        axis equal
    end
    nodes = nodal_geometry.nodes;
    conn = nodal_geometry.connections;

    [nodes, conn] = subdivide_nodes(nodes, conn, N_subdivisions);

    N_nodes = size(nodes, 1);
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

function [nodes, conn] = subdivide_nodes(nodes0, conn0, N_subdivisions)
    if N_subdivisions == 0
        nodes = nodes0;
        conn = conn0;
        return
    elseif N_subdivisions == 1
        nodes = nodes0;
        conn = [];
        for i = 1:size(conn0,1)
            idx = conn0(i,:);
            dx = (nodes0(idx(2),:) - nodes0(idx(1),:))/2;
            new_node = nodes0(idx(1),:) + dx;
            nodes = [nodes; new_node];
            k = size(nodes,1);
            conn = [conn; [idx(1),k]];
            conn = [conn; [k,idx(2)]];
        end
        return
    else
        nodes = nodes0;
        conn = [];
        for i = 1:size(conn0,1)
            idx = conn0(i,:);
            dx = (nodes0(idx(2),:) - nodes0(idx(1),:))/(N_subdivisions+1);
            % disp(dx)
            for j = 1:N_subdivisions
                new_node = nodes0(idx(1),:) + dx*j;
                nodes = [nodes; new_node];
                k = size(nodes,1);
                if j == 1
                    conn = [conn; [idx(1),k]];
                elseif j == N_subdivisions
                    conn = [conn; [k-1,k]];
                    conn = [conn; [k,idx(2)]];
                else
                    conn = [conn;[k-1,k]]
                end
            end
        end
    end
end