function mesh = generate_mesh(nodal_geometry, element_type, do_plots)
    if do_plots
        figure()
        scatter(nodal_geometry.nodes(:,1), nodal_geometry.nodes(:,2))
        axis padded
        axis equal
    end
    if nodal_geometry.type == "regular_grid"
        if element_type == "triangular"
            mesh = regular_tri_mesh(nodal_geometry.nodes);
        elseif element_type == "quadrilateral"
            mesh = regular_quad_mesh(nodal_geometry.nodes);
        else
            fprintf("Unsupported element type")
        end
    else
        fprintf("Only regular grid nodal geometry is supported")
    end
end



function mesh = regular_tri_mesh(nodes)
    N_nodes = size(nodes, 1);
    dim = size(nodes, 2);
    sorted_nodes = sortrows(nodes, [1 2]);
    keys = num2cell(sorted_nodes,2);
    vals = (1:N_nodes)';
    global_node_dict = dictionary(keys, vals);
    N_x = numel(unique(sorted_nodes(:,1)));
    N_y = numel(unique(sorted_nodes(:,2)));
    node_grid = pagetranspose(reshape(sorted_nodes, [N_y, N_x, dim]));
    N_elems = (N_x-1)*(N_y-1)*2;
    elem_struct = struct('element_type', cell(N_elems,1), 'local_nodes', cell(N_elems,1), 'global_nodes', cell(N_elems,1));
    [elem_struct.element_type] = deal('tri');
    [elem_struct.local_to_global] = deal(global_node_dict);
    [elem_struct.global_to_local] = deal(dictionary(vals, keys));
    k = 1;
    for i = 1:(N_x-1)
        for j = 1:(N_y-1)
            nodes = [
                    squeeze(node_grid(i,j,:))';
                    squeeze(node_grid(i+1,j,:))';
                    squeeze(node_grid(i,j+1,:))'
                    ];
            elem_struct(k).local_nodes = [nodes(1,:);nodes(2,:);nodes(3,:)];
            elem_struct(k).global_nodes = global_node_dict({nodes(1,:),nodes(2,:),nodes(3,:)});

            nodes = [
                    squeeze(node_grid(i+1,j+1,:))';
                    squeeze(node_grid(i,j+1,:))';
                    squeeze(node_grid(i+1,j,:))'
                    ];
            elem_struct(k+1).local_nodes = [nodes(1,:);nodes(2,:);nodes(3,:)];
            elem_struct(k+1).global_nodes = global_node_dict({nodes(1,:),nodes(2,:),nodes(3,:)});
            k = k + 2;
        end
    end

    mesh = elem_struct;
end

function mesh = regular_quad_mesh(node_grid)
    mesh = NaN;
end

