function plot_displacement(elem_struct,scale_factor)
    N_elems = numel(elem_struct);
    figure()
    hold on
    color1 = [0,0,1];
    color2 = [0,1,0];
    % colors = parula(N_elems);
    % colormap(parula(N_elems))
    clim([0.5 N_elems + 0.5])
    % Initial Mesh
    for i = 1:N_elems
        elem = elem_struct(i);
        if elem.element_type == "beam"
            plot(elem.local_nodes(:,1), elem.local_nodes(:,2),"Color",color1)
        elseif elem.element_type == "tri"
            patch(elem.local_nodes(:,1),elem.local_nodes(:,2),color1)
        end
        s = scatter(elem.local_nodes(:,1), elem.local_nodes(:,2),...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
        s.SizeData = 20;
        % s.AlphaData = 0.5;
        alpha(0.3)
    end
    % Deformed Mesh
    for i = 1:N_elems
        elem = elem_struct(i);
        if elem.element_type == "beam"
            plot(elem.local_nodes(:,1) + elem.node_deformations(:,1)*scale_factor,... 
            elem.local_nodes(:,2) + elem.node_deformations(:,2)*scale_factor,"Color",color2)
        elseif elem.element_type == "tri"
            patch(elem.local_nodes(:,1) + elem.node_deformations(:,1)*scale_factor,...
            elem.local_nodes(:,2) + elem.node_deformations(:,2)*scale_factor,color2)
        end
        s = scatter(elem.local_nodes(:,1) + elem.node_deformations(:,1)*scale_factor,...
        elem.local_nodes(:,2) + elem.node_deformations(:,2)*scale_factor,...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
        s.SizeData = 20;
        % s.AlphaData = 0.5;
        alpha(0.3)
    end
    hold off
    axis padded
    axis equal
    title(sprintf('Deformed Mesh, Scale Factor = %.1f', scale_factor))