function plot_displacement_and_force(elem_struct,scale_factor,F,R,F_scale)
    N_elems = numel(elem_struct);
    if elem_struct(1).element_type == "beam"
        N_dof = 2;
    elseif elem_struct(1).element_type == "tri"
        N_dof = 2;
    end
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
            p = patch(elem.local_nodes(:,1),elem.local_nodes(:,2),color1);
            p.FaceAlpha = 0.3;

        end
        s = scatter(elem.local_nodes(:,1), elem.local_nodes(:,2),...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
        s.SizeData = 20;
        s.AlphaData = 0.3;
        % alpha(0.3)
    end
    % Deformed Mesh
    for i = 1:N_elems
        elem = elem_struct(i);
        if elem.element_type == "beam"
            plot(elem.local_nodes(:,1) + elem.node_deformations(:,1)*scale_factor,... 
            elem.local_nodes(:,2) + elem.node_deformations(:,2)*scale_factor,"Color",color2)
        elseif elem.element_type == "tri"
            p = patch(elem.local_nodes(:,1) + elem.node_deformations(:,1)*scale_factor,...
            elem.local_nodes(:,2) + elem.node_deformations(:,2)*scale_factor,color2);
            p.FaceAlpha = 0.3;
        end
        s = scatter(elem.local_nodes(:,1) + elem.node_deformations(:,1)*scale_factor,...
        elem.local_nodes(:,2) + elem.node_deformations(:,2)*scale_factor,...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
        s.SizeData = 20;
        s.AlphaData = 0.3;
        % alpha(0.3)
    end
    global_to_local = elem_struct.global_to_local;

    % Force and reaction quiver
    F_array = reshape(F,[N_dof,numel(F)/N_dof])';
    R_array = reshape(R,[N_dof,numel(R)/N_dof])';
    for i = 1:numel(F)/N_dof
        f = F_array(i,:);
        r = R_array(i,:);
        if ~isequal(f,[0,0])
        pos = global_to_local{i};
        quiver(pos(1),pos(2),f(1)*F_scale,f(2)*F_scale,"Color","r")
        end
        if ~isequal(r,[0,0])
        pos = global_to_local{i};
        quiver(pos(1),pos(2),r(1)*F_scale,r(2)*F_scale,"Color","b")
        end
    end
    hold off
    axis padded
    axis equal
    title(sprintf('Deformed Mesh, Scale Factor = %.1f', scale_factor))