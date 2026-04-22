function plot_elements(elem_struct)
    N_elems = numel(elem_struct);
    figure()
    hold on
    colors = parula(N_elems);
    colormap(parula(N_elems))
    clim([0.5 N_elems + 0.5])
    for i = 1:N_elems
        elem = elem_struct(i);
        
        if elem.element_type == "beam"
            plot(elem.local_nodes(:,1), elem.local_nodes(:,2),"Color",colors(i,:))
        elseif elem.element_type == "tri"
            
            patch(elem.local_nodes(:,1),elem.local_nodes(:,2),colors(i,:))
        end
        [center_x, center_y] = deal(mean(elem.local_nodes(:,1)), mean(elem.local_nodes(:,2)));
        text(center_x, center_y, string(i),'FontWeight','bold','HorizontalAlignment','center')
        s = scatter(elem.local_nodes(:,1), elem.local_nodes(:,2),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor','k');
        s.SizeData = 200;
        for j = 1:numel(elem.global_nodes)
            text(elem.local_nodes(j,1),elem.local_nodes(j,2), string(elem.global_nodes(j)),'FontWeight','bold','HorizontalAlignment','center','Color','b')
        end
        cb = colorbar('Ticks',1:N_elems);
        cb.Label.String = 'Element Number';
        
    end
    hold off
    axis padded
    axis equal
    title('Elements')