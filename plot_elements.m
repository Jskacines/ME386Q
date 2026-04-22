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
            scatter(elem.local_nodes(:,1), elem.local_nodes(:,2),'k','filled')
            plot(elem.local_nodes(:,1), elem.local_nodes(:,2),"Color",colors(i,:))
        elseif elem.element_type == "tri"
            scatter(elem.local_nodes(:,1), elem.local_nodes(:,2),'k','filled')
            patch(elem.local_nodes(:,1),elem.local_nodes(:,2),colors(i,:))
        end
        [center_x, center_y] = deal(mean(elem.local_nodes(:,1)), mean(elem.local_nodes(:,2)));
        text(center_x, center_y, string(i),'FontWeight','bold','HorizontalAlignment','center')
        cb = colorbar('Ticks',1:N_elems);
        cb.Label.String = 'Element Number';
    end
    hold off
    axis padded
    axis equal
    title('Elements')