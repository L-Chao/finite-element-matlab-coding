function model_shape(mode, element_node, node_coordinate, element_displacement, v)
% v : model shape matrix
figure()
title(['mode = ', num2str(mode)])
hold on
for i = 1:size(element_node, 1)
    x = [node_coordinate(element_node(i, 1), 1), node_coordinate(element_node(i, 2), 1)];
    y = [node_coordinate(element_node(i, 1), 2), node_coordinate(element_node(i, 2), 2)];
    z = [node_coordinate(element_node(i, 1), 3), node_coordinate(element_node(i, 2), 3)];
    plot3(x, y, z, 'k-')
    mx = [0; 0];
    my = [0; 0];
    mz = [0; 0];
    for j = 1:2
        if element_displacement(element_node(i, j), 1) ~= 0
            mx(j) = x(j) + v(element_displacement(element_node(i, j), 1), mode)*10;
        else
            mx(j) = x(j);
        end
        if element_displacement(element_node(i, j), 2) ~=0
            my(j) = y(j) + v(element_displacement(element_node(i, j), 2), mode)*10;
        else
            my(j) = y(j);
        end
        if element_displacement(element_node(i, j), 3) ~=0
            mz(j) = z(j) + v(element_displacement(element_node(i, j), 3), mode)*10;
        else
            mz(j) = z(j);
        end
    end
    if i == 10
        plot3(mx, my, mz, 'g--')
    else
        plot3(mx, my, mz, 'r--')
    end
    view(3);
end
    hold off
end
