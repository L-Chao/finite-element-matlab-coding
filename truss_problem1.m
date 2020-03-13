% Matlab 有限元结构动力学分析与工程应用, P74 truss examole
clear;
E = 2.11e11;
A = 1e-4;
density = 7.3e3;
node_number = 5;
element_number = 7;
node_coordinate = [0, 0; 1, 0; 2, 0; 0, 1; 1, 1];
element_node = [1, 2; 2, 3; 1, 4; 2, 4; 2, 5; 3, 5; 4, 5];

element_dispalcement = ones(node_number, 2);
constraint = [1, 1; 1, 2; 3, 2];
for i = 1:size(constraint, 1)
    element_dispalcement(constraint(i, 1), constraint(i, 2)) = 0;
end

dof = sum(sum(element_dispalcement));
Gdof = 2*node_number;
% element stiffness matrix
K_e = E*A*[
    1, 0, -1, 0;
    0, 0, 0, 0;
    -1, 0, 1, 0;
    0, 0, 0, 0;];
% element mass matrix
% Concentrated mass matrix
M_e = density*A/2*eye(4);
% Consistent mass matrix
% M_e = density*A/4*[
%     2, 0, 1, 0;
%     0, 2, 0, 1;
%     1, 0, 2, 0;
%     0, 1, 0, 2];
% structure stiffness matrix
K_s = zeros(Gdof);
% structure mass matrix
M_s = zeros(Gdof);
for i = 1:element_number
    indice = element_node(i, :);
    element_dof = [2*indice(1) - 1, 2*indice(1), 2*indice(2) - 1, 2*indice(2)];
    x1 = node_coordinate(indice(1), 1);
    y1 = node_coordinate(indice(1), 2);
    x2 = node_coordinate(indice(2), 1);
    y2 = node_coordinate(indice(2), 2);
    element_length = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    C = (x2 - x1)/element_length;
    S = (y2 - y1)/element_length;
    lambda = [
        C, S;
        -S, C;];
    T = [
        lambda, zeros(2);
        zeros(2), lambda];
    K_s(element_dof, element_dof) = K_s(element_dof, element_dof) + T'*K_e*T/element_length;
    M_s(element_dof, element_dof) = M_s(element_dof, element_dof) + T'*M_e*T*element_length;
end
%%
% add constraint
prescribed_dof = 2*(constraint(:, 1) -1) + constraint(:, 2);
K_s(prescribed_dof, :) = [];
K_s(:, prescribed_dof) = [];
M_s(prescribed_dof, :) = [];
M_s(:, prescribed_dof) = [];
%%
[v, d] = eig(K_s, M_s);
% d: eigenvalue
% v: eigenvector
frequency = sqrt(diag(d))/(2*pi);
[frequency, indexf] = sort(frequency);
d = d(:, indexf);
%%
mode = 1;
for i = 1:element_number
    x = [node_coordinate(element_node(i, 1), 1), node_coordinate(element_node(i, 2), 1)];
    y = [node_coordinate(element_node(i, 1), 2), node_coordinate(element_node(i, 2), 2)];
    line(x, y, 'color', 'black')
    mx = [0; 0];
    my = [0; 0];
    for j = 1:2
        if element_dispalcement(element_node(i, j), 1) ~= 0
            mx(j) = x(j) + v(element_dispalcement(element_node(i, j), 1), mode)/8;
        else
            mx(j) = x(j);
        end
        if element_dispalcement(element_node(i, j), 2) ~=0
            my(j) = y(j) + v(element_dispalcement(element_node(i, j), 2), mode)/8;
        else
            my(j) = y(j);
        end
    end
    line(mx, my, 'color', 'black', 'linestyle', ':')
end
    
    
    