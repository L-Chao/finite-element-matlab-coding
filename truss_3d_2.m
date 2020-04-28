clear;
E = 6.89e10;

density = 2.77e3;

A = 2e-02;

node_coordinate = 1e-3*[
    952.5, 0, 5080;
    -952.5,0, 5080;
    -952.5, 952.5, 2540;
    952.5, 952.5, 2540;
    952.5, -952.5, 2540;
    -952.5, -952.5, 2540;
    -2540, 2540, 0;
    2540, 2540, 0;
    2540, -2540, 0;
    -2540, -2540, 0;];
element_node = [
    1, 2;
    1, 4;
    2, 3;
    1, 5;
    2, 6;
    2, 4;
    2, 5;
    1, 3;
    1, 6;
    6, 3;
    4, 5;
    3, 4;
    6, 5;
    3, 10;
    6, 7;
    4, 9;
    5, 8;
    4, 7;
    3, 8;
    5, 10;
    6, 9;
    6, 10;
    3, 7;
    4, 8;
    5, 9;];
node_number = size(node_coordinate, 1);
element_number = size(element_node, 1);
constraint = [
    7, 1;
    7, 2;
    7, 3;
    8, 1;
    8, 2;
    8, 3;
    9, 1;
    9, 2;
    9, 3;
    10, 1;
    10, 2;
    10, 3;];

element_displacement = ones(node_number, 3);
for i = 1:size(constraint, 1)
    element_displacement(constraint(i, 1), constraint(i, 2)) = 0;
end
dof = 0;
for i = 1:node_number
    for j = 1:3
        if element_displacement(i, j) ~= 0
            dof = dof + 1;
            element_displacement(i, j) = dof;
        end
    end
end
Gdof = 3*node_number;
% structure stiffness matrix
K_s = zeros(Gdof);
% structure mass matrix
M_s = zeros(Gdof);
for i = 1:element_number
    indice = element_node(i, :);
    element_dof = [3*indice(1) - 2, 3*indice(1) - 1, 3*indice(1), 3*indice(2) - 2, 3*indice(2) - 1, 3*indice(2)];
    x1 = node_coordinate(indice(1), 1);
    y1 = node_coordinate(indice(1), 2);
    z1 = node_coordinate(indice(1), 3);
    x2 = node_coordinate(indice(2), 1);
    y2 = node_coordinate(indice(2), 2);
    z2 = node_coordinate(indice(2), 3);
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
    l = (x2 - x1)/L;
    m = (y2 - y1)/L;
    n = (z2 - z1)/L;
    T = [
        l*l, l*m, l*n;
        l*m, m*m, m*n;
        l*n, m*n, n*n;];
    K_e = E*A/L*[
        T, -T;
        -T, T;];
    M_e = A*density*L/6*[
        2*eye(3), eye(3);
        eye(3), 2*eye(3);];
    K_s(element_dof, element_dof) = K_s(element_dof, element_dof) + K_e;
    M_s(element_dof, element_dof) = M_s(element_dof, element_dof) + M_e;
end
%%
% add constraint
prescribed_dof = 3*(constraint(:, 1) -1) + constraint(:, 2);
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
v = v(:, indexf);
%%
mode = 1;
model_shape(mode, element_node, node_coordinate, element_displacement, v)