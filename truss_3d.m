% Matlab 有限元结构动力学分析与工程应用, P74 truss examole
clear;
E = 210e6;
A = 0.02;
density = 7.3e3;

node_coordinate = [
    0, 0, 0;
    0, 0, 4;
    4, 0, 4;
    4, 0, 0;
    0, 5, 0;
    0, 5, 4;
    4, 5, 4;
    4, 5, 0;];
element_node = [
    1, 5;
    2, 6;
    3, 7;
    4, 8;
    5, 6;
    6, 7;
    7, 8;
    8, 5;];
node_number = size(node_coordinate, 1);
element_number = size(element_node, 1);
constraint = [
    1, 1;
    1, 2;
    1, 3;
    2, 1;
    2, 2;
    2, 3;
    3, 1;
    3, 2;
    3, 3;
    4, 1;
    4, 2;
    4, 3;];

% element_dispalcement = ones(node_number, 2);
% constraint = [1, 1; 1, 2; 3, 2];
% for i = 1:size(constraint, 1)
%     element_dispalcement(constraint(i, 1), constraint(i, 2)) = 0;
% end

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
        2*T, T;
        T, 2*T;];
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
d = d(:, indexf);