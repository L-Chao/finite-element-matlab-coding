% 3D truss example
clear;
E = 210000;
A = [100, 100, 100, 100];
nodeCoordinates = [
    4000, 4000, 3000;
    0, 4000, 0;
    0, 4000, 6000;
    4000, 0, 3000;
    8000, -1000, 1000];
elementNodes = [
    1, 2;
    1, 3;
    1, 4;
    1, 5];
numberElements = size(elementNodes, 1);
numberNodes = size(nodeCoordinates, 1);

xx = nodeCoordinates(1, :);
yy = nodeCoordinates(2, :);

GDof = 3*numberNodes;
U = zeros(GDof, 1);
force = zeros(GDof, 1);
force(2) = -10000;

stiffness = formStiffness3Dtruss(GDof, numberElements,...
    elementNodes, numberNodes, nodeCoordinates, E, A);

prescribedDof = (4:15)';
displacements = solution(GDof, prescribedDof, stiffness, force);

function stiffness = formStiffness3Dtruss(GDof, numberElements,...
    elementNodes, numberNodes, nodeCoordinates, E, A);
stiffness = zeros(GDof);
for e = 1:numberElements
    indice = elementNodes(e, :);
    elementDof = [3*indice(1) - 2, 3*indice(1) - 1, 3*indice(1)...
        3*indice(2) - 2, 3*indice(2) - 1, 3*indice(2)];
    x1 = nodeCoordinates(indice(1), 1);
    y1 = nodeCoordinates(indice(1), 2);
    z1 = nodeCoordinates(indice(1), 3);
    x2 = nodeCoordinates(indice(2), 1);
    y2 = nodeCoordinates(indice(2), 2);
    z2 = nodeCoordinates(indice(2), 3);
    L = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
    CXx = (x2 - x1)/L;
    CYx = (y2 - y1)/L;
    CZx = (z2 - z1)/L;
    T = [
        CXx*CXx, CXx*CYx, CXx*CZx;
        CYx*CXx, CYx*CYx, CYx*CZx;
        CZx*CXx, CZx*CYx, CZx*CZx];
    stiffness(elementDof, elementDof) =...
        stiffness(elementDof, elementDof) + E*A(e)/L*[T, -T; -T, T];    
end
end
function displacements=solution(GDof, prescribedDof, stiffness, force)
% function to find solution in terms of global displacements
activeDof=setdiff([1:GDof]', prescribedDof);
U=stiffness(activeDof,activeDof)\force(activeDof);
displacements=zeros(GDof,1);
displacements(activeDof)=U;
end