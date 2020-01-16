% MATLAB codes for Finite Element Analysis
% problem1.m
% antonio ferreira 2008
clear;

% ElementNodes: connections at elements
ElementNodes = [1, 2; 2, 3; 2, 4];

% NumberElements: number of Elements
NumberElements = size(ElementNodes, 1);

% NumberNodes: number of nodes
NumberNodes = 4;

% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
DispalceElements = zeros(NumberNodes, 1);
Force = zeros(NumberNodes, 1);
Stiffness = zeros(NumberNodes);

% applied load at node 2
Force(2) = 10.0;
% computation of the system stiffness matrix
for e = 1:NumberElements
    % elementDof: element degrees of freedom (Dof)
    ElementDof = ElementNodes(e, :);
    Stiffness(ElementDof, ElementDof) =...
        Stiffness(ElementDof, ElementDof) + [1, -1; -1, 1];
end

% boundary conditions and solution
% prescribed dofs
PrescribedDof = [1; 3; 4];
% free Dof: activeDof
ActiveDof = setdiff((1:NumberNodes)', PrescribedDof);

% solution
displacements = (Stiffness(ActiveDof, ActiveDof))^-1*Force(ActiveDof);

% positioning all displacements
displacements1 = zeros(NumberNodes, 1);
displacements1(ActiveDof) = displacements;


% reactions
F = Stiffness*displacements1;
reations = F(PrescribedDof);