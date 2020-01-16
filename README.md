# Matlab codes for Finite Element Analysis

# 1. Basic steps of  finite element problem

1. Define a set of elements connected at nodes
2. For each element, compute **stiffness matrix**  $\textbf{K^e}$, and **force vector**  $\textbf{f^e}$  
3. Assemble the contribution of all elements into the global system $\textbf{Ka=f}$
4. Modify the global system by imposing essential (displacements) boundary conditions
5. Solve the global system and obtain the global displacements $\textbf{a}$
6. For each element, evaluate the strains and stresses (post-processing)

