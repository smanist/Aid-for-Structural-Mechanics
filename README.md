Aid-for-Structural-Mechanics
============================

Mathematica package for variational operations, 1-D Hamilton's Principle, 1-D Lagrange's Equation and coefficient manipulations.

Example

<< StructuralMechanics`

T = 1/2 m wt^2 + \[Lambda\] (M/2 wt^2);
U = 1/2 EI wxx^2;
W = 1/2 P wx^2 + \[Lambda\] (-F w);
V = energyParser[{T, U, W}, {w}, x, {w[x, t]}];
MatrixForm@Transpose@hamiltonP[V, {w}, x]
phi = {(x/l)^2 (3 - x/l)/2, (x/l)^3 (3 - x/l)/2};
q = {q1[t], q2[t]};
V = energyParser[{T, U, W}, {w}, x, {phi.q}];
lagrangeP[V, {w}, x, {0, l}, q]

For energy expressions, terms like D[w[x,t],x,t] should be written as "wxt", so that the function can recognize that this term is 2nd mixed derivative of w[x,t]. t must be put as the last variable. Currently, the functions support up to wxxxxtttt. The \[Lambda\] should use the lambda symbol in Mathematica. It is used to differentiate concentrated terms (currently only supports terms at the end).

Description of each function is given in the package file.
