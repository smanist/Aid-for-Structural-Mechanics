Aid-for-Structural-Mechanics
============================

Mathematica package for variational operations, 1-D Hamilton's Principle, 1-D Lagrange's Equation and coefficient manipulations.
Some analytical and numerical methods are included.

For energy expressions, terms like D[w[x,t],x,t] should be written as "wxt", so that the function can recognize that this term is 2nd mixed derivative of w[x,t]. t must be put as the last variable. Currently, the functions support up to wxxxxtttt. The \[Lambda\] should use the lambda symbol in Mathematica. It is used to differentiate concentrated terms (currently only supports terms at the end).

Description of each function is given in the package file. A sample file is also provided.
