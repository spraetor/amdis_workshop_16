// Gmsh project created on Thu Sep 17 15:30:35 2015
Point(1) = {0, 0, 0, 1.0, 1.0};
Point(2) = {0, 1, 0, 1.0, 1.0};
Point(3) = {10, 1, 0, 1.0, 1.0};
Point(4) = {10, 0, 0, 1.0, 1.0};

gridsize = 0.01;
Point(5) = {2, 0.4, 0, 1.0, gridsize};
Point(6) = {2, 0.6, 0, 1.0, gridsize};
Point(7) = {2, 0.2, 0, 1.0, gridsize};
Point(8) = {1.8, 0.4, 0, 1.0, gridsize};
Point(9) = {2.2, 0.4, 0, 1.0, gridsize};

Circle(1) = {6, 5, 9};
Circle(2) = {9, 5, 7};
Circle(3) = {7, 5, 8};
Circle(4) = {8, 5, 6};

Line(5) = {2, 3};
Line(6) = {3, 4};
Line(7) = {4, 1};
Line(8) = {1, 2};
Line Loop(9) = {5, 6, 7, 8};
Line Loop(10) = {4, 1, 2, 3};
Plane Surface(11) = {9, 10};
Physical Surface(12) = {11};
Physical Line(5) = {4, 1, 2, 3};
