// Mesh.Algorithm = 8; 
gridsize = 6e-1;
ref_gridsize = 7.5e-4;
L = 4;
H = 4;
dx = 0.5;
dy = 0.0075;
//
Point(1) = { 0 , H/2  ,0.0, gridsize};
Point(2) = { 0 , -H/2  ,0.0, gridsize};
Point(3) = {  L , -H/2  ,0.0, gridsize};
Point(4) = {  L , H/2  ,0.0, gridsize};
// cuadrado donde estará la fractura
Point(5) = { 0.0 , dy ,0.0, ref_gridsize};
Point(6) = { 0.0 , -dy ,0.0, ref_gridsize};
Point(7) = { dx , -dy ,0.0, ref_gridsize};
Point(8) = { dx , dy ,0.0, ref_gridsize};
//
Line(1) = {1,5};
Line(2) = {6,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,1};
Line(10) = {5,8};
Line(11) = {8,7};
Line(12) = {7,6};
Line(13) = {6,5};
//
Curve Loop(1) = {1,10,11,12,2,3,4,5};
Curve Loop(2) = {-13, -12, -11, -10};
// EL curve loop 2 es mi crack inicial
Plane Surface(1) = {1};
Plane Surface(2) = {2};
// Colocacion de tags
Physical Line(1) = {4};
Physical Line(2) = {5};
Physical Line(3) = {1,13 ,2};
Physical Line(4) = {3};
Physical Surface(1) = {1, 2};
//
//Mesh.RecombineAll =1;
// Recombine Surface{1};
// Recombine Surface{2};