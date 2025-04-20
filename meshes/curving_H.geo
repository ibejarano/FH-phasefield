// Mesh.Algorithm = 6; 
gridsize = 4e-1;
ref_gridsize = 4e-3;
L = 4;
H_sup = 0.06;
H_inf = 2;
dx = 0.1;
y_crack = -0.02;
dy_crack = 0.04 - y_crack;
ref_gridsize2 = ref_gridsize*(H_sup/0.025);
//
Point(1) = { -L/2 , H_sup   ,0.0, gridsize};
Point(2) = { -L/2 , -H_inf  ,0.0, gridsize};
Point(3) = {  L/2 , -H_inf  ,0.0, gridsize};
Point(4) = {  L/2 , H_sup  ,0.0,  gridsize};
Point(5) = {  dx , H_sup  ,0.0, ref_gridsize2};
Point(6) = {  -dx , H_sup  ,0.0,  ref_gridsize2};
// cuadrado donde estará la fractura
Point(11) = { -dx , y_crack, 0.0, ref_gridsize};
Point(12) = {  dx , y_crack, 0.0, ref_gridsize};
Point(13) = {  dx , y_crack + dy_crack, 0.0, ref_gridsize};
Point(14) = { -dx , y_crack + dy_crack, 0.0, ref_gridsize};
//
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
//
Line(10) = {14,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13, 14};
//
Curve Loop(1) = {1,2,3,4,5,6};
Curve Loop(2) = {10,11,12,13};
// EL curve loop 2 es mi crack inicial
Plane Surface(1) = {1,2};
Plane Surface(2) = {2};
// Colocacion de tags
Physical Line(10) = {3};
Physical Line(20) = {2};
Physical Line(30) = {1};
Physical Line(40) = {4,5,6};
Physical Surface(1) = {2,1};
//
//Mesh.RecombineAll =1;
// Recombine Surface{1};
// Recombine Surface{2};