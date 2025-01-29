// Mesh.Algorithm = 8; 
gridsize = 10e-1;
ref_gridsize = 1e-3;
L = 10;
H = 10;
x_crack = 0.0;
dx = 0.1;
dy = 0.005;
y_crack = 0.0;
//
Point(1) = { -L/2 , H/2   ,0.0, gridsize};
Point(2) = { -L/2 , -H/2  ,0.0, gridsize};
Point(3) = {  L/2 , -H/2  ,0.0, gridsize};
Point(4) = {  L/2 , H/2  ,0.0,  gridsize};
// cuadrado donde estará la fractura
Point(11) = { x_crack - dx , y_crack - dy ,0.0, ref_gridsize};
Point(12) = { x_crack + dx , y_crack - dy ,0.0, ref_gridsize};
Point(13) = { x_crack + dx , y_crack + dy ,0.0, ref_gridsize};
Point(14) = { x_crack - dx , y_crack + dy ,0.0, ref_gridsize};
Point(15) = { x_crack - dx , 0.0 ,0.0, ref_gridsize};
Point(16) = { x_crack + dx , 0.0 ,0.0, ref_gridsize};
//
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
//
Line(10) = {11,12};
Line(11) = {12,16};
Line(12) = {16,15};
Line(13) = {15,11};
Line(14) = {16,13};
Line(15) = {13,14};
Line(16) = {14,15};
//
Curve Loop(1) = {1,2,3,4};
Curve Loop(2) = {14,15,16,-12};
Curve Loop(3) = {10,11,12,13};
// EL curve loop 2 es mi crack inicial
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(1) = {1,2,3};
// Colocacion de tags
Physical Line(10) = {3};
Physical Line(20) = {2};
Physical Line(30) = {1};
Physical Line(40) = {4};
Physical Surface(1) = {2,3,1};
//
//Mesh.RecombineAll =1;
// Recombine Surface{1};
// Recombine Surface{2};