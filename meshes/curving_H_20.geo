// Mesh.Algorithm = 6; 
gridsize = 10;
ref_gridsize = 1e-3;
L = 200;
H_sup = 0.2;
H_inf = 100;
dx = 0.3;
y_crack = -0.01;
dy_crack = 0.165 - y_crack;
point_surface = ref_gridsize*(H_sup/0.005);
offx = 0.075;
offy = ref_gridsize*20;
adhocx = 0.08;
//
Point(1) = { -L/2 , H_sup   ,0.0, gridsize};
Point(2) = { -L/2 , -H_inf  ,0.0, gridsize};
Point(3) = {  L/2 , -H_inf  ,0.0, gridsize};
Point(4) = {  L/2 , H_sup  ,0.0,  gridsize};
Point(5) = {  dx , H_sup  ,0.0, point_surface/2};
Point(6) = {  -dx , H_sup  ,0.0,  point_surface/2};
// cuadrado minimizado de la fractura
Point(11) = { -dx , y_crack, 0.0, ref_gridsize};
Point(12) = {  dx, y_crack, 0.0, ref_gridsize};
Point(13) = {  dx+adhocx , y_crack + dy_crack, 0.0, ref_gridsize};
Point(14) = { -(dx+adhocx) , y_crack + dy_crack, 0.0, ref_gridsize};
Point(15) = { -offx , y_crack + offy, 0.0, ref_gridsize};
Point(16) = {  offx , y_crack + offy, 0.0, ref_gridsize};
Point(17) = {  offx*2 + adhocx , y_crack + dy_crack, 0.0, ref_gridsize};
Point(18) = { -(offx*2 + adhocx) , y_crack + dy_crack, 0.0, ref_gridsize};
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
Line(13) = {13, 17};
Line(14) = {17, 16};
Line(15) = {16, 15};
Line(16) = {15, 18};
Line(17) = {18, 14};
//
Curve Loop(1) = {1,2,3,4,5,6};
Curve Loop(2) = {10,11,12,13,14,15,16,17};
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