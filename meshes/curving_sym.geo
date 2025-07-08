// Mesh.Algorithm = 6; 
gridsize = 8.0;
ref_gridsize = 0.001;
L = 10;
H_sup = 3;
H_inf = 9.7;
dx = 0.025;
dy_crack =0.02;
refb_gridsize = 0.8*H_sup/4;
//
Point(1) = { 0 , H_sup   ,0.0, refb_gridsize};
Point(2) = { 0 , -H_inf  ,0.0, gridsize};
Point(3) = {  L , -H_inf  ,0.0, gridsize};
Point(4) = {  L , H_sup  ,0.0,  gridsize};
Point(5) = {  dx*3 , H_sup  ,0.0, refb_gridsize};
// cuadrado donde estará la fractura
Point(11) = { 0 , -dy_crack/2, 0.0, ref_gridsize};
Point(12) = { dx, -dy_crack/2, 0.0, ref_gridsize};
Point(13) = { dx, dy_crack/2, 0.0, ref_gridsize};
Point(14) = { 0 , dy_crack/2, 0.0, ref_gridsize};
Point(15) = { dx*4, dy_crack*10, 0.0, ref_gridsize};
Point(16) = { dx*2.5 , -dy_crack/2, 0.0, ref_gridsize};
//
Line(1) = {1,14};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};
Line(6) = {11,2};
Line(7) = {12,16};
Line(8) = {16,15};
Line(9) = {15,13};
//
Line(10) = {14,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13, 14};
//
Curve Loop(1) = {1,-13,-9,-8,-7,-11,6,2,3,4,5};
Curve Loop(2) = {10,11,7,8,9,13};
// Curve Loop(3) = {7,8,9,-12};
// EL curve loop 2 es mi crack inicial
Plane Surface(1) = {1,2};
Plane Surface(2) = {2};
//
//Mesh.RecombineAll =1;
// Recombine Surface{1};
// Recombine Surface{2};