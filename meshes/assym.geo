// Mesh.Algorithm = 8; 
gridsize = 5e-3;
ref_gridsize = 6e-4;
L = 0.74/2;
x_crack = 0.0;
Hmax = 0.07;
Hmin = 0;
dx = 0.025;
dy = 0.01;
y_crack = 0.03;
distorsion = 0.019;
//
Point(1) = { -L , Hmax  ,0.0, gridsize};
Point(2) = { -L , Hmin  ,0.0, gridsize};
Point(3) = {  L , Hmin  ,0.0, gridsize};
Point(4) = {  L , Hmax  ,0.0, gridsize};
// cuadrado donde estará la fractura
Point(5) = { x_crack - dx , y_crack - 0.01 ,0.0, ref_gridsize};
Point(6) = { x_crack + dx*0.5 , y_crack - 0.01 - distorsion ,0.0, ref_gridsize};
Point(7) = { x_crack + dx*0.5 , y_crack + dy - distorsion,0.0   , ref_gridsize};
Point(8) = { x_crack - dx , y_crack + dy ,  0.0                 , ref_gridsize};
Point(9) = { x_crack - dx , y_crack      ,   0.0                , ref_gridsize};
Point(10) = { x_crack + dx*0.5 , y_crack - distorsion ,0.0, ref_gridsize};
//
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
//
Line(10) = {5,6};
Line(11) = {6,10};
Line(12) = {10,9};
Line(13) = {9,5};
Line(14) = {10,7};
Line(15) = {7,8};
Line(16) = {8,9};
//
Curve Loop(1) = {1,2,3,4};
Curve Loop(2) = {14,15,16,-12};
Curve Loop(3) = {10,11,12,13};
// EL curve loop 2 es mi crack inicial
Plane Surface(1) = {1,2,3};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
// Colocacion de tags
Physical Line(10) = {3};
Physical Line(20) = {2};
Physical Line(30) = {1};
Physical Surface(1) = {1,2,3};
//
//Mesh.RecombineAll =1;
// Recombine Surface{1};
// Recombine Surface{2};