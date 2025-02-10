Mesh.Algorithm = 6; 
h_coarse = 10e0;
h = 0.005;
L = 100;
H = 75*h;
W = 10*h;
//
Point(1) = { -L/2 , L/2   ,0.0, h_coarse};
Point(2) = { -L/2 , -L/2  ,0.0, h_coarse};
Point(3) = {  L/2 , -L/2  ,0.0, h_coarse};
Point(4) = {  L/2 , L/2  ,0.0,  h_coarse};
// cuadrado donde estará la fractura
Point(11) = { -H   ,   W/2 ,       0.0, h};
Point(12) = { -H   ,  0.0 ,          0.0, h};
Point(13) = { -H , -W/2 ,0.0, h};
Point(14) = { H  , -W/2 ,0.0, h};
Point(15) = { H  , 0.0 ,0.0, h};
Point(16) = { H  , W/2 ,0.0, h};

//
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
//
Line(10) = {11,12};
Line(11) = {12,13};
Line(12) = {13,14};
Line(13) = {14,15};
Line(14) = {15,16};
Line(15) = {16,11};
Line(16) = {12,15};
//
Curve Loop(1) = {1,2,3,4};
Curve Loop(2) = {10,16,14,15};
Curve Loop(3) = {11,12,13,-16};
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