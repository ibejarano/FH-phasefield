SetFactory("OpenCASCADE");

// --- Parámetros (ajustá a gusto)
lc = 0.40;        // tamaño general de malla
lc_path = 4e-3;   // tamaño fino cerca de la trayectoria
rin = 4e-2;       // radio interior: hasta aquí aplica lc_path
rout = 0.50;      // radio exterior: desde aquí vuelve a lc

Lx = 10; 
Ly = 10;
H = 0.5;

// --- Geometría: rectángulo
Point(1) = {0,  -Ly, 0, lc};
Point(2) = {Lx, -Ly, 0, lc};
Point(3) = {Lx, H, 0, lc};
Point(4) = {0,  H, 0, lc};
Line(11) = {1,2}; Line(12) = {2,3}; Line(13) = {3,4}; Line(14) = {4,1};
Curve Loop(21) = {11,12,13,14};
Plane Surface(31) = {21};

// --- Trayectoria (ejemplo: una spline dentro del dominio)
Point(101) = {0, 0.0, 0, lc};
Point(110) = {0.1, 0.0, 0, lc};
Point(111) = {0.3, 0.0, 0, lc};
Point(102) = {0.5, 0.05, 0, lc};
Point(103) = {0.81, 0.15, 0, lc};
Point(104) = {1.28, 0.38, 0, lc};
Spline(201) = {101,110, 111,102,103,104}; // <- ID de la trayectoria

// --- Campo de tamaño: distancia a la curva + umbral
Field[1] = Distance;
Field[1].CurvesList = {201};       // trayectoria
Field[1].NumPointsPerCurve = 100;  // muestreo (más = distancia más precisa)

Field[2] = Threshold;
Field[2].InField = 1;       // toma distancias del Field[1]
Field[2].SizeMin = lc_path; // tamaño fino dentro de rin
Field[2].SizeMax = lc;      // tamaño general más allá de rout
Field[2].DistMin = rin;     // transición empieza aquí
Field[2].DistMax = rout;    // y termina aquí

Background Field = 2;       // aplica el campo como tamaño de malla global

// (opcional) límites duros globales
Mesh.MeshSizeMin = lc_path;
Mesh.MeshSizeMax = lc;

