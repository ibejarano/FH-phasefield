L = 100.0;
H = 50.0;
H_sup = 0.2;
a = 0.1;
h_entalla = 0.001;

mesh_size_fractura= 1e-3;
mesh_size_gral = 2;
// Puntos del rectángulo
Point(1) = {-L/2, -H, 0, mesh_size_gral};
Point(2) = {L/2, -H, 0, mesh_size_gral};
Point(3) = {L/2, H_sup, 0, mesh_size_gral};
Point(4) = {-L/2, H_sup, 0, mesh_size_gral};

// Puntos de la entalla
Point(5) = {-a/2, -h_entalla, 0, mesh_size_fractura};
Point(6) = {a/2, -h_entalla, 0, mesh_size_fractura};
Point(7) = {a/2, h_entalla, 0, mesh_size_fractura};
Point(8) = {-a/2, h_entalla, 0, mesh_size_fractura};
Point(9) = {a*1.005/2, 0.0, 0, mesh_size_fractura};
Point(10) = {-a*1.02/2, 0.0, 0, mesh_size_fractura};

// Líneas
Line(1) = {1, 2}; // Inferior
Line(2) = {2, 3}; // Derecha
Line(3) = {3, 4}; // Superior
Line(4) = {4, 1}; // Izquierda

// Líneas Fractura
Line(5) = {8, 7}; // Inferior
Line(6) = {7, 9}; // Derecha
Line(7) = {9, 6}; // Derecha
Line(8) = {6, 5}; // Derecha
Line(9) = {5, 10}; // Izquierda
Line(10) = {10, 8}; // Izquierda

// Loop y superficie
Line Loop(10) = {1, 2, 3, 4};
Line Loop(11) = {5, 6, 7, 8, 9, 10};
Plane Surface(6) = {10, 11};

// Opciones de mallado 
// —– Campos de refinamiento —–
// Distancia a los puntos/líneas de la entalla
Field[1] = Distance;
Field[1].NodesList = {5, 6, 7, 8, 9, 10};
Field[1].EdgesList = {5, 6, 7, 8, 9, 10};

// Umbral de tamaños en función de la distancia
Field[2] = Threshold;
Field[2].IField      = 1;     // usa el Distance
Field[2].LcMin       = 0.001;  // tamaño mínimo en la entalla
Field[2].LcMax       = 1;   // tamaño “global” lejos de la entalla
Field[2].DistMin     = 0.005;     // a menos de 0 unidades → LcMin
Field[2].DistMax     = 1.2;   // a más de 0.2 unidades → LcMax

// Aplico el campo como fondo para la malla
Background Field = 2;

// —– Parámetros de generación y optimización —–
Mesh.Format = 1;            // MSH2
Mesh.Optimize = 1;          // activa smoothing
Mesh.OptimizeNetgen = 1;    // optimizador alternativo
Mesh.Algorithm = 6;         // por ejemplo: Delaunay


// Condiciones de borde (Physical Lines)
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4, 8};
Physical Line(5) = {5, 6, 10};
Physical Line(6) = {7, 8, 9};

// Superficie física
Physical Surface(1) = {6};
