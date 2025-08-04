// Define points of the rectangle
Lx = 1;
Ly = 2;
Lcrack = 0.3;
beta_angle = 0;
beta_rad = beta_angle * Pi/180 ;
x1 = Cos(beta_rad) * Lcrack ;
y1 = Sin(beta_rad) * Lcrack ;

lc = 0.2;
dist_min = 0.02;
dist_max = 0.3;
refinement_ratio = 400;

Point(1) = {-Lx, -Ly, 0, lc};
Point(2) = {Lx, -Ly, 0, lc};
Point(3) = {Lx, Ly, 0, lc};
Point(4) = {-Lx, Ly, 0, lc};


Point(11) = {-Lx, 0, 0, lc};
Point(12) = {-x1, -y1, 0, lc};
Point(13) = {x1, y1, 0, lc};
Point(14) = {Lx, 0, 0, lc};

// Define bottom rectangle
Line(1) = {1, 2};
Line(2) = {2, 14};
Line(3) = {14, 13};

// Linea bot crack
Line(4) = {13, 12};

Line(5) = {12, 11};
Line(6) = {11, 1};

// Definir top rectangle (provisorio)
Line(12) = {14, 3};
Line(13) = {3, 4};
Line(14) = {4, 11};

// Linea top crack
Line(15) = {12, 13};

// Loop bottom rect
Curve Loop(1) = {1, 2, 3, 4, 5, 6};

// Loop top rect
Curve Loop(2) = {-5, 15, -3 , 12, 13, 14};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Distancia a los puntos/líneas de la entalla
Field[1] = Distance;
Field[1].NodesList = {12, 13};

// Umbral de tamaños en función de la distancia
Field[2] = Threshold;
Field[2].IField      = 1;     // usa el Distance
Field[2].LcMin       = lc/refinement_ratio;  // tamaño mínimo en la entalla
Field[2].LcMax       = lc;   // tamaño “global” lejos de la entalla
Field[2].DistMin     = dist_min;     // a menos de 0 unidades → LcMin
Field[2].DistMax     = dist_max;   // a más de 0.2 unidades → LcMax

// Aplico el campo como fondo para la malla
Background Field = 2;

// —– Parámetros de generación y optimización —–
Mesh.Format = 1;            // MSH2
Mesh.Optimize = 1;          // activa smoothing
Mesh.OptimizeNetgen = 1;    // optimizador alternativo
Mesh.Algorithm = 6;         // por ejemplo: Delaunay

// lineas fisicas
Physical Line(1) = {14, 6};
Physical Line(2) = {1};
Physical Line(3) = {2, 12};
Physical Line(4) = {13};
Physical Line(10) = {15};
Physical Line(11) = {4};

// Superficie física
Physical Surface(1) = {1, 2};