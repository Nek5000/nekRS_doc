//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {2, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {1, 3} = 21 Using Progression 1;
//+
Transfinite Curve {2, 4} = 11 Using Progression 1;
//+
Transfinite Surface {1} = {1, 2, 4, 3};
//+
Recombine Surface {1};
//+
Extrude {0, 0, 0.1} {
  Surface{1}; Layers {3}; Recombine;
}
//+
Physical Volume("fluid", 27) = {1};
//+
Physical Surface("sideWalls", 3) = {13, 21};
//+
Physical Surface("walls", 4) = {17, 25};
//+
Physical Surface("p1", 1) = {1};
//+
Physical Surface("p2", 2) = {26};
