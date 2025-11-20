//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2, 0, 0, 1.0};
//+
Point(3) = {0, -0.2, 0, 1.0};
//+
Point(4) = {2, -0.2, 0, 1.0};
//+
Point(5) = {0, 1.0, 0, 1.0};
//+
Point(6) = {2, 1.0, 0, 1.0};
//+
Point(7) = {0, 1.2, 0, 1.0};
//+
Point(8) = {2, 1.2, 0, 1.0};
//+
Line(1) = {3, 4};
//+
Line(2) = {4, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 8};
//+
Line(6) = {8, 7};
//+
Line(7) = {7, 5};
//+
Line(8) = {1, 3};
//+
Curve Loop(1) = {3, 8, 1, 2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 4, 5, 6};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {1, 3, 4, 6} = 21 Using Progression 1;
//+
Transfinite Curve {8, 2, 5, 7} = 6 Using Progression 1;
//+
Transfinite Surface {1} = {3, 4, 2, 1};
//+
Transfinite Surface {2} = {5, 6, 8, 7};
//+
Recombine Surface {1};
//+
Recombine Surface {2};
//+
Extrude {0, 0, 0.1} {
  Surface{1}; Surface{2}; Layers {3}; Recombine;
}
//+
Physical Surface("p1", 11) = {30, 52};
//+
Physical Surface("p2", 12) = {2, 1};
//+
Physical Surface("walls", 4) = {43, 17};
//+
Physical Surface("bottomWall", 13) = {25};
//+
Physical Surface("topWall", 14) = {51};
//+
Physical Volume("solid", 53) = {2, 1};
//+
Physical Surface("sideWalls", 15) = {39, 21, 29, 47};
