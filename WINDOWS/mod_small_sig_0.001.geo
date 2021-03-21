SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthMin = 0.01;
Mesh.CharacteristicLengthMax = 1000;

Box(1)={-5000,-5000,-5000,10000,10000,5000};
Box(2)={-5000,-5000,    0,10000,10000,5000};

Coherence;

Physical Volume ( "air",11 ) = {1};
Physical Volume ( "layer",12 ) = {2};

//line source
l1=newreg; Circle(l1) = {0,0,1,2};
Curve{l1} In Volume {2};
Physical Curve ( "lsource1",51 ) = {l1};

// receivers
p71=newp;Point(p71)={0,0,1};
Coherence;

// field - source
Field[1] = Distance;
Field[1].NNodesByEdge = 60;
Field[1].EdgesList =  {l1};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.2;
Field[2].LcMax = 1000;
Field[2].DistMin = 0;
Field[2].DistMax = 5000;

// field - receivers
Field[3] = Distance;
Field[3].NodesList = {p71};
Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 0.01;
Field[4].LcMax = 1000;
Field[4].DistMin = 0;
Field[4].DistMax = 5000;

Field[11] = Min;
Field[11].FieldsList = {2,4};
Background Field = 11;