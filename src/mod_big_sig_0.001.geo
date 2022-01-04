SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthMin =0.5;
Mesh.CharacteristicLengthMax = 15000;

Box(1)={-50000,-50000,-50000,100000,100000,50000};
Box(2)={-50000,-50000,     0,100000,100000,50000};

Coherence;

Physical Volume ( "air",11 ) = {1};
Physical Volume ( "layer",12 ) = {2};

//line source
l1=newreg; Circle(l1) = {0,0,5,500};
Curve{l1} In Volume {2};
Physical Curve ( "lsource1",51 ) = {l1};

// receivers
p71=newp;Point(p71)={0,0,5};
Coherence;

// field - source
Field[1] = Distance;
Field[1].NNodesByEdge = 300;
Field[1].EdgesList =  {l1};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 10;
Field[2].LcMax =15000;
Field[2].DistMin = 0;
Field[2].DistMax = 50000;

// field - receivers
Field[3] = Distance;
Field[3].NodesList = {p71};
Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 0.5;
Field[4].LcMax =15000;
Field[4].DistMin = 0;
Field[4].DistMax = 50000;

Field[11] = Min;
Field[11].FieldsList = {2,4};
Background Field = 11;