eps:=1.*10^-6;
LoadPackage("GAPic");
Read("helper.gi");
Read("SelfIntersectingComplexes.gd");
Read("detection.gi");
Read("retriangulation.gi");
Read("chambers.gi");
Read("main.gi");
Read("nonmanifold.gi");
coordinates:=[ [ 0, 0, 1 ], [ 0, 0, 0 ], [ 0, 1, 0 ], [ 0, 1, 1 ], [ 1, 1, 1 ], [ 1, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 1 ] ];

faces:=[
    [1, 3, 2],
    [3, 4, 1],
    [7, 5, 8],
    [5, 6, 7],
    [1, 2, 7],
    [1, 8, 7],
    [3, 4, 5],
    [3, 6, 5],
    [7, 2, 3],
    [7, 3, 6],
    [1, 4, 8],
    [4, 5, 8]
];

cube:=SimplicialSurfaceByVerticesInFaces(faces);
faces_two:=Concatenation(faces,faces+8);
coordinates_two:=Concatenation(coordinates,coordinates+[0.5,0.5,0.])*1.;
cubes_two:=SimplicialSurfaceByVerticesInFaces(faces_two);
pr:=SetVertexCoordinates3D(cubes_two,coordinates_two);
DrawComplexToJavaScript(cubes_two,"Cubes/TwoCubesXY",pr);
eps:=1.*10^-6;
r := PrintableOuterHull(cubes_two,coordinates_two,"Cubes/TwoCubes",eps,0.1,Group(()),true,false);
surf := r[1];
pr2 := SetVertexCoordinates3D(surf,r[3]);
DrawComplexToJavaScript(surf,"Cubes/TwoCubesXY_intersc_free",pr2);

coordinates_two:=Concatenation(coordinates,coordinates+[0.5,0.5,0.5])*1.;
cubes_two:=SimplicialSurfaceByVerticesInFaces(faces_two);
pr:=SetVertexCoordinates3D(cubes_two,coordinates_two);
DrawComplexToJavaScript(cubes_two,"Cubes/TwoCubesXYZ",pr);
eps:=1.*10^-6;
PrintableOuterHull(cubes_two,coordinates_two,"Cubes/TwoCubesXYZ",eps,0.1,Group(()),true,false);


extra_points:=11;

coordinates:=[ [ 0, 0, 1 ], [ 0, 0, 0 ], [ 0, 1, 0 ], [ 0, 1, 1 ], [ 1, 1, 1 ], [ 1, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 1 ] ];

for i in [1..extra_points] do
    Add(coordinates,coordinates[5]+1.*i/(1+extra_points)*(coordinates[6]-coordinates[5]));
od;

faces:=[
    [1, 3, 2],
    [3, 4, 1],
    [1, 2, 7],
    [1, 8, 7],
    [7, 2, 3],
    [7, 3, 6],
    [1, 4, 8],
    [4, 5, 8]
];

Add(faces,[3,6,8+extra_points]);
Add(faces,[4,5,8+1]);
for i in [1..(extra_points-1)/2] do
    Add(faces,[3,8+extra_points-i,8+extra_points+1-i]);
    Add(faces,[4,8+i,8+i+1]);
od;
Add(faces,[3,4,(extra_points-1)/2+1+8]);

Add(faces,[7,6,8+extra_points]);
Add(faces,[8,5,8+1]);
for i in [1..(extra_points-1)/2] do
    Add(faces,[7,8+extra_points-i,8+extra_points+1-i]);
    Add(faces,[8,8+i,8+i+1]);
od;
Add(faces,[7,8,(extra_points-1)/2+1+8]);


faces_two:=Concatenation(faces,faces+Size(coordinates));

coordinates_two:=Concatenation(coordinates-[1,1,0],(coordinates-[1,1,0])*DiagonalMat([-1,-1,1]))*1.;
cubes_two:=SimplicialSurfaceByVerticesInFaces(faces_two);
pr:=SetVertexCoordinates3D(cubes_two,coordinates_two);

eps:=1.*10^-6;
data:=PrintableOuterHull(cubes_two,coordinates_two,"Cubes/TwoCubesEdge",eps,0.1,Group(()),true,false);;
outer_triangle:=FindOuterTriangle(data[1],data[3]);
data_outer:=OuterHull(data[1],data[3],outer_triangle[1],outer_triangle[2]);
DrawComplexToJavaScript(cubes_two,"Cubes/TwoCubesEdge",pr);

shift_param:=0.04;
f := RemedyNonManifold(data_outer,data[3],shift_param);
DrawSTLwithNormals(f[1],Concatenation("Cubes/TwoCubesEdgeFixed",String(extra_points)),f[2],data_outer[4],[]);



# two cubes 0 extra points
coordinates:=1.*[ [ 0, 0, 1 ], [ 0, 0, 0 ], [ 0, 1, 0 ], [ 0, 1, 1 ], [ 1, 1, 1 ], [ 1, 1, 0 ], [ 1, 0, 0 ], [ 1, 0, 1 ], [ 0, 0, -1], [ 0, 1, -1], [-1, 1, -1], [-1, 0, -1], [-1, 1, 0], [-1, 0, 0] ];
faces:=[
    [1, 3, 2],
    [3, 4, 1],
    [1, 2, 7],
    [1, 8, 7],
    [7, 2, 3],
    [7, 3, 6],
    [1, 4, 8],
    [4, 5, 8],
    [3, 4, 5],
    [3, 5, 6],
    [5, 6, 7],
    [5, 7, 8],
    [2, 3, 9],
    [3, 9, 10],
    [3, 10, 13],
    [10, 11, 13],
    [11, 13, 14],
    [11, 12, 14],
    [9, 12, 14],
    [2, 14, 9],
    [2, 3, 14],
    [3, 13, 14],
    [9, 10, 12],
    [10, 11, 12]
];
cubes_two:=TriangularComplexByVerticesInFaces(faces);
pr:=SetVertexCoordinates3D(cubes_two,coordinates);
data:=PrintableOuterHull(cubes_two,coordinates,"Cubes/TwoCubesEdge",eps,0.1,Group(()),true,false);;
outer_triangle:=FindOuterTriangle(data[1],data[3]);
data_outer:=OuterHull(data[1],data[3],outer_triangle[1],outer_triangle[2]);
DrawComplexToJavaScript(cubes_two,"Cubes/TwoCubesEdge0",pr);

shift_param:=0.04;
f := RemedyNonManifold(data_outer,data[3],shift_param);
DrawSTLwithNormals(f[1],"test",f[2],data_outer[4],[]);




    
