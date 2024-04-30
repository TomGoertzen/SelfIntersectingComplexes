LoadPackage("GAPic");
LoadPackage("SelfIntersectingComplexes");
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
r := PrintableOuterHull(cubes_two,coordinates_two,"Cubes/TwoCubes");
surf := r[1];
pr2 := SetVertexCoordinates3D(surf,r[3]);
DrawComplexToJavaScript(surf,"Cubes/TwoCubesXY_intersc_free",pr2);


ico := Icosahedron();

Coord3_2:= [
                    [  0.9510565160,  0.0000000000,  0.0000000000 ],
                    [  0.4253254040,  0.8506508085,  0.0000000000 ],
                    [  0.4253254040,  0.2628655560,  0.8090169940 ],
                    [ -0.0449027976, -0.0277514551,  0.0854101965 ],
                    [  0.4253254040, -0.6881909604, -0.4999999998 ],
                    [  0.4253254040, -0.6881909604,  0.4999999998 ],
                    [ -0.4253254040,  0.6881909604,  0.4999999998 ],
                    [ -0.4253254040,  0.6881909604, -0.4999999998 ],
                    [ -0.4253254040, -0.2628655560, -0.8090169940 ],
                    [ -0.4253254040, -0.8506508085,  0.0000000000 ],
                    [  0.0449027976,  0.0277514551, -0.0854101965 ],
                    [ -0.9510565160,  0.0000000000,  0.0000000000 ],
                    ];
                
pr3 := SetVertexCoordinates3D(ico,Coord3_2);

DrawComplexToJavaScript(ico,"ico_3_2_raw",pr3);
data := PrintableOuterHull(ico,Coord3_2,"ico_3_2");

surf := data[1];
pts := data[3];

pr3 := SetVertexCoordinates3D(surf,pts);
DrawComplexToJavaScript(surf,"ico_3_2",pr3);


# figure
DeactivateEdges(surf,pr3);
for e in RamifiedEdges(surf) do
	SetEdgeColour(surf,e,"0xFF0000",pr3);
	ActivateEdge(surf,e,pr3);
od;
#ActivateLineWidth(surf,pr3);


pr3.edgeThickness := 0.15;
DrawComplexToJavaScript(surf,"ico_3_2_red",pr3);

outer_triangle:=FindOuterTriangle(data[1],data[3]);

data_outer:=ExtractChamber(data[1],data[3],outer_triangle[1],outer_triangle[2]);

shift_param:=0.1;

f := RemedyNonManifold(data_outer,data[3],shift_param);

unram_surf := f[1];
pts := f[2];
pr4 := SetVertexCoordinates3D(unram_surf,pts);

DeactivateEdges(unram_surf,pr4);
old_ram_edges := [25,31,60,77,78,90,49,50,29,24,
		  86,84,87,61,62,75,72,74,45,46];
for e in old_ram_edges do
	SetEdgeColour(unram_surf,e,"0xFF0000",pr4);
	ActivateEdge(unram_surf,e,pr4);
od;
DeactivateEdge(unram_surf,76,pr4);
DrawComplexToJavaScript(unram_surf,"ico_3_2_unram_red",pr4);


coords := [[0.,0.,0.],[0.,0.5,0.],[-0.5,0.3,0.]]

