# Merges multiple sets of surface data into a single set
BindGlobal("_MergeSurfacesCoordinates", function(data)
    local vof, coordinates, i;
    vof := []; # Initialize list for vertices of faces
    coordinates := []; # Initialize list for coordinates
    for i in [1..Size(data)] do
        # Concatenate vertices of faces, adjusting indices by the current size of coordinates
        vof := Concatenation(vof, VerticesOfFaces(data[i][1]) + Size(coordinates));
        # Concatenate coordinates
        coordinates := Concatenation(coordinates, data[i][2]);
    od;
    # Return the new structure comprising a triangular complex (assumed from context) and combined coordinates
    return [TriangularComplexByVerticesInFaces(vof), coordinates];
end);


# Create Simplicial Surface from Coordinate List
# Input is of type [IsList,IsFloat]
BindGlobal( "_SimplicialSurfaceFromCoordinates", function(params,eps)
        local Coords, faces, f, i, j, l, pos, VerticesInFaces, VerticesCoords, verts, surf;
        Coords := params[1];
        faces := params[2];
        
        # vertex counter
        i := 1;

        # face counter
        j := 1;
        
        VerticesInFaces := [];
        VerticesCoords := [];
        
        for f in faces do
            verts := [];

            for l in [1,2,3] do 
                pos := _NumericalPosition(VerticesCoords,Coords[f][l],eps);
            
                if pos = fail then
                    # vertex coord is new
                    VerticesCoords[i] := Coords[f][l];
                    
                    verts[l] := i;

                    i := i + 1;
                
                else
                    verts[l] := pos;
                fi;

                
            od;

            VerticesInFaces[j] := verts;
            j := j + 1;
        od;
        
        surf := TriangularComplexByVerticesInFaces(VerticesInFaces);

        return [surf,VerticesCoords];
end);

InstallGlobalFunction(PrintableSymmetricOuterHull, function(t, points, name,group,group_ort)
    local data, f, n, order_data, unramified_data, unram_surf, unram_points,coords;
    data := SymmetricRetriangulation(t, points, group, group_ort);
    points := data[1];
    t := TriangularComplexByVerticesInFaces(data[2]);
    data := FindOuterTriangle(t, _Sublist(points, Vertices(t)));
    f := data[1];
    n := data[2];
    data := ExtractChamber(t, points, f, n);
    t := ShallowCopy(data[2]);
    
    
    f := RemedyNonManifold(data,points,_SelfIntersectingComplexesParameters.non_manifold_edge_shift);
    unram_surf := f[1];
    unram_points := f[2];
    coords := f[3];
    order_data := f[4];
    
    DrawSTLwithNormals(unram_surf, name, unram_points, -data[4], []);
    # TODO: normals wont be completely correct after we remedy the mfd
    return [data[1], unram_surf, unram_points, data[4]];
end);

InstallGlobalFunction(PrintableOuterHull, function(t, points, name)
    local data, f, n, coords, order_data, unramified_data, unram_surf, unram_points;
    data := Retriangulation(t, points);
    points := data[1];
    t := TriangularComplexByVerticesInFaces(data[2]);
    data := FindOuterTriangle(t, _Sublist(points, Vertices(t)));
    f := data[1];
    n := data[2];
    data := ExtractChamber(t, points, f, n);
    t := ShallowCopy(data[2]);
    
    
    f := RemedyNonManifold(data,points,_SelfIntersectingComplexesParameters.non_manifold_edge_shift);
    unram_surf := f[1];
    unram_points := f[2];
    coords := f[3];
    order_data := f[4];
    
    DrawSTLwithNormals(unram_surf, name, unram_points, -data[4], []);
    # TODO: normals wont be completely correct after we remedy the mfd
    return [data[1], unram_surf, unram_points, data[4]];
end);


InstallGlobalFunction(JoinComponents, function(t, coordinates)
    local map, i, new_faces, new_vertices, map2, data,eps;
    data := [coordinates, VerticesOfFaces(t)];
    data[1] := 1. * data[1]; # Ensure coordinates are in floating-point format
    # Delete multiple occurrences of vertices
    eps:=_SelfIntersectingComplexesParameters.eps;
    map := [];
    for i in [1..Size(data[1])] do
        map[i] := _MyNumericalPosition(data[1], data[1][i], eps);
    od;
    new_vertices := [];
    map2 := [];
    for i in [1..Size(data[1])] do
        if i in map then
            Add(new_vertices, data[1][i]);
            map2[i] := Size(new_vertices);
        else
            map2[i] := map2[map[i]];
        fi;
    od;
    new_faces := List(data[2], j -> [map2[j[1]], map2[j[2]], map2[j[3]]]);
    data[2] := new_faces;
    data[1] := new_vertices;
    # Delete multiple occurrences of faces
    data[2] := Set(List(data[2], i -> Set(i)));
    return [data[1], TriangularComplexByVerticesInFaces(data[2])];
end);



InstallGlobalFunction(MergeOuterHull, function(info, name)
    local data, f, n, order_data, unramified_data, unram_surf, unram_points, t, points,eps;
    eps:=_SelfIntersectingComplexesParameters.eps;
    # Parallelizable first step
    t := info[1];
    points := info[2];
    info := JoinComponents(t, points, eps);
    points := info[1];
    t := info[2];
    data := Retriangulation(t, points);
    points := data[1];
    t := TriangularComplexByVerticesInFaces(data[2]);
    data := FindOuterTriangle(t, _Sublist(points, Vertices(t)));
    f := data[1];
    n := data[2];
    data := ExtractChamber(t, points, f, n);
    t := ShallowCopy(data[2]);
    DrawSTLwithNormals(t, name, points, -data[4], []);
    return [data[1], t, points, data[4]];
end);

InstallGlobalFunction(ComponentsOuterHull, function(t,points, name)
    local data, f, n, order_data, unramified_data, unram_surf, unram_points, normals, new_vof, c,assembly;
    # Parallelizable first step
    normals := [];
    assembly := ComponentsRetriangulation(t,points);
    data := _MergeSurfacesCoordinates(assembly);
    points := data[2];
    t := data[1];
    new_vof := [];
    for c in ConnectedComponents(t) do
        data := FindOuterTriangle(c, _Sublist(points, Vertices(c)));
        f := data[1];
        n := data[2];
        data := ExtractChamber(c, points, f, n);
        for f in Faces(ShallowCopy(data[2])) do
            new_vof[f] := VerticesOfFaces(ShallowCopy(data[2]))[f];
            normals[f] := data[4][f];
        od;
    od;
    t := TriangularComplexByVerticesInFaces(new_vof);
    DrawSTLwithNormals(t, name, points, -normals, []);
    return [t, points, -normals];
end);


