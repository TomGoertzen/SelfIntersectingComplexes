# Merges multiple sets of surface data into a single set
BindGlobal("MergeSurfacesCoordinates", function(data)
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
BingGlobal( "SimplicialSurfaceFromCoordinates", function(params,eps)
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
                pos := NumericalPosition(VerticesCoords,Coords[f][l],eps);
            
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
    end
);

InstallGlobalFunction(PrintableSymmetricOuterHull, function(t, points, name, eps, shift_param, group, group_ort, has_intersections, has_ram_edges)
    local data, f, n, order_data, unramified_data, unram_surf, unram_points;
    if has_intersections then
        data := calculate_intersections_groups(t, points, false, group, group_ort);
        points := data[1];
        t := TriangularComplexByVerticesInFaces(data[2]);
    fi;
    data := FindOuterTriangle(t, Sublist(points, Vertices(t)));
    f := data[1];
    n := data[2];
    data := OuterHull(t, points, f, n);
    t := ShallowCopy(data[2]);
    if has_ram_edges and Size(RamifiedEdges(t)) > 0 then
        order_data := OrderRamifiedEdges(t, ShallowCopy(data));
        unramified_data := FixRamPath(t, order_data, data, ShallowCopy(points), shift_param);
        t := unramified_data[1];
        points := unramified_data[2];
        DrawSTLwithNormals(t, name, points, -data[4], []);
        return [data[1], t, points, data[4]];
    fi;
    DrawSTLwithNormals(t, name, points, -data[4], []);
    return [data[1], t, points, data[4]];
end);

InstallGlobalFunction(PrintableOuterHull, function(t, points, name, eps, shift_param, group, has_intersections, has_ram_edges)
    local data, f, n, order_data, unramified_data, unram_surf, unram_points;
    if has_intersections then
        data := calculate_intersections(Set(VerticesOfFaces(t)), points, false, group);
        points := data[1];
        t := TriangularComplexByVerticesInFaces(data[2]);
    fi;
    data := FindOuterTriangle(t, Sublist(points, Vertices(t)));
    f := data[1];
    n := data[2];
    data := OuterHull(t, points, f, n);
    t := ShallowCopy(data[2]);
    if has_ram_edges and Size(RamifiedEdges(t)) > 0 then
        order_data := OrderRamifiedEdges(t, ShallowCopy(data));
        unramified_data := FixRamPath(t, order_data, data, ShallowCopy(points), shift_param);
        t := unramified_data[1];
        points := unramified_data[2];
        DrawSTLwithNormals(t, name, points, -data[4], []);
        return [data[1], t, points, data[4]];
    fi;
    DrawSTLwithNormals(t, name, points, -data[4], []);
    return [data[1], t, points, data[4]];
end);


InstallGlobalFunction(JoinComponents, function(t, coordinates, eps)
    local map, i, new_faces, new_vertices, map2, data;
    data := [coordinates, VerticesOfFaces(t)];
    data[1] := 1. * data[1]; # Ensure coordinates are in floating-point format
    # Delete multiple occurrences of vertices
    map := [];
    for i in [1..Size(data[1])] do
        map[i] := MyNumericalPosition(data[1], data[1][i], eps);
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



InstallGlobalFunction(MergeOuterHull, function(info, name, eps, has_self)
    local data, f, n, order_data, unramified_data, unram_surf, unram_points, t, points;
    # Parallelizable first step
    t := info[1];
    points := info[2];
    info := JoinComponents(t, points, eps);
    points := info[1];
    t := info[2];
    if has_self then
        data := calculate_intersections(Set(VerticesOfFaces(t)), points, false, Group(()));
        points := data[1];
        t := TriangularComplexByVerticesInFaces(data[2]);
    fi;
    data := FindOuterTriangle(t, Sublist(points, Vertices(t)));
    f := data[1];
    n := data[2];
    data := OuterHull(t, points, f, n);
    t := ShallowCopy(data[2]);
    DrawSTLwithNormals(t, name, points, -data[4], []);
    return [data[1], t, points, data[4]];
end);

InstallGlobalFunction(ComponentsOuterHull, function(info, name, eps, has_self)
    local data, f, n, order_data, unramified_data, unram_surf, unram_points, t, points, normals, new_vof, c,assembly;
    # Parallelizable first step
    t := info[1];
    points := info[2];
    normals := [];
    if has_self then
        assembly := calculate_intersections_comp(t, points, false);
        assembly := List(assembly, i -> [TriangularComplexByVerticesInFaces(i[2]), i[1]]);
        data := MergeSurfacesCoordinates(assembly);
        points := data[2];
        t := data[1];
    fi;
    new_vof := [];
    for c in ConnectedComponents(t) do
        data := FindOuterTriangle(c, Sublist(points, Vertices(c)));
        f := data[1];
        n := data[2];
        data := OuterHull(c, points, f, n);
        for f in Faces(ShallowCopy(data[2])) do
            new_vof[f] := VerticesOfFaces(ShallowCopy(data[2]))[f];
            normals[f] := data[4][f];
        od;
    od;
    t := TriangularComplexByVerticesInFaces(new_vof);
    DrawSTLwithNormals(t, name, points, -normals, []);
    return [info[1], t, points, -normals];
end);

InstallGlobalFunction(ComponentsOuterHullCombinatorics, function(info, eps, has_self)
    local data, f, n, order_data, unramified_data, unram_surf, unram_points, t, points, normals, new_vof, c, assembly, FixCoordinateFormat, InvertList;
    # Parallelizable first step
    t := info[1];
    points := info[2];
    normals := [];
    if has_self then
        assembly := calculate_intersections_comp(t, points, false);
        assembly := List(assembly, i -> [TriangularComplexByVerticesInFaces(i[2]), i[1]]);
        data := MergeSurfacesCoordinates(assembly);
        points := data[2];
        t := data[1];
    fi;
    new_vof := [];
    for c in ConnectedComponents(t) do
        data := FindOuterTriangle(c, Sublist(points, Vertices(c)));
        f := data[1];
        n := data[2];
        data := OuterHull(c, points, f, n);
        for f in Faces(ShallowCopy(data[2])) do
            new_vof[f] := VerticesOfFaces(ShallowCopy(data[2]))[f];
            normals[f] := data[4][f];
            # Reordering if necessary
            if Dot(Crossproduct(points[new_vof[f][2]] - points[new_vof[f][1]], points[new_vof[f][3]] - points[new_vof[f][1]]), normals[f]) < 0. then
                new_vof[f] := [new_vof[f][2], new_vof[f][1], new_vof[f][3]];
            fi;
        od;
    od;

    # FixCoordinateFormat and InvertList functions are used here
    
    return FixCoordinateFormat(new_vof, points);
end);



