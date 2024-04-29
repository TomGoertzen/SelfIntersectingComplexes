InstallGlobalFunction(DiscTriangulation, function(data)
    local map, i, edge, poss_edges, no_poss_edges, new_edges, ok, e, res, data_tri, t, boundary_vertices, count_new, j, count, v, eps;
    eps:=_SelfIntersectingComplexesParameters.eps;
    # fix intersections of intersection
    # Print("Fix planar intersections \n");
    data := RectifyDiscIntersections(data);
    poss_edges := [];
    count := 0;
    data_tri := triangulate_comb(data);
    t := TriangularComplexByVerticesInFaces(data_tri);
    # find number of boundary edges = boundary vertices
    boundary_vertices := [];
    for i in [1..Size(data[1])] do
        for j in [1..3] do
            if j = 1 then
                e := [1, 2];
            elif j = 2 then
                e := [2, 3];
            elif j = 3 then
                e := [3, 1];
            fi;
            if i in e then
                Add(boundary_vertices, i);
                break;
            fi;
            if LineSegmentIntersectionColinear([data[1][i], data[1][e[1]]], [data[1][e[1]], data[1][e[2]]], eps)[1] and LineSegmentIntersectionColinear([data[1][i], data[1][e[2]]], [data[1][e[1]], data[1][e[2]]], eps)[1] then
                Add(boundary_vertices, i);
                break;
            fi;
        od;
    od;
    boundary_vertices := Set(boundary_vertices);
    v := Size(data[1]);
    if v - Size(data[2]) + (2 * Size(data[2]) - Size(boundary_vertices)) / 3 = 1 then
        return data;
    fi;
    # first calculate all possible edges and sort them by their distance
    # next we will add these edges until we receive a triangulated disc
    for i in [1..Size(data[1])] do
        for j in [i + 1..Size(data[1])] do
            if not Set([i, j]) in data[2] and not (i in InnerVertices(t) and j in InnerVertices(t)) then
                Add(poss_edges, [Set([i, j]), MyNorm(data[1][i] - data[1][j])]);
            fi;
        od;
    od;
    SortBy(poss_edges, i -> -i[2]);
    count_new := 0;
    # as long as we do not have a simplicial disc add edges
    # Print("Worst case estimate of loops: ", Size(poss_edges) * Size(data[2]), "\n");
    while v - Size(data[2]) + (2 * Size(data[2]) - Size(boundary_vertices)) / 3 <> 1 do
        ok := true;
        edge := Remove(poss_edges);
        i := edge[1][1];
        j := edge[1][2];
        for e in [1..Size(data[2])] do
            if Set([i, j]) in data[2] or (i in InnerVertices(t) and j in InnerVertices(t)) then
                ok := false;
                break;
            fi;
            count := count + 1;
            res := LineSegmentIntersection([data[1][i], data[1][j]], [data[1][data[2][e][1]], data[1][data[2][e][2]]], eps);
            if (res[1] and (not MyNumericalPosition([data[1][i], data[1][j]], res[2], eps) <> fail or not MyNumericalPosition([data[1][data[2][e][1]], data[1][data[2][e][2]]], res[2], eps) <> fail)) then
                ok := false;
                break;
            fi;
            if LineSegmentIntersectionColinear([data[1][i], data[1][j]], [data[1][data[2][e][1]], data[1][data[2][e][2]]], eps)[1] then
                ok := false;
                break;
            fi;
        od;
        if ok then
            Add(data[2], Set([i, j]));
            count_new := count_new + 1;
            if count_new = 5 then
                t := TriangularComplexByVerticesInFaces(triangulate_comb(data));
                count_new := 0;
            fi;
        fi;
    od;
    # Print("In the end we got ", count, " loops \n");
    return CleanData(data, eps);
end);
