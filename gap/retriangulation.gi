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



InstallGlobalFunction(Retriangulation,function(t,coordinates)
    local l,i,data_triangulated;
    l:=ComputeSelfIntersections(t,coordinates);
    data_triangulated:=[];
    for i in Faces(t) do
            data_triangulated[i]:=triangulate(DiscTriangulation(l[i]));
    od;
    return join_triangles(data_triangulated);
end);

InstallGlobalFunction(ComponentsRetriangulation,function(t,coordinates)
    local l,components,data_triangulated,c,i;
    l:=ComputeSelfIntersections(t,coordinates);
    components:=[];
    data_triangulated:=[];
    for c in ConnectedComponents(t) do
        for i in Faces(c) do
            data_triangulated[i]:=triangulate(DiscTriangulation(l[i]));
        od;
        Add(components,join_triangles(data_triangulated{Faces(c)}));
    od;
    return components;
end);

# Gets args[1] triangular complex
# args[2] coordinates
# args[3] combinatorial group on vertices
# args[4] symmetry group respecting coordinates as subgroup of O(3)
InstallGlobalFunction(SymmetricRetriangulation,function(args...)
    local l,k,i,j,intersection,intersection2,h,orbs_vof,data_triangulated,orb_rep,res,cur_line,inside_points,m,l1,l2,n,vof,coordinates,group,group_orthogonal,orbs,shift,new_data,g,eps,t,pairs_reps,orbs_pairs,p,orb;
    if Size(args)=3 then
        t:=args[1];
        vof:=VerticesOfFaces(args[1]);
        coordinates:=args[2];
        group:=args[3];
        orbs_vof:=Orbits(group,vof,OnSets);
        l:=[];
        for i in [1..Size(vof)] do
            l[i]:=[List([1..3],j->coordinates[vof[i][j]]),Set([Set([1,2]),Set([2,3]),Set([3,1])])];
        od;
        #for each representative find self intersections
        for k in [1..Size(orbs_vof)] do
            orb_rep:=orbs_vof[k][1];
            i:=Position(vof,orb_rep);
            for j in [1..Size(vof)] do
                # output of the form [vertices,intersection_edges]
                if i<>j then
                    # call by reference: add intersection points
                    # return values add edges
                    l[i][2]:=Concatenation(l[i][2],TwoTriangleIntersection(l[i][1],l[j][1]));
                fi;
            od;
        od;
        # add stabilizer
        for k in [1..Size(orbs)] do
            i:=orbs[k][1];
            new_data:=[StructuralCopy(l[i][1]),StructuralCopy(l[i][2])];
            shift:=Size(l[1]);
            j:=1;
            for g in Stabilizer(group,i) do
                if Size(Stabilizer(group,i))=1 then
                    break;
                fi;
                if g<>() then
                    l[i][1]:=Concatenation(l[i][1],new_data[1]*group_orthogonal[Position(Elements(group),g)]);
                    l[i][2]:=Concatenation(l[i][2],new_data[2]+j*Size(new_data[1]));
                    j:=j+1;
                fi;
            od;
        od;
        data_triangulated:=[];
        # transfer triangulations to all other triangles
        for k in [1..Size(orbs_vof)] do
            i:=Position(vof,orbs_vof[k][1]);
            data_triangulated[i]:=triangulate(DiscTriangulation(l[i]));
            for h in [2..Size(orbs_vof[k])] do
                j:=Position(vof,orbs_vof[k][h]);
                data_triangulated[j]:=SymmetryAction(data_triangulated[i],vof[i],vof[j],coordinates,group);
            od;
        od;
        return join_triangles(data_triangulated);
    elif Size(args)=4 then
        t:=args[1];
        vof:=VerticesOfFaces(args[1]);
        coordinates:=args[2];
        group:=args[3];
        group_orthogonal:=args[4];
        eps:=_SelfIntersectingComplexesParameters.eps;
        vof:=VerticesOfFaces(t);
        orbs:=Orbits(group,Faces(t));
        # calculate representatives for faces that intersect with face representatives
        orbs_pairs:=Orbits(group,Combinations(Faces(t),2),OnSets);
        pairs_reps:=List(Faces(t),f->[]);
        for k in [1..Size(orbs)] do
            i:=orbs[k][1];
            for orb in orbs_pairs do
                for p in Orbits(Stabilizer(group,i),orb,OnSets) do
                    if i in p[1] then
                        Add(pairs_reps[i],Other(p[1],i));
                    fi;
                od;
            od; 
        od;
        l:=[];
        for i in [1..Size(vof)] do
            l[i]:=[List([1..3],j->coordinates[vof[i][j]]),Set([Set([1,2]),Set([2,3]),Set([3,1])])];
        od;
        #for each representative find self intersections
        for k in [1..Size(orbs)] do
            i:=orbs[k][1];
            for j in pairs_reps[i] do
                # output of the form [vertices,intersection_edges]
                if i<>j then
                    l[i][2]:=Concatenation(l[i][2],TwoTriangleIntersection(l[i][1],l[j][1]));
                fi;
            od;
        od;
        # add stabilizer
        for k in [1..Size(orbs)] do
            i:=orbs[k][1];
            new_data:=[StructuralCopy(l[i][1]),StructuralCopy(l[i][2])];
            shift:=Size(l[1]);
            j:=1;
            for g in Stabilizer(group,i) do
                if Size(Stabilizer(group,i))=1 then
                    break;
                fi;
                if g<>() then
                    l[i][1]:=Concatenation(l[i][1],new_data[1]*group_orthogonal[Position(Elements(group),g)]);
                    l[i][2]:=Concatenation(l[i][2],new_data[2]+j*Size(new_data[1]));
                    j:=j+1;
                fi;
            od;
        od;
        data_triangulated:=[];
        # transfer triangulations to all other triangles
        for k in [1..Size(orbs)] do
            i:=orbs[k][1];
            #Print("Triangulation for face ",String(i),"\n");
            data_triangulated[i]:=triangulate(DiscTriangulation(l[i]));
            while not IsSimplicialSurface(TriangularComplexByVerticesInFaces(data_triangulated[i][2])) do
                Print("Triangulation of face ",String(i)," failed! Try to triangulate it again.\n");
                data_triangulated[i]:=triangulate(DiscTriangulation(l[i]));
            od;
            for h in [2..Size(orbs[k])] do
                j:=orbs[k][h];
                data_triangulated[j]:=SymmetryActionDirect(data_triangulated[i],i,j,coordinates,group,group_orthogonal);
            od;
        od;
        return join_triangles(data_triangulated);
    else
        Print("Unknown Number of arguments. Should be either 3 or 4.");
        return fail;
    fi;
end);