# Dot Product
BindGlobal("Dot",function(a,b)
        return a[1]*b[1]+a[2]*b[2]+a[3]*b[3];
end);

BindGlobal("MyNorm",function(a)
        return Sqrt(Dot(a,a));
end);

BindGlobal("SurfaceTriangle",function(x,y,z)
    	local a,b,c;
        # surface of triangle
        a := MyNorm(x-y);
        b := MyNorm(x-z);
        c := MyNorm(y-z);      
        if AbsoluteValue((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c))<eps then
        	# bc of numerical errors, expression can be negative instead of 0
        	return 0.;
        else
        	return 0.25 * CubeRoot( (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c) )^(3./2.);
        fi;
end);

# Rounding Vector to Nullify Components Smaller than Epsilon
BindGlobal("MyRoundVector", function(n, eps)
    local i;
    n := 1.0*n; # Ensure n is treated as a float
    for i in [1..3] do
        if Sqrt(n[i]^2) < eps then
            n[i] := 0.;
        fi;
    od;
    return n;
end);

# Angle Between Two Vectors
BindGlobal("VectorAngle", function(a, b)
    return Acos(Dot(a, b) / (MyNorm(a) * MyNorm(b)));
end);

# Angle Between Two Vectors in the Plane Defined by a Normal
BindGlobal("VectorAnglePlane", function(a, b, n)
    return Atan2(Determinant(1.0*[a, b, n]), Dot(a, b));
end);

# Rotate a Vector About an Axis by an Angle Alpha
BindGlobal("Rotate_Vector_Axis", function(v, u, alpha)
    local r_alpha_u;
    u := u / MyNorm(u); # Normalize the axis
    r_alpha_u := Cos(alpha) * IdentityMat(3) + Sin(alpha) * [[0, -u[3], u[2]], [u[3], 0, -u[1]], [-u[2], u[1], 0]] + (1 - Cos(alpha)) * (TransposedMat([u]) * [u]);
    return (v * TransposedMat(r_alpha_u));
end);

# Check if three points are collinear
BindGlobal("PointsInOneLine",function(v1, v2, v3, eps)
    local area;
    area := SurfaceTriangle(v1,v2,v3);
    return AbsoluteValue(area)<eps;
end);

# Cross product of two vectors
BindGlobal("Crossproduct",function(x,y)
    return [x[2]*y[3]-x[3]*y[2],x[3]*y[1]-x[1]*y[3],x[1]*y[2]-x[2]*y[1]];
end);

# Returns a list of unique vectors from a list, considering numerical tolerance
BindGlobal("NumericalUniqueListOfLists",function(list, eps)
    local n, I, i, unique;
    
    n := Length(list);
    I := [2..n];
    unique := [];
    unique[1] := list[1];
    
    for i in I do
        if ForAll(unique, x-> not MyNorm(x-list[i])<eps) then
            unique[Length(unique)+1] := list[i];
        fi;
    od;
    return unique;
end);

# Projects a line on a triangle
BindGlobal("ProjectLineOnTriangle",function(x,l,d,normal,points, eps)
    local lambda, p, orthog, coords, o, o_true, i;
    if not AbsoluteValue(d - x*normal)<eps then
        # the point is not on the triangle plane already
        lambda := (d - x*normal) / (l*normal);
        
    elif AbsoluteValue(d - x*normal)<eps and not AbsoluteValue(l*normal)<eps then
        # the edge is not orthogonal to the triangle and the point already on the triangle plane
        lambda := 0.;
    else
        # then the edge is orthogonal to the triangle and the point already on the triangle plane
        orthog := points[1];
        coords := points[2];
        o := [(coords[1]*orthog[1] - x*orthog[1]) / (l*orthog[1]), (coords[2]*orthog[2] - x*orthog[2]) / (l*orthog[2]), (coords[3]*orthog[3] - x*orthog[3]) / (l*orthog[3])];
        o_true := [];
        for i in [1..3] do
            if AbsoluteValue(o[i])>eps then
                o_true[i] := o[i];
            else
                o_true[i] := Maximum(o)+ 1.;
            fi;
        od;
        lambda := Minimum(o_true);
        if AbsoluteValue(lambda-1.)>eps then
            # if the face is fully inside the triangle then we always need to triangulate
            lambda := 1.;
        fi;
    fi;
    p := x+lambda*l;
    return [p,lambda];
end);

# Floating point equality check
BindGlobal("FlEq",function(x,y,epsilon)
    return AbsoluteValue(x-y)<epsilon;
end);

# Vector equality check with floating point precision
BindGlobal("FlVEq",function(x,y,epsilon)
    return MyNorm(x-y)<epsilon;
end);

# Floating point less than or equal check
BindGlobal("FlLeq",function(x,y,epsilon)
    return x-y <= epsilon;
end);

# Floating point greater than or equal check
BindGlobal("FlGeq",function(x,y,epsilon)
    return y-x <= epsilon;
end);

# Check if two lines lie on the same plane
BindGlobal("OnSamePlane", function(l1, l2, eps)
    local v1, v2, v3, c1, c2, cross_product, norm_cross_product;
    
    v1 := l1[1] - l1[2]; # Vector of first line
    v2 := l1[1] - l2[2]; # Vector connecting a point from the first line to a point on the second line
    v3 := l1[1] - l2[1]; # Vector connecting a point from the first line to another point on the second line
    
    c1 := Crossproduct(v1, v2); # Cross product of v1 and v2
    c2 := Crossproduct(v1, v3); # Cross product of v1 and v3
    
    cross_product := Crossproduct(c1, c2); # Cross product of c1 and c2
    norm_cross_product := MyNorm(cross_product); # Norm of the cross product
    
    return norm_cross_product < eps; # If the norm is less than epsilon, lines are on the same plane
end);

# Check if a point lies inside a given triangle
BindGlobal("MyPointInTriangle", function(a, b, c, p, eps)
    local n1, n2, n3;
    # Unit normal for triangles (A, B, P), (B, C, P), and (C, A, P)
    n1 := Crossproduct(p-a, b-a);
    if MyNorm(n1) < eps then
        if MyNorm(p-a) + MyNorm(p-b) <= MyNorm(a-b) + eps then
            return true;
        else
            return false;
        fi;
    else
        n1 := n1 / MyNorm(n1);
    fi;

    n2 := Crossproduct(p-b, c-b);
    if MyNorm(n2) < eps then
        if MyNorm(p-b) + MyNorm(p-c) <= MyNorm(b-c) + eps then
            return true;
        else
            return false;
        fi;
    else
        n2 := n2 / MyNorm(n2);
    fi;

    n3 := Crossproduct(p-c, a-c);
    if MyNorm(n3) < eps then
        if MyNorm(p-a) + MyNorm(p-c) <= MyNorm(a-c) + eps then
            return true;
        else
            return false;
        fi;
    else
        n3 := n3 / MyNorm(n3);
    fi;

    # Check if dot products of normals are close to 1
    if AbsoluteValue(Dot(n1, n2) - 1) <= eps and AbsoluteValue(Dot(n2, n3) - 1) <= eps then
        return true;
    fi;
    return false;
end);

# Calculate the orthogonal transformation that maps one triangle onto another
BindGlobal("OrthogonalTransformation", function(triangle_i, triangle_j)
    return triangle_i^-1 * triangle_j;
end);

# Transfer data from one triangulated object to another via group action
BindGlobal("SymmetryAction", function(data_triangulated_i, vof_i, vof_j, coordinates, group)
    local group_elem, triangle_j, triangle_i, ort_elem, first, second;
    group_elem := RepresentativeAction(group, vof_i, vof_j, OnSets);    
    first := vof_i{[1, 2, 3]};
    second := List(first, i -> i^group_elem);
    triangle_i := coordinates{first};
    triangle_j := coordinates{second};
    ort_elem := OrthogonalTransformation(triangle_i, triangle_j);
    return [List(data_triangulated_i[1], i -> i * ort_elem), data_triangulated_i[2]];
end);

# Direct symmetry action transfer using a predefined orthogonal group
BindGlobal("SymmetryActionDirect", function(data_triangulated_i, i, j, coordinates, group, group_orthogonal)
    local group_elem, triangle_j, triangle_i, ort_elem, first, second;
    group_elem := RepresentativeAction(group, i, j);    
    ort_elem := group_orthogonal[Position(Elements(group), group_elem)];
    return [List(data_triangulated_i[1], i -> i * ort_elem), data_triangulated_i[2]];
end);

# Remove a specified item from a list and return the remaining item
BindGlobal("Other", function(list, i)
    list := StructuralCopy(ShallowCopy(list)); # Create a modifiable copy of the list
    Remove(list, Position(list, i)); # Remove the specified item
    return list[1]; # Return the remaining item
end);

# Finds the position of an entry in a list using a numerical tolerance
BindGlobal("MyNumericalPosition", function(list, entry, eps)
    local i;
    for i in [1..Size(list)] do
        if MyNorm(list[i] - entry) < eps then
            return i;
        fi;
    od;
    return fail; # Return 'fail' if the entry is not found within the numerical tolerance
end);

# Creates a sublist from the given list at the specified positions
BindGlobal("Sublist", function(list, positions)
    local new_list, p;
    new_list := []; # Initialize the new list
    for p in positions do
        new_list[p] := list[p]; # Assign elements from the original list to the new list based on positions
    od;
    return new_list; # Return the newly formed sublist
end);



BindGlobal("NumericalPosition", function(list,entry,epsilon)
        local i, n;

        n := Length(list);
        i := 1;
        
        while i <= n do
            if FlVEq(list[i],entry, epsilon) then
                return i;
            fi;
            i := i + 1;
        od;
        
        return fail;
    end
);

# creates a triangular complex from coordinate data
BindGlobal( "TriangularComplexFromCoordinates", function(params,eps)
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
                    # vertex coord. is new
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
);;


# Create a simplicial surface from changed coordinates and a previous simplicial surface (that determines the faces of the new simplicial surface
# Input is of type[IsList,IsFloat], where List = [Coordinates,SimplicialSurface]
BindGlobal( "SimplicialSurfaceFromChangedCoordinates", function(params,eps)
        local Coords, faces, f, i, j, l, pos, old_surf, VerticesInFaces, VerticesCoords, verts, surf;
        Coords := params[1];
        old_surf := params[2];
        faces := ShallowCopy(Faces(old_surf));
        
        VerticesInFaces := [];
        
        
        for f in faces do
            verts := Coords[f][5];

            VerticesInFaces[f] := verts;

        od;
        
        surf := TriangularComplexByVerticesInFaces(VerticesInFaces);

        return [surf];
    end
);



BindGlobal("UmbrellaPathFaceVertex",function(t,f,v,data,points)
	local faces,new_face,edges,edge,old_face;
	faces:=[];
	new_face:=f;
	edge:=Intersection(EdgesOfVertex(t,v),EdgesOfFace(t,new_face))[2];
	while not new_face in faces do
		Add(faces,new_face);
		# this way you only find one butterfly part of the non-manifold edge
		old_face:=new_face;
		if IsRamifiedEdge(t,edge) then
			new_face:=UpwardContinuation(t,edge,points,new_face,data[4][new_face])[1];
		else
			new_face:=Difference(FacesOfEdge(t,edge),[new_face])[1];
		fi;
		edge:=Difference(Intersection(EdgesOfVertex(t,v),EdgesOfFace(t,new_face)),EdgesOfFace(t,old_face))[1];
	od;
	return faces;
end);;

# outputs local umbrellas at vertex v
BindGlobal("UmbrellaComponents",function(t,v)
	local u_comps, available_edges, v_faces, v_edges, edges_of_faces, f, e, fa, edges, avaiable_edges, comp, connection, cur_edges;;
	
	u_comps := [];
	v_faces := ShallowCopy(FacesOfVertex(t,v));
	v_edges := ShallowCopy(EdgesOfVertex(t,v));
	edges_of_faces := [];
	for f in v_faces do
		 edges := ShallowCopy(EdgesOfFace(t,f));
		 edges_of_faces[f] := Intersection(edges,v_edges);
	od;
	while v_faces <> [] do
		fa := v_faces[1];
		Remove(v_faces,Position(v_faces,fa));
		comp := [fa];
		cur_edges := [edges_of_faces[fa][1]];
		
		available_edges := Unique(Flat(ShallowCopy(edges_of_faces)));
		# so that we dont land in an infinite loop
		Remove(available_edges,Position(available_edges,edges_of_faces[fa][2]));
		
		while Intersection(cur_edges,available_edges) <> [] do
			# while we can still reach a face on our side of the umbrella, continue
			for f in v_faces do
				connection := Intersection(edges_of_faces[f],cur_edges);
				if connection <> [] then
					Add(comp,f);
					cur_edges := Union(cur_edges,edges_of_faces[f]);
					Remove(v_faces,Position(v_faces,f));
					for e in connection do
						if e in available_edges then
							Remove(available_edges,Position(available_edges,e));
						fi;
					od;
				fi;
			od;
		od;
		Add(u_comps,comp);
	od;
	return u_comps;
end);;

# reads a stl file and creates a complex / simplicial surface (depending on the regularity of the object)
BindGlobal("ReadSTL", function(fileName)
	# reads a file from the current dir
	local surf, file, name, r, r2, eps, filesepr, faces, endsign, normal, data, normals, points, test, i,j, index, verts, coords, input, Coords;
	eps := 1./10^6;
	filesepr := SplitString(fileName, ".");
        name := filesepr[1];
        file := Filename( DirectoryCurrent(), Concatenation(name,".stl") );
        points := [];
        Coords:=[];
        i := 1;
       	# test file name
	if IsReadableFile(file) then
		
        	input := InputTextFile(file);
		r := ReadLine(input);
		Print(r);
		endsign := SplitString(ShallowCopy(r), " ")[1];
		
		while not endsign = "endsolid" do
			
			
			r := ReadLine(input);
			r2 := SplitString(ShallowCopy(r), " ");
			endsign := r2[1];
			if not endsign = "endsolid" then
				Coords[i]:=[];
				# TODO:  maybe find way to round less?
				
				normal := [Float(r2[3]),Float(r2[4]),Float(r2[5])];
				Coords[i][4]:=normal;
				
				r := ReadLine(input);
				
				j := 1;
				verts := [];
				while j < 4 do
					r := ReadLine(input);
					r2 := SplitString(ShallowCopy(r), " ");
					coords := [Float(r2[2]),Float(r2[3]),Float(r2[4])];
					
					test := ShallowCopy(points);
					Add(test,coords);
					if Length(NumericalUniqueListOfLists(test,eps)) > Length(points) then
						Add(points,coords);
						index := Length(points); 
					else
						index := NumericalPosition(points,coords,eps);
					fi;
					verts[j] := index;
					Coords[i][j] := coords;
					j := j+1;
				od;
				Coords[i][5] := verts;
				r := ReadLine(input);
				r := ReadLine(input);
				i := i + 1;
			fi;
			
		od;
		
	else
		Print("file does not exist");
		return 0;
	fi;
	
	faces := [1..Length(Coords)];
	data := TriangularComplexFromCoordinates([Coords,faces],eps);
	return [data[1],Coords,points];
end);;



BindGlobal("ConvertDataFormatPtC", function(surf,points,normals)
    local faces, f, Coords, verts, c_verts;

    faces := Faces(surf);
    Coords := [];
    for f in faces do
        Coords[f] := [];
        verts := ShallowCopy(VerticesOfFace(surf,f));
        c_verts := [points[verts[1]],points[verts[2]],points[verts[3]]];

        Coords[f][1] := c_verts[1];
        Coords[f][2] := c_verts[2];
        Coords[f][3] := c_verts[3];
        Coords[f][4] := normals[f];
        Coords[f][5] := verts;
    od;
    return Coords;
end);;


BindGlobal("ConvertDataFormatCtP", function(surf,Coords)
    local f, faces, points, normals, c_p, normal, v;

    faces := Faces(surf);
    points := [];
    normals := [];

    for f in faces do
        c_p := [Coords[f][1],Coords[f][2],Coords[f][3]];
        normal := Coords[f][4];
        v := Coords[f][5];
        
        points[v[1]] := c_p[1];
        points[v[2]] := c_p[2];
        points[v[3]] := c_p[3];
        normals[f] := normal;
    od;
    return [points,normals];
end);;
