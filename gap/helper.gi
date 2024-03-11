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
    if x-y <= epsilon then
        return true;
    fi;
    return false;
end);

# Floating point greater than or equal check
BindGlobal("FlGeq",function(x,y,epsilon)
    if y-x <= epsilon then
        return true;
    fi;
    return false;
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


