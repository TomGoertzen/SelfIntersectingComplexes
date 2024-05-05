# Triangulate the data by finding all triangles formed by the vertices and edges
BindGlobal("_triangulate", function(data)
    local e, faces, is_triangle, i, j;
    faces := [];
    data[2] := Set(List(data[2], i -> Set(i))); # Ensure edges are unique
    for e in data[2] do
        for i in [1..Size(data[1])] do
            if Set([i, e[1]]) in data[2] and Set([i, e[2]]) in data[2] then
                is_triangle := true;
                for j in [1..Size(data[1])] do
                    if j <> i and j <> e[1] and j <> e[2] and _MyPointInTriangle(data[1][e[1]], data[1][e[2]], data[1][i], data[1][j], _SelfIntersectingComplexesParameters.eps) then
                        is_triangle := false;
                    fi;
                od;
                if is_triangle then
                    Add(faces, [e[1], e[2], i]);
                fi;
            fi;
        od;
    od;
    return [data[1], Set(List(faces, i -> Set(i)))];
end);

# Join the data of all triangles into a single structure
BindGlobal("_join_triangles", function(list)
    local vertices, VerticesOfFaces, map, i, pos, data;
    vertices := [];
    VerticesOfFaces := [];
    for data in list do
        map := [];
        for i in [1..Size(data[1])] do
            pos := _MyNumericalPosition(vertices, data[1][i], _SelfIntersectingComplexesParameters.eps);
            if pos = fail then
                Add(vertices, data[1][i]);
                map[i] := Size(vertices);
            else
                map[i] := pos;
            fi;
        od;
        VerticesOfFaces := Concatenation(VerticesOfFaces, List(data[2], triangle -> Set([map[triangle[1]], map[triangle[2]], map[triangle[3]]])));
    od;
    return [vertices, Set(VerticesOfFaces)];
end);

# Computes intersection of two Line Segments in 3D
# https://math.stackexchange.com/questions/270767/find-intersection-of-two-3d-lines
BindGlobal("_LineSegmentIntersection",function(l1,l2,eps)
	local C,D,e,f,g,sol,dotsol1,dotsol2,_MyNorm_cross_e_f,cross_e_f,cross_f_g;
	e:=l1[2]-l1[1];
	f:=l2[2]-l2[1];
	cross_e_f:=_Crossproduct(e,f);
	_MyNorm_cross_e_f:=_MyNorm(cross_e_f);
	if _MyNorm_cross_e_f > eps then
		if _OnSamePlane(l1,l2,eps) then
			C:=l1[1];
			D:=l2[1];
			g:=D-C;
			cross_f_g:=_Crossproduct(f,g);
			if _Dot(cross_f_g,-cross_e_f) >= 0. then
				sol := C + _MyNorm(cross_f_g)/_MyNorm_cross_e_f*e;
			else
				sol := C - _MyNorm(cross_f_g)/_MyNorm_cross_e_f*e;
			fi;
			dotsol1 := _Dot(l1[1]-l1[2],sol-l1[2]);
			dotsol2 := _Dot(l2[1]-l2[2],sol-l2[2]);
			if _MyNorm(l1[1]-sol)+_MyNorm(l1[2]-sol)<=_MyNorm(l1[1]-l1[2])+eps and _MyNorm(l2[1]-sol)+_MyNorm(l2[2]-sol)<=_MyNorm(l2[1]-l2[2])+eps then
				if not _MyNumericalPosition(l1,sol,eps)<>fail or not _MyNumericalPosition(l2,sol,eps)<> fail then
					return [true,sol];
				fi;
			fi;
		fi;
	fi;
	return [false];
end);


BindGlobal("_LineSegmentIntersectionColinear",function(l1,l2,eps)
	local v1,v2,x1,x2,y1,y2,r1,r2,s1,s2,i,j;
	l1:=1.*l1;
	l2:=1.*l2;
	v1:=l1[2]-l1[1];
	v2:=l2[2]-l2[1];
	if not _OnSamePlane(l1,l2,eps) then
		return [false];
	fi;
	if not _MyNorm(_Crossproduct(v1,v2)) < eps then
		return [false];
	fi;
	if not _PointsInOneLine(l1[1],l1[2],l2[1],eps) then
		return [false];
	fi;
	if not _PointsInOneLine(l1[1],l1[2],l2[2],eps) then
		return [false];
	fi;
			x1:=l1[1];
			x2:=l2[1];
			y1:=l1[2];
			y2:=l2[2];
			r1:=[];
			r2:=[];
			s1:=[];
			s2:=[];
			for i in [1..3] do
				if Sqrt(v1[i]^2)>eps then
					r1[i]:=(x2[i]-x1[i])/v1[i];
					r2[i]:=(y2[i]-x1[i])/v1[i];
				fi;
				if Sqrt(v2[i]^2)>eps then
					s1[i]:=(x1[i]-x2[i])/v2[i];
					s2[i]:=(y1[i]-x2[i])/v2[i];
				fi;
			od;
			# check if r1 has same non-zero entries
			for i in BoundPositions(r1) do
				for j in BoundPositions(r1) do
					if r1[i]^2>eps and r1[j]^2>eps and (r1[i]-r1[j])^2>eps then
						return [false];
					fi;
				od;
			od;
			for i in BoundPositions(r1) do
				if r1[i]^2>eps then
					r1:=r1[i];
					break;
				fi;
			od;
			if IsList(r1) then
				r1:=0.;
			fi;
			for i in BoundPositions(r2) do
				if r2[i]^2>eps then
					r2:=r2[i];
					break;
				fi;
			od;
			if IsList(r2) then
				r2:=0.;
			fi;
			for i in BoundPositions(s1) do
				if s1[i]^2>eps then
					s1:=s1[i];
					break;
				fi;
			od;
			if IsList(s1) then
				s1:=0.;
			fi;
			for i in BoundPositions(s2) do
				if s2[i]^2>eps then
					s2:=s2[i];
					break;
				fi;
			od;
			if IsList(s2) then
				s2:=0.;
			fi;
			if eps<= r1 and r1 <=1-eps then
				if eps<=s1 and s1<=1-eps then
					return [true,[[x1,x1+r1*v1],[x1+r1*v1,y2-s1*v2],[y2-s1*v2,x2]],1];
				elif eps<=s2 and s2 <=1-eps then
					return [true,[[x1,x1+r1*v1],[x1+r1*v1,x2+s2*v2],[x2+s2*v2,y2]],2];
				fi;
			fi;
			if eps<=r2 and r2 <=1-eps then
				if eps<=s1 and s1 <=1-eps then
					return [true,[[y1,y1-r2*v1],[y1-r2*v1,y2-s1*v2],[y2-s1*v2,x2]],3];
				elif eps<=s2 and s2<=1-eps then
					return [true,[[y1,y1-r2*v1],[y1-r2*v1,x2+s2*v2],[x2+s2*v2,y2]],4];
				fi;
			fi;
			# test if l1 contains l2 and a point of l2

			if _MyNorm(x1-x2)<eps then
				#r1=0 and s1=0
				if eps<=r2 and r2 <=1-eps then
					return [true,[[x2,y2],[y2,y1]],5];
				fi;
				if eps<=s2 and s2<=1-eps then
					return [true,[[x1,y1],[y1,y2]],6];
				fi;
			elif _MyNorm(y1-x2)<eps then
				#s2=0 and r1=1
				if eps<=r2 and r2 <=1-eps then
					return [true,[[x1,y2],[y2,x2]],7];
				fi;
				if eps<=s1 and s1<=1-eps then
					return [true,[[y1,x1],[x1,y2]],8];
				fi;
			elif _MyNorm(y2-x1)<eps then
				#r2=0 and s1=1
				if eps<=r1 and r1 <=1-eps then
					return [true,[[y2,x2],[x2,y1]],9];
				fi;
				if eps<=s2 and s2<=1-eps then
					return [true,[[x2,y1],[y1,x1]],10];
				fi;
			elif _MyNorm(y2-y1)<eps then
				#r2=1 and s2=1
				if eps<=r1 and r1 <=1-eps then
					return [true,[[x1,x2],[x2,y2]],11];
				fi;
				if eps<=s1 and s1<=1-eps then
					return [true,[[x2,x1],[x1,y1]],12];
				fi;
			fi;
	return [false];
end);

BindGlobal("_TwoTriangleIntersectionNonPlanar",function(coords_j,coords_k,eps)
        local l,ccoords,x,y,normal,t_coords,t_normal,ints,different,Compared,orthog,c_coords,numb_int,not_on_edge,ints_points,o,p,res,c_normal,d1,diff,not_on_edges,possb_intersec,lambda,edge_param,vOedge,d2;
        ccoords:=coords_j;
        x := ccoords[2]-ccoords[1];
        y := ccoords[3]-ccoords[1];
        normal := _Crossproduct(x,y);
        normal := normal / Sqrt(normal*normal);
        coords_j[4] := normal;
        ccoords:=coords_k;
        x := ccoords[2]-ccoords[1];
        y := ccoords[3]-ccoords[1];
        normal := _Crossproduct(x,y);
        normal := normal / Sqrt(normal*normal);
        coords_k[4] := normal;
        ints := false;
        different := true;            
        ints_points := [];
        Compared := [];
        possb_intersec := 0;
        # set coords again since the coordinate matrix may change every iteration
        different := true;
        if coords_j <> [] then
        c_normal := coords_j[4];
        d1 := coords_j[1]*c_normal;
        c_coords := coords_j{[1..3]};
        fi;
        t_coords := ShallowCopy(coords_k{[1..3]});
        # cant intersect if incident                       
    t_normal := coords_k[4];
    d2 := coords_k[1]*t_normal;
    orthog := [_Crossproduct(c_coords[2]-c_coords[1],c_normal),_Crossproduct(c_coords[3]-c_coords[2],c_normal),_Crossproduct(c_coords[1]-c_coords[3],c_normal)];
    
    orthog[1] := orthog[1] / Sqrt(orthog[1]*orthog[1]);
    orthog[2] := orthog[2] / Sqrt(orthog[2]*orthog[2]);
    orthog[3] := orthog[3] / Sqrt(orthog[3]*orthog[3]);
    
    # must be right of planes, need to have the orthogonal vectors point to the inside of the triangle
    orthog[1] := orthog[1] * ((orthog[1]*c_coords[3]-orthog[1]*c_coords[1]) / AbsoluteValue(orthog[1]*c_coords[3]-orthog[1]*c_coords[1]));
    orthog[2] := orthog[2] * ((orthog[2]*c_coords[1]-orthog[2]*c_coords[2]) / AbsoluteValue(orthog[2]*c_coords[1]-orthog[2]*c_coords[2]));
    orthog[3] := orthog[3] * ((orthog[3]*c_coords[2]-orthog[3]*c_coords[3]) / AbsoluteValue(orthog[3]*c_coords[2]-orthog[3]*c_coords[3]));
    # check if triangle k intersects with j
    numb_int := 0;
    ints_points := [];
    # use this to compute if both of the intersection points lie on an edge
    not_on_edge := [true, true, true];
    for l in [1..3] do
        vOedge := [coords_k[l],coords_k[l mod 3 + 1]];
        # if two faces share a edge then this edge cant be an intersection
        if Length(_NumericalUniqueListOfLists([c_coords[1],c_coords[2],c_coords[3],vOedge[1],vOedge[2]],eps)) >  3 then
            edge_param := coords_k[l mod 3 + 1]-coords_k[l];                                  
            # find point that intersects with plane
            res := _ProjectLineOnTriangle(coords_k[l],edge_param, d1, c_normal,[orthog,c_coords],eps);                                   
            p := res[1];
            lambda := res[2];                                  
            o := [orthog[1]*p-orthog[1]*c_coords[1], orthog[2]*p-orthog[2]*c_coords[2], orthog[3]*p-orthog[3]*c_coords[3]];           
            # check if point inside triangle             
            if _FlGeq(Minimum(o[1],o[2],o[3]),0.,eps) and _FlGeq(lambda,0.,eps) and _FlLeq(lambda,1.,eps) then           
                not_on_edge[1] := not_on_edge[1] and _FlEq(0.,o[1],eps); 
                not_on_edge[2] := not_on_edge[2] and _FlEq(0.,o[2],eps); 
                not_on_edge[3] := not_on_edge[3] and _FlEq(0.,o[3],eps);           
                ints := true;
                possb_intersec := possb_intersec + 1;
                if EqFloat(p[1]+0.0,p[1]) and EqFloat(p[2]+0.0,p[2]) and EqFloat(p[3]+0.0,p[3]) then
                  numb_int := numb_int + 1;
                  ints_points[numb_int] := p; 
                fi; 
            fi;                                   
        fi;                            
    od;
    not_on_edges := not( not_on_edge[1] or not_on_edge[2] or not_on_edge[3]);                        
    if Length(ints_points) > 0 then
        ints_points := _NumericalUniqueListOfLists(ints_points,eps);
        diff := numb_int - Length(ints_points);
        possb_intersec := possb_intersec - diff;
    fi;                        
    numb_int := Length(ints_points);
    return ints_points;
end);

# output of the form [vertices,intersection_edges]
InstallGlobalFunction(TwoTriangleIntersection,function(triangle1,triangle2)
	local edges,n,m,inside_points,cur_line,l1,l2,res,intersection,intersection2,eps;
	eps:=_SelfIntersectingComplexesParameters.eps;
	edges:=[];
	#check if triangles are coplanar and fix cases
	n:=_Crossproduct(triangle1[2]-triangle1[1],triangle1[3]-triangle1[1]);
	if (_Dot(n,triangle2[1])-_Dot(n,triangle1[1]))^2<eps and (_Dot(n,triangle2[2])-_Dot(n,triangle1[1]))^2<eps and (_Dot(n,triangle2[3])-_Dot(n,triangle1[1]))^2<eps then
		inside_points:=[];
		for m in [1..3] do
			if _MyPointInTriangle(triangle1[1],triangle1[2],triangle1[3],triangle2[m],eps) then
				Add(inside_points,m);
			fi;
		od;
		if Size(inside_points)=2 then
			Add(triangle1,triangle2[inside_points[1]]);
			Add(triangle1,triangle2[inside_points[2]]);
			AddSet(edges,[Size(triangle1),Size(triangle1)-1]);
		elif Size(inside_points)=3 then
			Add(triangle1,triangle2[inside_points[1]]);
			Add(triangle1,triangle2[inside_points[2]]);
			Add(triangle1,triangle2[inside_points[3]]);
			AddSet(edges,[Size(triangle1),Size(triangle1)-1]);
			AddSet(edges,[Size(triangle1)-2,Size(triangle1)-1]);
			AddSet(edges,[Size(triangle1),Size(triangle1)-2]);
		fi;
		#if Size(inside_points)=0 then
			# no inside points check for line intersections
			for l1 in Combinations([1..3],2) do
				cur_line:=[];
				for l2 in Combinations([1..3],2) do
					res:=_LineSegmentIntersection(triangle1{l2},triangle2{l1},eps);
					if res[1] then
						if Size(res[2])=3 then
							if _MyNumericalPosition(triangle1{[1,2,3]},res[2],eps)=fail then
								Add(cur_line,res[2]);
							fi;
						fi;
					fi;
				od;
				if Size(cur_line)=1 then
					#test if point of j from edge l1 lies inside i
						if l1[1] in inside_points and Size(_NumericalUniqueListOfLists(Concatenation([triangle2[l1[1]]],cur_line),eps))=2 then
							Add(triangle1,cur_line[1]);
							Add(triangle1,triangle2[l1[1]]);
							AddSet(edges,[Size(triangle1),Size(triangle1)-1]);
						elif l1[2] in inside_points and Size(_NumericalUniqueListOfLists(Concatenation([triangle2[l1[2]]],cur_line),eps))=2 then
							Add(triangle1,cur_line[1]);
							Add(triangle1,triangle2[l1[2]]);
							AddSet(edges,[Size(triangle1),Size(triangle1)-1]);
						else
							Add(triangle1,cur_line[1]);
						fi;
				elif Size(cur_line)=2 then
					Add(triangle1,cur_line[1]);
					Add(triangle1,cur_line[2]);
					AddSet(edges,[Size(triangle1),Size(triangle1)-1]);
				fi;
			od;
	else
		intersection:=_TwoTriangleIntersectionNonPlanar(triangle1{[1..3]},triangle2{[1..3]},eps);
		#test if both intersection and intersection2 are the same
		#otherwise join for lines
		intersection2:=_TwoTriangleIntersectionNonPlanar(triangle2{[1..3]},triangle1{[1..3]},eps);
		if (Size(intersection)=1 or Size(intersection)=0) and Size(intersection2)>0 then
			if Size(_NumericalUniqueListOfLists(Concatenation(intersection,intersection2),eps))=2 then
				intersection:=_NumericalUniqueListOfLists(Concatenation(intersection,intersection2),eps);
			fi;
		fi;
		if Size(intersection)=2 then
			# add intersection informations only to triangle1
			# triangle2 will be obtained later
			Add(triangle1,intersection[1]);
			Add(triangle1,intersection[2]);
			AddSet(edges,[Size(triangle1),Size(triangle1)-1]);
		fi;
		if Size(intersection)=1 then
			Add(triangle1,intersection[1]);
		fi;
	fi;
	return edges;
end);


InstallGlobalFunction(ComputeSelfIntersections,function(t,coordinates)
	local vof,i,l,j;
	vof:=VerticesOfFaces(t);
	l:=[];
	for i in Faces(t) do
		l[i]:=[List([1..3],j->coordinates[vof[i][j]]),Set([Set([1,2]),Set([2,3]),Set([3,1])])];
	od;
	#for each representative find self intersections
	for i in Faces(t) do
		for j in Faces(t) do
			# output of the form [vertices,intersection_edges]
			if i<>j then
				# call by reference: add intersection points
				# return values add edges
				l[i][2]:=Concatenation(l[i][2],TwoTriangleIntersection(l[i][1],l[j][1]));
			fi;
		od;
	od;
	return l;
end);



# Clean the data by removing duplicate vertices and edges
BindGlobal("_CleanData", function(data, eps)
    local map, i, new_edges, new_vertices, map2;
    data[1] := 1. * data[1]; # Ensure vertices are in floating-point format
    # Delete multiple occurrences of vertices
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
    new_edges := List(data[2], j -> [map2[j[1]], map2[j[2]]]);
    data[2] := new_edges;
    data[1] := new_vertices;
    # Delete multiple occurrences of edges
    data[2] := Set(List(data[2], i -> Set(i)));
    return data;
end);




InstallGlobalFunction(RectifyDiscIntersections,function(data)
	local i,j,res,entry,entry1,entry2,entry_in_i,entry_in_j,addy,check_edges,orig_j,orig_i,cur,k,remove,e,l,n,eps;
	eps:=_SelfIntersectingComplexesParameters.eps;
	data:=_CleanData(data,eps);
	data[2]:=List(data[2]);
	for i in [1..Size(data[1])] do
		check_edges:=ShallowCopy(data[2]);
		while check_edges<>[] do
			e:=Remove(check_edges);
			if _LineSegmentIntersectionColinear([data[1][i],data[1][e[1]]],[data[1][e[1]],data[1][e[2]]],eps)[1] and _LineSegmentIntersectionColinear([data[1][i],data[1][e[2]]],[data[1][e[1]],data[1][e[2]]],eps)[1] and not i in e then
				orig_j:=ShallowCopy(e);
				if not Set([e[1],i]) in data[2] or not Set([e[2],i]) in data[2] then
					data[2]:=Concatenation(data[2],Set([Set([e[1],i]),Set([e[2],i])]));
					Remove(data[2],Position(data[2],orig_j));
					check_edges:=Concatenation(check_edges,Set([Set([e[1],i]),Set([e[2],i])]));
				fi;
			fi;
		od;
	od;
	# check for non-paralel intersecting lines
	data:=_CleanData(data,eps);
	data[2]:=List(data[2]);
	check_edges:=Combinations([1..Size(data[2])],2);;
	while check_edges<>[] do
		cur:=Remove(check_edges,1);
		i:=cur[1];
		j:=cur[2];
		if i <=Size(data[2]) and j<=Size(data[2]) and i in BoundPositions(data[2]) and j in BoundPositions(data[2]) then
			res:=_LineSegmentIntersection([data[1][data[2][i][1]],data[1][data[2][i][2]]],[data[1][data[2][j][1]],data[1][data[2][j][2]]],eps);
			if res[1] then
				entry:=_MyNumericalPosition(data[1],res[2],eps);
				if entry = fail then
					Add(data[1],res[2]);
					entry:=Size(data[1]);
					orig_i:=data[2][i];
					orig_j:=data[2][j];
					n:=Size(data[2]);
					data[2]:=Concatenation(data[2],[Set([data[2][i][1],entry]),Set([data[2][i][2],entry]),Set([data[2][j][1],entry]),Set([data[2][j][2],entry])]);
					Unbind(data[2][Position(data[2],orig_i)]);
					Unbind(data[2][Position(data[2],orig_j)]);
					for k in [1..4] do
						for l in [1..n] do
							Add(check_edges,[l,n+k]);
						od;
					od;
				fi;
			fi;
		fi;
		for k in BoundPositions(data[2]) do
			if Size(data[2][k])=1 then
				Error();
			fi;
		od;
	od;
	data:=_CleanData(data,eps);
	data[2]:=List(data[2]);
	for i in [1..Size(data[1])] do
		check_edges:=ShallowCopy(data[2]);
		while check_edges<>[] do
			e:=Remove(check_edges);
			if _LineSegmentIntersectionColinear([data[1][i],data[1][e[1]]],[data[1][e[1]],data[1][e[2]]],eps)[1] and _LineSegmentIntersectionColinear([data[1][i],data[1][e[2]]],[data[1][e[1]],data[1][e[2]]],eps)[1] and not i in e then
				orig_j:=ShallowCopy(e);
				Remove(data[2],Position(data[2],orig_j));
				if not Set([e[1],i]) in data[2] or not Set([e[2],i]) in data[2] then
					data[2]:=Concatenation(data[2],Set([Set([e[1],i]),Set([e[2],i])]));
					check_edges:=Concatenation(check_edges,Set([Set([e[1],i]),Set([e[2],i])]));
				fi;
			fi;
		od;

	od;
	return _CleanData(data,eps);;
end);


BindGlobal("_triangulate_comb",function(data)
	local e,faces,is_triangle,i,j;
	faces:=[];
	data[2]:=Set(List(data[2],i->Set(i)));
	for e in data[2] do
		for i in [1..Size(data[1])] do
			if Set([i,e[1]]) in data[2] and Set([i,e[2]]) in data[2] then
				Add(faces,[e[1],e[2],i]);
			fi;
		od;
	od;
	return Set(List(faces,i->Set(i)));
end);