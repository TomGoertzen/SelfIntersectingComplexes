# Triangulate the data by finding all triangles formed by the vertices and edges
BindGlobal("triangulate", function(data)
    local e, faces, is_triangle, i, j;
    faces := [];
    data[2] := Set(List(data[2], i -> Set(i))); # Ensure edges are unique
    for e in data[2] do
        for i in [1..Size(data[1])] do
            if Set([i, e[1]]) in data[2] and Set([i, e[2]]) in data[2] then
                is_triangle := true;
                for j in [1..Size(data[1])] do
                    if j <> i and j <> e[1] and j <> e[2] and MyPointInTriangle(data[1][e[1]], data[1][e[2]], data[1][i], data[1][j], eps) then
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
BindGlobal("join_triangles", function(list)
    local vertices, VerticesOfFaces, map, i, pos, data;
    vertices := [];
    VerticesOfFaces := [];
    for data in list do
        map := [];
        for i in [1..Size(data[1])] do
            pos := MyNumericalPosition(vertices, data[1][i], eps);
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

InstallGlobalFunction(TwoTriangleIntersection,function(coords_j,coords_k,eps)
        local l,ccoords,x,y,normal,t_coords,t_normal,ints,different,Compared,orthog,c_coords,numb_int,not_on_edge,ints_points,o,p,res,c_normal,d1,diff,not_on_edges,possb_intersec,lambda,edge_param,vOedge,d2;
        ccoords:=coords_j;
        x := ccoords[2]-ccoords[1];
        y := ccoords[3]-ccoords[1];
        normal := Crossproduct(x,y);
        normal := normal / Sqrt(normal*normal);
        coords_j[4] := normal;
        ccoords:=coords_k;
        x := ccoords[2]-ccoords[1];
        y := ccoords[3]-ccoords[1];
        normal := Crossproduct(x,y);
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

                        orthog := [Crossproduct(c_coords[2]-c_coords[1],c_normal),Crossproduct(c_coords[3]-c_coords[2],c_normal),Crossproduct(c_coords[1]-c_coords[3],c_normal)];
                        
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
                                # if two faces share a edge than this edge cant be an intersection
                                if Length(NumericalUniqueListOfLists([c_coords[1],c_coords[2],c_coords[3],vOedge[1],vOedge[2]],eps)) >  3 then
                                    edge_param := coords_k[l mod 3 + 1]-coords_k[l];                                  
                                    # find point that intersects with plane
                                    res := ProjectLineOnTriangle(coords_k[l],edge_param, d1, c_normal,[orthog,c_coords],eps);                                   
                                    p := res[1];
                                    lambda := res[2];                                  
                                    o := [orthog[1]*p-orthog[1]*c_coords[1], orthog[2]*p-orthog[2]*c_coords[2], orthog[3]*p-orthog[3]*c_coords[3]];
                                    
                                    # check if point inside triangle
                                     
                                    if FlGeq(Minimum(o[1],o[2],o[3]),0.,eps) and FlGeq(lambda,0.,eps) and FlLeq(lambda,1.,eps) then
                                    
                                        not_on_edge[1] := not_on_edge[1] and FlEq(0.,o[1],eps); 
                                        not_on_edge[2] := not_on_edge[2] and FlEq(0.,o[2],eps); 
                                        not_on_edge[3] := not_on_edge[3] and FlEq(0.,o[3],eps); 
                                        
                                        ints := true;
                                        possb_intersec := possb_intersec + 1;
                                        numb_int := numb_int + 1;
                                        ints_points[numb_int] := p; 	
                                    fi;                                   
                                fi;                            
                        od;
                        not_on_edges := not( not_on_edge[1] or not_on_edge[2] or not_on_edge[3]);                        
                        if Length(ints_points) > 0 then
                            ints_points := NumericalUniqueListOfLists(ints_points,eps);
                            diff := numb_int - Length(ints_points);
                            possb_intersec := possb_intersec - diff;
                        fi;                        
                        numb_int := Length(ints_points);
                        return ints_points;
end);

# Computes intersection of two Line Segments in 3D
# https://math.stackexchange.com/questions/270767/find-intersection-of-two-3d-lines
BindGlobal("LineSegmentIntersection",function(l1,l2,eps)
	local C,D,e,f,g,sol,dotsol1,dotsol2,MyNorm_cross_e_f,cross_e_f,cross_f_g;
	e:=l1[2]-l1[1];
	f:=l2[2]-l2[1];
	cross_e_f:=Crossproduct(e,f);
	MyNorm_cross_e_f:=MyNorm(cross_e_f);
	if MyNorm_cross_e_f > eps then
		if OnSamePlane(l1,l2,eps) then
			C:=l1[1];
			D:=l2[1];
			g:=D-C;
			cross_f_g:=Crossproduct(f,g);
			if Dot(cross_f_g,-cross_e_f) >= 0. then
				sol := C + MyNorm(cross_f_g)/MyNorm_cross_e_f*e;
			else
				sol := C - MyNorm(cross_f_g)/MyNorm_cross_e_f*e;
			fi;
			dotsol1 := Dot(l1[1]-l1[2],sol-l1[2]);
			dotsol2 := Dot(l2[1]-l2[2],sol-l2[2]);
			if MyNorm(l1[1]-sol)+MyNorm(l1[2]-sol)<=MyNorm(l1[1]-l1[2])+eps and MyNorm(l2[1]-sol)+MyNorm(l2[2]-sol)<=MyNorm(l2[1]-l2[2])+eps then
				if not MyNumericalPosition(l1,sol,eps)<>fail or not MyNumericalPosition(l2,sol,eps)<> fail then
					return [true,sol];
				fi;
			fi;
		fi;
	fi;
	return [false];
end);


BindGlobal("LineSegmentIntersectionColinear",function(l1,l2,eps)
	local v1,v2,x1,x2,y1,y2,r1,r2,s1,s2,i,j;
	l1:=1.*l1;
	l2:=1.*l2;
	v1:=l1[2]-l1[1];
	v2:=l2[2]-l2[1];
	if not OnSamePlane(l1,l2,eps) then
		return [false];
	fi;
	if not MyNorm(Crossproduct(v1,v2)) < eps then
		return [false];
	fi;
	if not PointsInOneLine(l1[1],l1[2],l2[1],eps) then
		return [false];
	fi;
	if not PointsInOneLine(l1[1],l1[2],l2[2],eps) then
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

			if MyNorm(x1-x2)<eps then
				#r1=0 and s1=0
				if eps<=r2 and r2 <=1-eps then
					return [true,[[x2,y2],[y2,y1]],5];
				fi;
				if eps<=s2 and s2<=1-eps then
					return [true,[[x1,y1],[y1,y2]],6];
				fi;
			elif MyNorm(y1-x2)<eps then
				#s2=0 and r1=1
				if eps<=r2 and r2 <=1-eps then
					return [true,[[x1,y2],[y2,x2]],7];
				fi;
				if eps<=s1 and s1<=1-eps then
					return [true,[[y1,x1],[x1,y2]],8];
				fi;
			elif MyNorm(y2-x1)<eps then
				#r2=0 and s1=1
				if eps<=r1 and r1 <=1-eps then
					return [true,[[y2,x2],[x2,y1]],9];
				fi;
				if eps<=s2 and s2<=1-eps then
					return [true,[[x2,y1],[y1,x1]],10];
				fi;
			elif MyNorm(y2-y1)<eps then
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

InstallGlobalFunction(calculate_intersections_groups,function(t,coordinates,intersections_only,group,group_orthogonal)
	local l,k,i,j,intersection,intersection2,h,orbs_vof,data_triangulated,orb_rep,vof,inside_points,orbs,g,shift,new_data,cur_line,res,pairs_reps,m,l1,l2,p,n,orb,orbs_pairs;
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
				#check if triangles are coplanar
				n:=Crossproduct(l[i][1][2]-l[i][1][1],l[i][1][3]-l[i][1][1]);
				if (Dot(n,l[j][1][1])-Dot(n,l[i][1][1]))^2<eps and (Dot(n,l[j][1][2])-Dot(n,l[i][1][1]))^2<eps and (Dot(n,l[j][1][3])-Dot(n,l[i][1][1]))^2<eps then
					continue;
					inside_points:=[];
					for m in [1..3] do
						if MyPointInTriangle(l[i][1][1],l[i][1][2],l[i][1][3],l[j][1][m],eps) then
							Add(inside_points,m);
						fi;
					od;
					if Size(inside_points)=0 then
						# no inside points check for line intersections
						for l1 in Combinations([1..3],2) do
							cur_line:=[];
							for l2 in Combinations([1..3],2) do
								res:=LineSegmentIntersection(l[i][1]{l2},l[j][1]{l1},eps);
								if res[1] then
									if Size(res[2])=3 then		
										Add(cur_line,res[2]);
									elif Size(res[2])=2 then
										Error();
										Add(l[i][1],res[2][1]);
										Add(l[i][1],res[2][2]);
										AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
										continue;
									fi;
								fi;
							od;
							if Size(cur_line)=1 then
								Error();
								Add(l[i][1],cur_line[1]);
							elif Size(cur_line)=2 then
								Add(l[i][1],cur_line[1]);
								Add(l[i][1],cur_line[2]);
								AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
							fi;
						od;
					elif Size(inside_points)=1 then
						# one inside point
						Add(l[i][1],l[j][1][inside_points[1]]);
					elif Size(inside_points)=2 then
						Add(l[i][1],l[j][1][inside_points[1]]);
						Add(l[i][1],l[j][1][inside_points[2]]);
					elif Size(inside_points)=3 then
						Add(l[i][1],l[j][1][inside_points[1]]);
						Add(l[i][1],l[j][1][inside_points[2]]);
						Add(l[i][1],l[j][1][inside_points[3]]);
					fi;
					continue;
				fi;
				intersection:=TwoTriangleIntersection(l[i][1]{[1..3]},l[j][1]{[1..3]},eps);
				#test if both intersection and intersection2 are the same
				#otherwise join for lines
				intersection2:=TwoTriangleIntersection(l[j][1]{[1..3]},l[i][1]{[1..3]},eps);
				if (Size(intersection)=1 or Size(intersection)=0) and Size(intersection2)>0 then
					if Size(NumericalUniqueListOfLists(Concatenation(intersection,intersection2),eps))=2 then
						intersection:=NumericalUniqueListOfLists(Concatenation(intersection,intersection2),eps);
					fi;
				fi;
				if Size(intersection)=2 then
					# add intersection informations only to l[i]
					# l[j] will be obtained later
					Add(l[i][1],intersection[1]);
					Add(l[i][1],intersection[2]);
					AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
				fi;
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
	if intersections_only then
		return l;
	else
		#Print("Retriangulation of triangles. \n");
		data_triangulated:=[];
		# transfer triangulations to all other triangles
		for k in [1..Size(orbs)] do
			i:=orbs[k][1];
			#Print("Triangulation for face ",String(i),"\n");
			data_triangulated[i]:=triangulate(my_triangle_fix(l[i]));
			while not IsSimplicialSurface(TriangularComplexByVerticesInFaces(data_triangulated[i][2])) do
				Print("Triangulation of face ",String(i)," failed! Try to triangulate it again.\n");
				data_triangulated[i]:=triangulate(my_triangle_fix(l[i]));
			od;
			for h in [2..Size(orbs[k])] do
				j:=orbs[k][h];
				data_triangulated[j]:=SymmetryActionDirect(data_triangulated[i],i,j,coordinates,group,group_orthogonal);
			od;
		od;
		return join_triangles(data_triangulated);
	fi;
end);

# version without group in orthogonal group
InstallGlobalFunction(calculate_intersections,function(vof,coordinates,intersections_only,group)
	local l,k,i,j,intersection,intersection2,h,orbs_vof,data_triangulated,orb_rep,res,cur_line,inside_points,m,l1,l2,n;
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
				#check if triangles are coplanar and fix cases
				n:=Crossproduct(l[i][1][2]-l[i][1][1],l[i][1][3]-l[i][1][1]);
				if (Dot(n,l[j][1][1])-Dot(n,l[i][1][1]))^2<eps and (Dot(n,l[j][1][2])-Dot(n,l[i][1][1]))^2<eps and (Dot(n,l[j][1][3])-Dot(n,l[i][1][1]))^2<eps then
					inside_points:=[];
					for m in [1..3] do
						if MyPointInTriangle(l[i][1][1],l[i][1][2],l[i][1][3],l[j][1][m],eps) then
							Add(inside_points,m);
						fi;
					od;
					if Size(inside_points)=2 then
						Add(l[i][1],l[j][1][inside_points[1]]);
						Add(l[i][1],l[j][1][inside_points[2]]);
						AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
					elif Size(inside_points)=3 then
						Add(l[i][1],l[j][1][inside_points[1]]);
						Add(l[i][1],l[j][1][inside_points[2]]);
						Add(l[i][1],l[j][1][inside_points[3]]);
						AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
						AddSet(l[i][2],[Size(l[i][1])-2,Size(l[i][1])-1]);
						AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-2]);
					fi;
					#if Size(inside_points)=0 then
						# no inside points check for line intersections
						for l1 in Combinations([1..3],2) do
							cur_line:=[];
							for l2 in Combinations([1..3],2) do
								res:=LineSegmentIntersection(l[i][1]{l2},l[j][1]{l1},eps);
								if res[1] then
									if Size(res[2])=3 then
										if MyNumericalPosition(l[i][1]{[1,2,3]},res[2],eps)=fail then
											Add(cur_line,res[2]);
										fi;
									elif Size(res[2])=2 then
										Error();
										Add(l[i][1],res[2][1]);
										Add(l[i][1],res[2][2]);
										AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
										continue;
									fi;
								fi;
							od;
							if Size(cur_line)=1 then
								#Error();
								#test if point of j from edge l1 lies inside i
									if l1[1] in inside_points and Size(NumericalUniqueListOfLists(Concatenation([l[j][1][l1[1]]],cur_line),eps))=2 then
										Add(l[i][1],cur_line[1]);
										Add(l[i][1],l[j][1][l1[1]]);
										AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
									elif l1[2] in inside_points and Size(NumericalUniqueListOfLists(Concatenation([l[j][1][l1[2]]],cur_line),eps))=2 then
										Add(l[i][1],cur_line[1]);
										Add(l[i][1],l[j][1][l1[2]]);
										AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
									else
										Add(l[i][1],cur_line[1]);
									fi;
							elif Size(cur_line)=2 then
								Add(l[i][1],cur_line[1]);
								Add(l[i][1],cur_line[2]);
								AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
							fi;
						od;
					#elif Size(inside_points)=1 then
						# one inside point
					#	Add(l[i][1],l[j][1][inside_points[1]]);
					#if Size(inside_points)=2 then
					#	Add(l[i][1],l[j][1][inside_points[1]]);
					#	Add(l[i][1],l[j][1][inside_points[2]]);
					#elif Size(inside_points)=3 then
					#	Add(l[i][1],l[j][1][inside_points[1]]);
					#	Add(l[i][1],l[j][1][inside_points[2]]);
					#	Add(l[i][1],l[j][1][inside_points[3]]);
					#fi;
					continue;
				fi;
				intersection:=TwoTriangleIntersection(l[i][1]{[1..3]},l[j][1]{[1..3]},eps);
				#test if both intersection and intersection2 are the same
				#otherwise join for lines
				intersection2:=TwoTriangleIntersection(l[j][1]{[1..3]},l[i][1]{[1..3]},eps);
				if (Size(intersection)=1 or Size(intersection)=0) and Size(intersection2)>0 then
					if Size(NumericalUniqueListOfLists(Concatenation(intersection,intersection2),eps))=2 then
						intersection:=NumericalUniqueListOfLists(Concatenation(intersection,intersection2),eps);
					fi;
				fi;
				if Size(intersection)=2 then
					# add intersection informations only to l[i]
					# l[j] will be obtained later
					Add(l[i][1],intersection[1]);
					Add(l[i][1],intersection[2]);
					AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
				fi;
				if Size(intersection)=1 then
					Add(l[i][1],intersection[1]);
				fi;
			fi;
		od;
	od;
	if intersections_only then
		return l;
	else
		data_triangulated:=[];
		# transfer triangulations to all other triangles
		for k in [1..Size(orbs_vof)] do
			i:=Position(vof,orbs_vof[k][1]);
			data_triangulated[i]:=triangulate(my_triangle_fix(l[i]));
			for h in [2..Size(orbs_vof[k])] do
				j:=Position(vof,orbs_vof[k][h]);
				data_triangulated[j]:=SymmetryAction(data_triangulated[i],vof[i],vof[j],coordinates,group);
			od;
		od;
		return join_triangles(data_triangulated);
	fi;
end);

InstallGlobalFunction(calculate_intersections_comp,function(s,coordinates,intersections_only)
	local l,k,i,j,data,check,c,intersection,intersection2,h,orbs_vof,data_triangulated,orb_rep,inside_points,cur_line,res,m,vof,l2,l1,n;
	vof:=VerticesOfFaces(s);
	l:=[];
	for i in Faces(s) do
		l[i]:=[List([1..3],j->coordinates[vof[i][j]]),Set([Set([1,2]),Set([2,3]),Set([3,1])])];
	od;
	#for each representative find self intersections
	for i in Faces(s) do
		for j in Faces(s) do
			if j in NeighbourFacesOfFace(s,i) then
				inside_points:=[];
					for m in [1..3] do
						if MyPointInTriangle(l[i][1][1],l[i][1][2],l[i][1][3],l[j][1][m],eps) then
							Add(inside_points,m);
						fi;
					od;
				if Size(inside_points)=2 then
					continue;
				fi;
			fi;
			# output of the form [vertices,intersection_edges]
			if i<>j then
				#check if triangles are coplanar
				n:=Crossproduct(l[i][1][2]-l[i][1][1],l[i][1][3]-l[i][1][1]);
				if (Dot(n,l[j][1][1])-Dot(n,l[i][1][1]))^2<eps and (Dot(n,l[j][1][2])-Dot(n,l[i][1][1]))^2<eps and (Dot(n,l[j][1][3])-Dot(n,l[i][1][1]))^2<eps then
					inside_points:=[];
					for m in [1..3] do
						if MyPointInTriangle(l[i][1][1],l[i][1][2],l[i][1][3],l[j][1][m],eps) then
							Add(inside_points,m);
						fi;
					od;
					if Size(inside_points)=0 then
						# no inside points check for line intersections
						for l1 in Combinations([1..3],2) do
							cur_line:=[];
							for l2 in Combinations([1..3],2) do
								res:=LineSegmentIntersection(l[i][1]{l2},l[j][1]{l1},eps);
								if res[1] then
									if Size(res[2])=3 then
										Add(cur_line,res[2]);
									elif Size(res[2])=2 then
										Error();
										Add(l[i][1],res[2][1]);
										Add(l[i][1],res[2][2]);
										AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
										continue;
									fi;
								fi;
							od;
							if cur_line<>[] then
							if Size(NumericalUniqueListOfLists(cur_line,eps))=1 then
								Add(l[i][1],cur_line[1]);
							elif Size(NumericalUniqueListOfLists(cur_line,eps))=2 then
								Add(l[i][1],cur_line[1]);
								Add(l[i][1],cur_line[2]);
								AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
							fi;
						fi;
						od;
					elif Size(inside_points)=1 then
						# one inside point
						if Size(NumericalUniqueListOfLists(Concatenation([l[j][1][inside_points[1]]],l[i][1]{[1,2,3]}),eps))=3 then
							continue;
						fi;
						Add(l[i][1],l[j][1][inside_points[1]]);
						# check edges from inside point
						check:=[1,2,3];
						Remove(check,inside_points[1]);
						for c in check do
							l1:=[inside_points[1],c];
							for l2 in Combinations([1,2,3],2) do
								res:=LineSegmentIntersection(l[i][1]{l2},l[j][1]{l1},eps);
								if res[1] then
									Add(l[i][1],res[2]);
								fi;
							od;
						od;
						if Size(NumericalUniqueListOfLists(l[i][1]{[Size(l[i][1]),Size(l[i][1])-2]},eps))=2 then
							AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-2]);
						fi;
						if Size(NumericalUniqueListOfLists(l[i][1]{[Size(l[i][1])-1,Size(l[i][1])-2]},eps))=2 then
							AddSet(l[i][2],[Size(l[i][1])-1,Size(l[i][1])-2]);
						fi;
						# check other edge
						l1:=[1,2,3];
						Remove(l1,inside_points[1]);
						cur_line:=[];
						for l2 in Combinations([1..3],2) do
							res:=LineSegmentIntersection(l[i][1]{l2},l[j][1]{l1},eps);
							if res[1] then
								if Size(res[2])=3 then	
									Add(cur_line,res[2]);
								elif Size(res[2])=2 then
									Add(l[i][1],res[2][1]);
									Add(l[i][1],res[2][2]);
									if Size(NumericalUniqueListOfLists(l[i][1]{[Size(l[i][1])-1,Size(l[i][1])]},eps))=2 then
										AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
									fi;
									continue;
								fi;
							fi;
						od;
						if cur_line<>[] then
							if Size(NumericalUniqueListOfLists(cur_line,eps))=1 then
								Add(l[i][1],cur_line[1]);
							elif Size(NumericalUniqueListOfLists(cur_line,eps))=2 then
								Add(l[i][1],cur_line[1]);
								Add(l[i][1],cur_line[2]);
								AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
							fi;
						fi;
					elif Size(inside_points)=2 then
						Add(l[i][1],l[j][1][inside_points[1]]);
						Add(l[i][1],l[j][1][inside_points[2]]);
						AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
						check:=[1,2,3];
						Remove(check,inside_points[1]);
						Remove(check,inside_points[2]);
						for c in check do
							l1:=[inside_points[1],c];
							for l2 in Combinations([1,2,3],2) do
								res:=LineSegmentIntersection(l[i][1]{l2},l[j][1]{l1},eps);
								if res[1] then
									Add(l[i][1],res[2]);
								fi;
							od;
						od;
						if Size(NumericalUniqueListOfLists(l[i][1]{[Size(l[i][1]),Size(l[i][1])-2]},eps))=2 then
							AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-2]);
						fi;
						check:=[1,2,3];
						Remove(check,inside_points[1]);
						Remove(check,inside_points[2]);
						for c in check do
							l1:=[inside_points[2],c];
							for l2 in Combinations([1,2,3],2) do
								res:=LineSegmentIntersection(l[i][1]{l2},l[j][1]{l1},eps);
								if res[1] then
									Add(l[i][1],res[2]);
								fi;
							od;
						od;
						if Size(NumericalUniqueListOfLists(l[i][1]{[Size(l[i][1]),Size(l[i][1])-2]},eps))=2 then
							AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-2]);
						fi;
					elif Size(inside_points)=3 then
						Add(l[i][1],l[j][1][inside_points[1]]);
						Add(l[i][1],l[j][1][inside_points[2]]);
						Add(l[i][1],l[j][1][inside_points[3]]);
						AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
						AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-2]);
						AddSet(l[i][2],[Size(l[i][1])-2,Size(l[i][1])-1]);
					fi;
					continue;
				fi;
				intersection:=TwoTriangleIntersection(l[i][1]{[1..3]},l[j][1]{[1..3]},eps);
				#test if both intersection and intersection2 are the same
				#otherwise join for lines
				intersection2:=TwoTriangleIntersection(l[j][1]{[1..3]},l[i][1]{[1..3]},eps);
				if (Size(intersection)=1 or Size(intersection)=0) and Size(intersection2)>0 then					
					if Size(NumericalUniqueListOfLists(Concatenation(intersection,intersection2),eps))=2 then
						intersection:=NumericalUniqueListOfLists(Concatenation(intersection,intersection2),eps);
					fi;
				fi;
				if Size(intersection)=2 then
					# add intersection informations only to l[i]
					# l[j] will be obtained later
					Add(l[i][1],intersection[1]);
					Add(l[i][1],intersection[2]);
					AddSet(l[i][2],[Size(l[i][1]),Size(l[i][1])-1]);
				fi;
			fi;
		od;
	od;
	data_triangulated:=[];
	for i in Faces(s) do
		#Print(i,"\n");
		data_triangulated[i]:=triangulate(my_triangle_fix(l[i]));
	od;
	data:=[];
	for c in ConnectedComponents(s) do
		Add(data,join_triangles(data_triangulated{Faces(c)}));
	od;
	return data;
end);

# Clean the data by removing duplicate vertices and edges
BindGlobal("clean_data", function(data, eps)
    local map, i, new_edges, new_vertices, map2;
    data[1] := 1. * data[1]; # Ensure vertices are in floating-point format
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
    new_edges := List(data[2], j -> [map2[j[1]], map2[j[2]]]);
    data[2] := new_edges;
    data[1] := new_vertices;
    # Delete multiple occurrences of edges
    data[2] := Set(List(data[2], i -> Set(i)));
    return data;
end);




InstallGlobalFunction(fix_intersections_planar,function(data)
	local i,j,res,entry,entry1,entry2,entry_in_i,entry_in_j,addy,check_edges,orig_j,orig_i,cur,k,remove,e,l,n;
	data:=clean_data(data,eps);
	data[2]:=List(data[2]);
	for i in [1..Size(data[1])] do
		check_edges:=ShallowCopy(data[2]);
		while check_edges<>[] do
			e:=Remove(check_edges);
			if LineSegmentIntersectionColinear([data[1][i],data[1][e[1]]],[data[1][e[1]],data[1][e[2]]],eps)[1] and LineSegmentIntersectionColinear([data[1][i],data[1][e[2]]],[data[1][e[1]],data[1][e[2]]],eps)[1] and not i in e then
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
	data:=clean_data(data,eps);
	data[2]:=List(data[2]);
	check_edges:=Combinations([1..Size(data[2])],2);;
	while check_edges<>[] do
		cur:=Remove(check_edges,1);
		i:=cur[1];
		j:=cur[2];
		if i <=Size(data[2]) and j<=Size(data[2]) and i in BoundPositions(data[2]) and j in BoundPositions(data[2]) then
			res:=LineSegmentIntersection([data[1][data[2][i][1]],data[1][data[2][i][2]]],[data[1][data[2][j][1]],data[1][data[2][j][2]]],eps);
			if res[1] then
				entry:=MyNumericalPosition(data[1],res[2],eps);
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
	data:=clean_data(data,eps);
	data[2]:=List(data[2]);
	for i in [1..Size(data[1])] do
		check_edges:=ShallowCopy(data[2]);
		while check_edges<>[] do
			e:=Remove(check_edges);
			if LineSegmentIntersectionColinear([data[1][i],data[1][e[1]]],[data[1][e[1]],data[1][e[2]]],eps)[1] and LineSegmentIntersectionColinear([data[1][i],data[1][e[2]]],[data[1][e[1]],data[1][e[2]]],eps)[1] and not i in e then
				orig_j:=ShallowCopy(e);
				Remove(data[2],Position(data[2],orig_j));
				if not Set([e[1],i]) in data[2] or not Set([e[2],i]) in data[2] then
					data[2]:=Concatenation(data[2],Set([Set([e[1],i]),Set([e[2],i])]));
					check_edges:=Concatenation(check_edges,Set([Set([e[1],i]),Set([e[2],i])]));
				fi;
			fi;
		od;

	od;
	return clean_data(data,eps);;
end);


BindGlobal("triangulate_comb",function(data)
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