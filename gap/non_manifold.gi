Read("triangulation.gd");
Read("triangulation.gi");
Read("my_triangle_fix.g"); 

# now shift the right faces
# find faces with same vertex v as f (3)
UmbrellaPathFaceVertex:=function(t,f,v,data,points)
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
end;;

UmbrellaComponents := function(t,v)
	local u_comps, v_faces, v_edges, edges_of_faces, f, e, fa, edges, avaiable_edges, comp, connection, cur_edges;;
	
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
end;;

ChooseIsolatedVertex:=function(t,e,points,data)
	local voe,foe;
	voe:=VerticesOfEdge(t,e);
	foe:=FacesOfEdge(t,e);
	if IsSubset(Set(foe),Set(UmbrellaPathFaceVertex(t,foe[1],voe[1],data,points))) then
		return voe[1];
	fi;
	return voe[2];
end;;


FixNormals:=function(Coords)
    local f, faces, norm, c_verts, eps;
    eps := 1./10^6;
    for f in Coords do
        c_verts := ShallowCopy([f[1],f[2],f[3]]);
        norm := Crossproduct(c_verts[1]-c_verts[2],c_verts[1]-c_verts[3])/Norm2(Crossproduct(c_verts[1]-c_verts[2],c_verts[1]-c_verts[3]));
        if FlGeq(norm*f[4],0.,eps) then
        	f[4] := norm;
        else
        	f[4] := -norm;
        fi; 
    od;
    return Coords;
end;;



ReadSTL:=function(fileName)
	# reads a file from the current dir
	local surf, file, name, r, r2, eps, filesepr, endsign, normal, data, normals, points, test, i,j, index, verts, coords, input, Coords;
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
	data := SimplicialSurfaceFromCoordinates([Coords,faces],eps);
	return [data[1],Coords,points];
end;


ConvertDataFormatPtC:=function(surf,points,normals)
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
end;;


ConvertDataFormatCtP:=function(surf,Coords)
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
end;;

FixVertOfOuterNMEdge:=function(t,e,Coords,points,data,shift_param,not_split)
	local v,data_fan,data_fix,alpha,verts_e, vC,n,v_p, p_f, n_f, n_ff,f, ff, vec,v_alpha,w,index_f, verts_f, verts_ff, index_ff,points_fix,t_fix,MyNormal,VoE,VoF,fa, u, p, v_index,w_index,v_faces,w_faces,i;
	f:=FacesOfEdge(t,e)[1];

	
	data_fan:=CalculateFan(t,e,points,f,data[4][f]);
    

	index_f:=0;
	for i in [1..Size(data_fan)] do
		if data_fan[i][1]=f then
			index_f:=i;
		fi;
	od;
	index_ff:=(index_f mod (Size(data_fan))) + 1;

    	ff := data_fan[index_ff][1];
    	verts_e := ShallowCopy(VerticesOfEdge(t,e));
    	verts_f := VerticesOfFace(t,f);       
    	verts_ff := VerticesOfFace(t,ff);

    	# have to check direction -> so pick correct faces to change
	
    	n_f := Filtered(verts_f,x-> not x in verts_e)[1];
   	n_ff := Filtered(verts_ff,x-> not x in verts_e)[1];
	
    	v := verts_e[1];
    	if v in not_split then
    		v := verts_e[2];
    	fi;
    	
    	points_fix:=points;
    	
	p := points[v];
	v_p := p + 0.5*shift_param*(points[n_f]-p) + 0.5*shift_param*(points[n_ff]-p);

	    	
	Add(points_fix,v_p);
	v_index:=Size(points_fix);

	    
    	# need to change all verts on umbrella
	    	
	    
    	u := UmbrellaPathFaceVertex(t,f,v,data,points);

    	for fa in u do
    		if v in Coords[fa][5] then
		    	p_f := Position(Coords[fa][5],v);
		    	Coords[fa][p_f] := v_p;
		    	Coords[fa][5][p_f] := v_index;
		fi;
    	od;
    	
	return [Coords,points_fix,[v]];
end;;

FixNMVert:= function(surf,p,Coords,points,data,shift_param)
	local coord_p,f, fp, pu, verts, umbrella_verts, points_fix, umb, p_index, coord_p_new, vec, p_f, fa, u, i;
	
	# fix non manifold vertex
	
	coord_p := points[p];
	umb := UmbrellaComponents(surf,p);
	while umb <> [] do 
		umbrella_verts := [];
		u := umb[1];
		f := u[1];
		
		i := 1;
		
		# find vertices in local umbrella that are not p
		for fp in u do
			verts := ShallowCopy(VerticesOfFace(surf,fp));
			Remove(verts,p);
			umbrella_verts[i] := verts;
			i := i + 1;
		od;
		umbrella_verts := DuplicateFreeList(Flat(umbrella_verts));
		
		# take averaged direction of these vertices
		vec := [0.,0.,0.];
		for pu in umbrella_verts do
			vec := vec - coord_p + points[pu];
		od;
		vec := 1./Length(umbrella_verts)*shift_param*vec;
		coord_p_new := coord_p + vec;	
		Add(points,coord_p_new);
		p_index := Length(points);
		
		# Change the coords of all vertices on the local umbrella to remedy p
		
		
		for fa in u do
	    		if p in Coords[fa][5] then
			    	p_f := Position(Coords[fa][5],p);
			    	Coords[fa][p_f] := coord_p_new;
			    	Coords[fa][5][p_f] := p_index;
			fi;
		od;
		
		
		
		Remove(umb,Position(umb,u));
	od;
	points_fix := ShallowCopy(points);
	return [Coords,points_fix,[p]];
end;;


FixVertOfInnerNMEdge:=function(t,e,Coords,points,data,shift_param,not_split)
	local v,data_fan,data_fix, verts_e, n_f, n_ff, p_f, v_p, fa, alpha,vC,n,vec, verts_f, verts_ff, v_alpha,f,w,u,p, ff, initial_f, initial_ff, normal, index_f,index_ff,points_fix,t_fix,MyNormal,VoE,VoF,v_index,w_index,v_faces,w_faces,i;
	initial_f:=FacesOfEdge(t,e)[1];

	#halber winkel
	data_fan:=CalculateFan(t,e,points,initial_f,data[4][initial_f]);
    

	index_f:=0;
	for i in [1..Size(data_fan)] do
		if data_fan[i][1]=initial_f then
			index_f:=i;
		fi;
	od;
	index_ff:=(index_f mod (Size(data_fan))) + 1;

    	initial_ff := data_fan[index_ff][1];
    	verts_e := VerticesOfEdge(t,e);
    	

    	# we have to check direction
	
	normal:=data[4][index_f];
	n:=points[verts_e[2]]-points[verts_e[1]];
	n:=n/MyNorm(n);
	vec:=data_fan[index_f][3];
    	
	
    	v := verts_e[1];
    	if v in not_split then
    		v := verts_e[2];
    	fi;
    	
    	points_fix:=points;
	
	
	if Determinant([vec,n,normal]) < 0. then
			f := initial_f;
			ff := initial_ff;
	else
		f := initial_f;
		ff := initial_ff;
	fi;
	p := points[v];
	
	verts_f := VerticesOfFace(t,f);       
	verts_ff := VerticesOfFace(t,ff);
	n_f := Filtered(verts_f,x-> not x in verts_e)[1];
	n_ff := Filtered(verts_ff,x-> not x in verts_e)[1];
	  	
	v_p := p + 0.5*shift_param*(points[n_f]-p) + 0.5*shift_param*(points[n_ff]-p);
	

	Add(points_fix,v_p);
	v_index:=Length(points_fix);
	    
    	# need to change all verts on umbrella path
    	
    	u := UmbrellaPathFaceVertex(t,f,v,data,points);

    	for fa in u do
    		if v in Coords[fa][5] then
		    	p_f := Position(Coords[fa][5],v);
		    	Coords[fa][p_f] := v_p;
		    	Coords[fa][5][p_f] := v_index;
		fi;
	od;

    	

	return [Coords,points_fix,[v]];
end;;



InnerRemEdge:= function(e,edges,verts_of_es)
	local in1, in2, vert1, vert2, q;
	in1 := false;
	in2 := false;
	vert1 := verts_of_es[e][1];
	vert2 := verts_of_es[e][2];
	for q in edges do
		if q <> e then
			if vert1 in verts_of_es[q] then
				in1 := true;
			elif vert2 in verts_of_es[q] then
				in2 := true;
			fi;
		fi;
	od;
	if in1 and in2 then
		return true;
	else
		return false;
	fi;
end;;

# only works with simple paths
ChooseStartOfNMPath:=function(e,verts_of_e)
	local delet_e, inner, l, q, delet1, delet2, vert1, vert2;
    inner := [];
    i := 1;
	delet_e := ShallowCopy(e);
	for l in e do
        inner[l] := false;
		if InnerRemEdge(l,e,verts_of_e) then
			Remove(delet_e,Position(delet_e,l));
            inner[l] := true;
		fi;
        i := i+1;
	od;
	
	Flat(delet_e);
	if delet_e = [] then
		# NM edges form circle
		return [e,true,inner];
	else
		return [delet_e,false,inner];
	fi;
end;;

DetectCircle:=function(Edges,VertsOfEdges,start,next)
	local l, traversedVerts, leftVerts, leftEdges, reset, start_e, cur_e, cur_vs, n_v;
	reset := false;
	traversedVerts := [];
	leftVerts := Unique(Flat(ShallowCopy(VertsOfEdges)));
	leftEdges := ShallowCopy(Edges);
	start_e := start[1];
	cur_e := ShallowCopy(start_e);
	cur_vs := start[2];
	
	while leftVerts <> [] and leftEdges <> [] do
		for l in Edges do
			n_v := Intersection(cur_vs,VertsOfEdges[l],leftVerts);
			if n_v <> [] and l <> cur_e then
				if l = start_e then
					return true;
				fi;
				
				Remove(leftVerts,Position(leftVerts,n_v[1]));
				if cur_e in leftEdges then
					Remove(leftEdges,Position(leftEdges,cur_e));
				fi;
				cur_e := l;
				cur_vs := VertsOfEdges[l];
				reset := false;
			fi;
		od;
		if reset then
			return false;
		fi;
		l := start_e;
		reset := true;
	od;
	
	return false;
end;; 

OrderPath:=function(NMEdges,VertsOfNMEdges,start,path)
	local found, cur_verts, coincident_vert, circle_in_direction, q, c_numb, p, l, i, numb, next, next_verts, next_edges, next_edge_verts, comp;
	numb := 0;
	c_numb := 0;
	next := [];
	l := start[1];
	cur_verts := start[2];
	coincident_vert := start[3];
	
	# find next edges, make sure to not go backwards
	# TODO: fix circles
	for q in NMEdges do
		if ((cur_verts[1] in VertsOfNMEdges[q] and cur_verts[1] <> coincident_vert) or (cur_verts[2] in VertsOfNMEdges[q] and cur_verts[2] <> coincident_vert)) and not q in Flat(path) and not q = l then
				
				
				circle_in_direction := DetectCircle(NMEdges,VertsOfNMEdges,start,q);
				if not (circle_in_direction and c_numb > 0) then
					numb := numb + 1;
					if circle_in_direction then
						c_numb := c_numb + 1;
					fi;
					next[numb] := q;
				fi;
				
		fi;
		
	od;
	# check along which vert we are going
	if numb > 0 then
		if cur_verts[1] in VertsOfNMEdges[next[1]] then
			coincident_vert := cur_verts[1];
		else
			coincident_vert := cur_verts[2];
		fi;
	else
		coincident_vert := 0;
	fi;
	i := 0;
	comp := [];
	if numb = 1 then
		p := next[1];
		next_edges := ShallowCopy(Flat(DuplicateFreeList(NMEdges)));
		Unbind\[\](next_edges,Position(next_edges,p));
			
		next_verts := VertsOfNMEdges[p];
		next_edge_verts := ShallowCopy(VertsOfNMEdges);
		Unbind\[\](next_edge_verts,p);
		
		path[Length(path)+1] := l;
		OrderPath(next_edges,next_edge_verts,[p,next_verts,coincident_vert],path);
	elif numb = 0 then
		path[Length(path)+1] := l;
	else
		path[Length(path)+1] := l;
		for p in next do
			i := i + 1;
			
			next_edges := ShallowCopy(Flat(DuplicateFreeList(NMEdges)));
			Unbind\[\](next_edges,Position(next_edges,p));
			
			next_verts := VertsOfNMEdges[p];
			next_edge_verts := ShallowCopy(VertsOfNMEdges);
			Unbind\[\](next_edge_verts,p);
			
			if next_edges <> [] then
				comp[i] := OrderPath(next_edges,next_edge_verts,[p,next_verts,coincident_vert],[]);
			else
				comp[i] := [p];
			fi;
		od;
		
		if Length(Filtered(Flat(ShallowCopy(comp)),x-> not x in path)) = Length(Flat(ShallowCopy(comp))) then
			path[Length(path)+1] := comp;
		fi;
	fi;
	
	return path;
end;;




OrderNMEdges:=function(surf, data)
	local e,VertsOfNMEdges, coincident_vert, path, isolated_edges, cur_path, found, l, q, i, cur_verts, info, pos, curr_verts;
	e:=ShallowCopy(RamifiedEdges(surf));
	VertsOfNMEdges := [];
	
	for l in e do
		VertsOfNMEdges[l] := VerticesOfEdge(surf,l);
	od;
	
	info := ChooseStartOfNMPath(e,VertsOfNMEdges);
	
	l := info[1][1];
	isolated_edges := [];
	path := [];
	cur_verts := [];
	i := 1;
	j := 1;
	cur_path := [l];
	coincident_vert := 0;
	
	while e <> [] do
		found := false;
		cur_verts := ShallowCopy(VertsOfNMEdges[l]);
		cur_path := OrderPath(e,VertsOfNMEdges,[l,cur_verts,0],[]);
		path[i] := ShallowCopy(cur_path);
		
		cur_path := Flat(cur_path);
		for k in cur_path do
			Remove(e,Position(e,k));
		od;
		if e <> [] then
			l := e[1];
		fi;
		i := i + 1;
	od;
    return [path,info[2],info[3]];
end;;


FixNMIntersection := function(surf,order,info,data,points,Coords,shift_param)
	local l, int_v, branches, not_split, data_fix, b;
	int_v := info[1];
	branches := info[2];
	not_split := [];
	for b in branches do
		data_fix := FixVertOfOuterNMEdge(surf,b[1],Coords,points,data,shift_param,not_split);
		Coords:=data_fix[1];
		points:=data_fix[2];
		if Length(b[2]) > 1 then
			# the path continues in this direction
			data_fix := FixNMPathRec(surf,[b[2],order[2],order[3]],data,points,Coords,shift_param);
			points := data_fix[1];
			Coords := data_fix[2];
		fi;
	od;
	return [points,Coords];
end;;

FixNMPathRec := function(surf,order,data,points,Coords,shift_param)
	local l, comp, path, not_split, data_fix, inner, is_circle, points_fix, s_data, branches, branch, verts_current, verts_comp, same, info,j, e, v,len, s,t;

	is_circle := order[2];
	inner := order[3];
	for comp in order[1] do
		path := comp;
		not_split := [];
		if IsInt(path) then
			if path in Edges(surf) then
				if inner[path] or is_circle then			
						data_fix:=FixVertOfInnerNMEdge(surf,path,Coords,points,data,shift_param,not_split);
						Coords:=data_fix[1];
						points:=data_fix[2];
						not_split := data_fix[3];
				else
						verts_current := VerticesOfEdge(surf,path);
						
						same := ChooseIsolatedVertex(surf,path,points,data);
						Add(not_split,same);
						data_fix := FixVertOfOuterNMEdge(surf,path,Coords,points,data,shift_param,not_split);
			      			Coords:=data_fix[1];
						points:=data_fix[2];
			
				fi;
			else 
				Print("error");
			fi;
			
		else
			# recursive case
			for i in [1..Length(path)] do
				e:=path[i];
				len := Length(path);
				
				# if l>1 then we are at a crossing
				if len > 1 then
					v := Intersection(VerticesOfEdge(surf,e[1]),VerticesOfEdge(surf,path[2][1]))[1];
					branches := [];
					
					for j in [1..len] do
						branch := Flat(ShallowCopy(path[j]))[1];
						branches[j]:= [branch,path[j]];
					od;
					info := [v,branches];
					data_fix := FixNMIntersection(surf,[e,is_circle,inner],info,data,points,Coords,shift_param);
					points := data_fix[1];
					Coords := data_fix[2];
				else
					if not IsInt(e) then
						e := e[1];
					fi;
					if e in Edges(surf) then
						if inner[e] or is_circle then
							data_fix:=FixVertOfInnerNMEdge(surf,e,Coords,points,data,shift_param,not_split);
							Coords:=data_fix[1];
							points:=data_fix[2];
							not_split := [];
						elif Length(path) = 1 then
							verts_current := VerticesOfEdge(surf,e);
							
							same := ChooseIsolatedVertex(surf,e,points,data);
							data_fix := FixVertOfOuterNMEdge(surf,e,Coords,points,data,shift_param,not_split);
				      			Coords:=data_fix[1];
							points:=data_fix[2];
						else
						    	verts_current := VerticesOfEdge(surf,e);
						    	if i = 1 then
						     		verts_comp :=  VerticesOfEdge(surf,Flat(ShallowCopy(path))[i+1]);
						  	else
						   		verts_comp :=  VerticesOfEdge(surf,Flat(ShallowCopy(path))[i-1]);
							fi;
						 	# vertex which does not need to be split
						 	same := Filtered(verts_current,x-> not x in verts_comp)[1];
							data_fix := FixVertOfOuterNMEdge(surf,e,Coords,points,data,shift_param,not_split);
				      			Coords:=data_fix[1];
							points:=data_fix[2];
							not_split := data_fix[3];
						fi;
					else
						Print("error");
					fi;
					
				fi;
		    	od;
	    	fi;
	    	
	od;
	points_fix := points;
    	return [points_fix,Coords];
end;;


FixNMEdgePath := function(surf,order,data,points,Coords,shift_param)
	local l, comp, path, not_split, data_fix, inner, is_circle, NM_verts, points_fix, s_data, verts_current, verts_comp, same, e, s,t; 
	# TODO: calculate order data here
	not_split := [];
	is_circle := order[2];
	inner := order[3];
	
	for comp in order[1] do
		path := comp;
		data_fix := FixNMPathRec(surf,[comp,is_circle,inner],data,points,Coords,shift_param);
		points := data_fix[1];
		Coords:= data_fix[2];
		
	od;
	
	s_data := SimplicialSurfaceFromChangedCoordinates([Coords,surf],1./10^6);
	surf := s_data[1];
	points_fix := points;
	
	
    	return [surf,points_fix,Coords];
end;;



FixNMVerts := function(surf, data, points, Coords, shift_param)
	local NM_verts, v, new_data, s_data;
	
	# to be executed after fixing the non-manifold edges
	
	NM_verts := ShallowCopy(RamifiedVertices(surf));
	
	for v in NM_verts do
		new_data := FixNMVert(surf,v,Coords,points,data,shift_param);
		Coords := new_data[1];
		points := new_data[2];
		
	od;
	
	s_data := SimplicialSurfaceFromChangedCoordinates([Coords,surf],1./10^6);
	surf := s_data[1];
	
	Coords := FixNormals(Coords);
	return [surf,points,Coords];
end;;


Remedy_NonManifold := function(data,points, shift_param)
	local surf, shift_param, normals, Coords, order_data, m_data, m_surf, m_points, m_coords, fully_m_data;
	#
	# input structure is data=[surf1,surf2,faces,normals_coordinates] (output of outer hull function), points: coordinates of the points of surf2
	# surf1 is a intersection free complex 
	# surf2 its outer hull
	# faces are the faces of surf2
	# normals_coordinates are the coordinates of the normals of the faces of surf2, indexed by the resp. face number
	#
	# output has structure [fully_m_surf, fully_m_points, fully_m_coords] where fully_m_surf is the manifold version of surf2 and the other two its data
	
	
	
	surf := ShallowCopy(data[2]);
	normals := ShallowCopy(data[4]);
	
	Coords := ConvertDataFormatPtC(surf,points,normals);
	
	order_data := OrderNMEdges(surf,ShallowCopy(data));
	
	m_data := FixNMEdgePath(surf,order_data,ShallowCopy(data),points,Coords,shift_param);
	
	m_surf := m_data[1];
	m_points := m_data[2];
	m_coords := m_data[3];


	fully_m_data := FixNMVerts(m_surf,data,ShallowCopy(m_points),m_coords,shift_param);
	
	# TODO: output order data?
	return fully_m_data;
end;;


##############################################################################################
##########
########
######## exemplary use case
Coord3_1:= [
                        [  0.5877852523,  0.0000000000,  0.0000000000 ],
                        [ -0.2628655561,  0.5257311121,  0.0000000000 ],
                        [ -0.2628655561, -0.4253254041,  0.3090169943 ],
                        [  0.4979796570,  0.8057480107,  0.5854101964 ],
                        [ -0.2628655561,  0.1624598481,  0.4999999999 ],
                        [ -0.2628655561,  0.1624598481, -0.4999999999 ],
                        [  0.2628655561, -0.1624598481, -0.4999999999 ],
                        [  0.2628655561, -0.1624598481,  0.4999999999 ],
                        [  0.2628655561,  0.4253254041, -0.3090169943 ],
                        [  0.2628655561, -0.5257311121,  0.0000000000 ],
                        [ -0.4979796570, -0.8057480107, -0.5854101964 ],
                        [ -0.5877852523,  0.0000000000,  0.0000000000 ],
                        ];
#### this block should be handeled by the previous functions  
#
                   
calculate_intersections(VerticesOfFaces(Icosahedron()),Coord3_1,false,Group(()));
datas:=last;
points:=datas[1];
points_fix := ShallowCopy(points);
t:=TriangularComplexByVerticesInFaces(datas[2]); 
VerticesOfFace(t,1); 
n:=Crossproduct(points[2]-points[1],points[3]-points[1]);
n:=-n/Norm2(n);
datas:=OuterHull(t,points,1,-n);

#
####

shift_param:=0.04;
f := Remedy_NonManifold(datas,points_fix);
DrawSTLwithNormals(f[1],"ico_3_1_path_repaired",f[2],datas[4],[]);




