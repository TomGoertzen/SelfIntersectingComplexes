# output umbrella around vertex v (with possible non-manifold edges present)
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
			new_face:=_UpwardContinuation(t,edge,points,new_face,data[4][new_face])[1];
		else
			new_face:=Difference(FacesOfEdge(t,edge),[new_face])[1];
		fi;
		edge:=Difference(Intersection(EdgesOfVertex(t,v),EdgesOfFace(t,new_face)),EdgesOfFace(t,old_face))[1];
	od;
	return faces;
end);;

# outputs local umbrellas at vertex v
# TODO: what should output be exactly?
BindGlobal("UmbrellaComponents",function(t,v)
	local u_comps, available_edges, v_faces, v_edges, edges_of_faces, f, e, fa, edges, avaiable_edges, comp, connection, cur_edges, used_edges;;
	
	u_comps := [];
	v_faces := ShallowCopy(FacesOfVertex(t,v));
	v_edges := ShallowCopy(EdgesOfVertex(t,v));
	edges_of_faces := [];
	used_edges := [];
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
				if connection <> [] and not connection in used_edges then
					used_edges := Concatenation(used_edges,connection);
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

BindGlobal("_ChooseIsolatedVertex", function(t,e,path,points,data)
	local voe,foe, at_end, pos, next, vonext;
	voe:=VerticesOfEdge(t,e);
	foe:=FacesOfEdge(t,e);
	pos := Position(path,e);
	at_end := false;
	if Length(path) = 1 then
		Print("error");
	else
		if not pos = Length(path) then
			next := path[pos+1];
			
		else
			next := path[pos-1];
			at_end := true;
		fi;
		while not IsInt(next) do
				next := next[1];
		od;
	fi;
	vonext := VerticesOfEdge(t,next);
	
	return [Difference(voe,Intersection(vonext,voe))[1],at_end];
end);;


BindGlobal("_UpdateNormals", function(Coords)
    local f, faces, norm, c_verts, eps;
    eps := 1./10^6;
    for f in Coords do
        c_verts := ShallowCopy([f[1],f[2],f[3]]);
        norm := _Crossproduct(c_verts[1]-c_verts[2],c_verts[1]-c_verts[3])/_MyNorm(_Crossproduct(c_verts[1]-c_verts[2],c_verts[1]-c_verts[3]));
        if _FlGeq(norm*f[4],0.,eps) then
        	f[4] := norm;
        else
        	f[4] := -norm;
        fi; 
    od;
    return Coords;
end);;



BindGlobal("_FixVertOfOuterNMEdge", function(t,e,Coords,points,data,shift_param,not_split)
	local v,data_fan,data_fix,alpha,verts_e, vC,n,v_p, p_f, n_f, n_ff,f, ff, vec,v_alpha,w,index_f, verts_f, verts_ff, index_ff,points_fix,t_fix,MyNormal,VoE,VoF,fa, u, p, v_index,w_index,v_faces,w_faces,i;
	f:=FacesOfEdge(t,e)[1];

	
	data_fan:=_CalculateFan(t,e,points,f,data[4][f]);
    

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

	    
    	# need to change all verts on (partial) umbrella
    	# TOOD: does this fail if we have three or more partial umbrellas?
	    	
	    
    	u := UmbrellaPathFaceVertex(t,f,v,data,points);

    	for fa in u do
    		if v in Coords[fa][5] then
		    	p_f := Position(Coords[fa][5],v);
		    	Coords[fa][p_f] := v_p;
		    	Coords[fa][5][p_f] := v_index;
		fi;
    	od;
    	
	return [Coords,points_fix,[v]];
end);;

BindGlobal("_FixNMVert", function(surf,p,Coords,points,data,shift_param)
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
end);;


BindGlobal("_FixVertOfInnerNMEdge", function(t,e,Coords,points,data,shift_param,not_split)
	local v,data_fan,data_fix, verts_e, n_f, n_ff, p_f, v_p, fa, alpha,vC,n,vec, verts_f, verts_ff, v_alpha,f,w,u,p, ff, initial_f, initial_ff, normal, index_f,index_ff,points_fix,t_fix,MyNormal,VoE,VoF,v_index,w_index,v_faces,w_faces,i;
	initial_f:=FacesOfEdge(t,e)[1];

	#halber winkel
	data_fan:=_CalculateFan(t,e,points,initial_f,data[4][initial_f]);
    

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
	n:=n/_MyNorm(n);
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
end);;



BindGlobal("_InnerRemEdge", function(e,edges,verts_of_es)
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
end);;

# only works with simple paths
BindGlobal("_ChooseStartOfNMPath", function(e,verts_of_e)
	local delet_e, inner, i,l, q, delet1, delet2, vert1, vert2;
    	inner := [];
   	i := 1;
	delet_e := ShallowCopy(e);
	for l in e do
        inner[l] := false;
		if _InnerRemEdge(l,e,verts_of_e) then
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
end);;

BindGlobal("_DetectCircle", function(Edges,VertsOfEdges,start,next)
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
end);; 

InstallGlobalFunction(_OrderPath, function(NMEdges,VertsOfNMEdges,start,path)
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
				
				
				circle_in_direction := _DetectCircle(NMEdges,VertsOfNMEdges,start,q);
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
		_OrderPath(next_edges,next_edge_verts,[p,next_verts,coincident_vert],path);
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
				comp[i] := _OrderPath(next_edges,next_edge_verts,[p,next_verts,coincident_vert],[]);
			else
				comp[i] := [p];
			fi;
		od;
		
		if Length(Filtered(Flat(ShallowCopy(comp)),x-> not x in path)) = Length(Flat(ShallowCopy(comp))) then
			path[Length(path)+1] := comp;
		fi;
	fi;
	
	return path;
end);;




InstallGlobalFunction(OrderNMEdges, function(surf, data)
	local e,VertsOfNMEdges, coincident_vert, path, isolated_edges, cur_path, found, l, q, i, j, k, cur_verts, info, pos, curr_verts;
	e:=ShallowCopy(RamifiedEdges(surf));
	VertsOfNMEdges := [];
	
	for l in e do
		VertsOfNMEdges[l] := VerticesOfEdge(surf,l);
	od;
	
	info := _ChooseStartOfNMPath(e,VertsOfNMEdges);
	
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
		cur_path := _OrderPath(e,VertsOfNMEdges,[l,cur_verts,0],[]);
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
end);;

#34 outer  [8,9]
#37 inner  [9,10]
#28 inner 
#29 inner 
#45 inner 
#42 outer

BindGlobal("_FixNMPathRec", function(surf,order,data,points,Coords,shift_param)
	local l, comp, path, not_split, data_fix, inner, is_circle, points_fix, s_data, branches, branch, verts_current, verts_comp, same, info,j, e, v, i, len, s,t;
	# wrong somewhere: not_split is issue
	is_circle := order[2];
	inner := order[3];
	not_split := [];
	#Print(order[1]);
	for comp in order[1] do
		path := comp;
		
		if IsInt(path) then
			if path in Edges(surf) then
				if inner[path] or is_circle then
						data_fix:=_FixVertOfInnerNMEdge(surf,path,Coords,points,data,shift_param,not_split);
						Coords:=data_fix[1];
						points:=data_fix[2];
						not_split := data_fix[3];
				else
						verts_current := VerticesOfEdge(surf,path);
						
						same := _ChooseIsolatedVertex(surf,path,order[1],points,data);
						if not same[2] then
							Add(not_split,same[1]);
							data_fix := _FixVertOfOuterNMEdge(surf,path,Coords,points,data,shift_param,not_split);
							#Print("split");							
				      		Coords:=data_fix[1];
							points:=data_fix[2];
							not_split := data_fix[3];
						else
							#Print("End\n");
						fi;
			
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
					data_fix := _FixNMIntersection(surf,[e,is_circle,inner],info,data,points,Coords,shift_param);
					points := data_fix[1];
					Coords := data_fix[2];
				else
					while not IsInt(e) do
						e := e[1];
					od;
					if e in Edges(surf) then
						if inner[e] or is_circle then
							data_fix:=_FixVertOfInnerNMEdge(surf,e,Coords,points,data,shift_param,not_split);
							Coords:=data_fix[1];
							points:=data_fix[2];
							not_split := [];
						elif Length(path) = 1 then
							# this might have unwanted behavior
							verts_current := VerticesOfEdge(surf,e);
							
							same := _ChooseIsolatedVertex(surf,e,path,points,data);
							data_fix := _FixVertOfOuterNMEdge(surf,e,Coords,points,data,shift_param,not_split);
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
							data_fix := _FixVertOfOuterNMEdge(surf,e,Coords,points,data,shift_param,not_split);
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
end);;


InstallGlobalFunction(_FixNMIntersection,function(surf,order,info,data,points,Coords,shift_param)
	local l, int_v, branches, not_split, data_fix, b;
	int_v := info[1];
	branches := info[2];
	not_split := [];
	for b in branches do
		data_fix := _FixVertOfOuterNMEdge(surf,b[1],Coords,points,data,shift_param,not_split);
		Coords:=data_fix[1];
		points:=data_fix[2];
		if Length(b[2]) > 1 then
			# the path continues in this direction
			data_fix := _FixNMPathRec(surf,[b[2],order[2],order[3]],data,points,Coords,shift_param);
			points := data_fix[1];
			Coords := data_fix[2];
		fi;
	od;
	return [points,Coords];
end);;




InstallGlobalFunction(FixNMEdgePath, function(surf,data,points,Coords,shift_param)
	local l, comp, path, order_data, not_split, data_fix, inner, is_circle, NM_verts, points_fix, s_data, verts_current, verts_comp, same, e, s,t; 
	if RamifiedEdges(surf) <> [] then
		order_data := OrderNMEdges(surf,ShallowCopy(data));
		
		not_split := [];
		is_circle := order_data[2];
		inner := order_data[3];
		
		for comp in order_data[1] do
			path := comp;
			data_fix := _FixNMPathRec(surf,[comp,is_circle,inner],data,points,Coords,shift_param);
			points := data_fix[1];
			Coords:= data_fix[2];
			
		od;
		
		s_data := _SimplicialSurfaceFromChangedCoordinates([Coords,surf],1./10^6);
		surf := s_data[1];
		points_fix := points;
	else
		points_fix := points;
		order_data := [false];
	fi;
	
    	return [surf,points_fix,Coords,order_data];
end);;



InstallGlobalFunction(FixNMVerts, function(surf, data, points, Coords, shift_param)
	local NM_verts, v, new_data, s_data;
	
	# to be executed after fixing the non-manifold edges
	
	NM_verts := ShallowCopy(RamifiedVertices(surf));
	
	for v in NM_verts do
		new_data := _FixNMVert(surf,v,Coords,points,data,shift_param);
		Coords := new_data[1];
		points := new_data[2];
		
	od;
	
	s_data := _SimplicialSurfaceFromChangedCoordinates([Coords,surf],1./10^6);
	surf := s_data[1];
	
	Coords := _UpdateNormals(Coords);
	return [surf,points,Coords];
end);;


InstallGlobalFunction(RemedyNonManifold, function(data,points, shift_param)
	local surf, normals, Coords, order_data, m_data, m_surf, m_points, m_coords, fully_m_data;
	#
	# input structure is data=[surf1,surf2,faces,normals_coordinates] (output of outer hull function), points: coordinates of the points of surf2, shift_param: norm of the vectors that will be added to shift vertices
	# surf1 is a intersection free complex 
	# surf2 its outer hull
	# faces are the faces of surf2
	# normals_coordinates are the coordinates of the normals of the faces of surf2, indexed by the resp. face number
	#
	# output has structure [fully_m_surf, fully_m_points, fully_m_coords, order_data] where fully_m_surf is the manifold version of surf2 and the other two its data. order_data is the information about the non-manifold edge paths present in the original complex
	
	
	
	surf := ShallowCopy(data[2]);
	normals := ShallowCopy(data[4]);
	
	Coords := _ConvertDataFormatPtC(surf,points,normals);
	# investigate path based
	m_data := FixNMEdgePath(surf, ShallowCopy(data),points,Coords,shift_param);
	m_surf := m_data[1];
	m_points := m_data[2];
	m_coords := m_data[3];
	order_data := m_data[4];
	fully_m_data := FixNMVerts(m_surf,data,ShallowCopy(m_points),m_coords,shift_param);
	fully_m_data[4] := order_data;
	

	return fully_m_data;
end);;

