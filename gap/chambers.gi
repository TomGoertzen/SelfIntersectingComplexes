# Given an edge e, from a triangular complex s, with vertexCoordinates vC and a face f with MyNormal vector n
# calculate the fan [(f,n),(f2,n2),...] in direction of MyNormal
BindGlobal("CalculateFan",function(s,e,vC,f,MyNormal)
	local FacesOfEdge_e, ThirdPoint, VoE, a, p, n, t, vec, face; 
	# for all faces with edge e in s arrange them first
	FacesOfEdge_e:=FacesOfEdge(s,e);
	ThirdPoint:=[];
	VoE:=VerticesOfEdge(s,e);
	for face in FacesOfEdge_e do
		Add(ThirdPoint,[face,DifferenceLists(VerticesOfFace(s,face),VoE)[1]]);
	od;
	a:=vC[VoE[1]];
	n:=vC[VoE[2]]-vC[VoE[1]];
	n:=n/MyNorm(n);
	for t in ThirdPoint do
		p:=vC[t[2]];
		t[3]:=p-a-Dot((p-a),n)*n;
		# if we consider face f save the vector
		if t[1]=f then
			vec:=t[3];
		fi;
	od;
	# now compute the angle of all derived vectors with vec and rotate MyNormal vector n according to the angle
	for t in ThirdPoint do
		if Determinant(1.*[vec,n,MyNormal]) < 0. then
			t[4]:=VectorAnglePlane(vec,t[3],n);
			t[5]:=Rotate_Vector_Axis(MyNormal,n,t[4]);
		else
			t[4]:=VectorAnglePlane(vec,t[3],-n);
			t[5]:=Rotate_Vector_Axis(MyNormal,-n,t[4]);
		fi;
	od;
	SortBy(ThirdPoint,x->x[4]);
	return ThirdPoint;
end);

BindGlobal("UpwardContinuation",function(s,e,vC,f,MyNormal)
	local Fan, index, i;
	Fan:=CalculateFan(s,e,vC,f,MyNormal);
	for i in [1..Size(Fan)] do 
		if Fan[i][1]=f then
			index:=(i mod (Size(Fan))) + 1;
			return [Fan[index][1],-Fan[index][5]];
		fi;
	od;
	return fail;
end);



# Find outer triangle
# for a triangular complex with points find an outer face and
# a normal vector pointing outwards
BindGlobal("FindOuterTriangle",function(t,points)
    local i_max,px_max,i,eov,j_min,x_max,j_max,x_min,e,f,n,foe,j,index,indices;
    # find point with maximal x coordinate
    indices:=Filtered([1..Size(points)],i->IsBound(points[i]));
    i_max:=indices[1];
    px_max:=points[i_max][1];
    for i in [2..Size(indices)] do
        index:=indices[i];
        if points[index][1]>px_max then
            px_max:=points[index][1];
            i_max:=index;
        fi;
    od;
    eov:=EdgesOfVertex(t,i_max);
    eov:=List(eov,e->[e,DifferenceLists(VerticesOfEdge(t,e),[i_max])]);
    eov:=List(eov,e->[e[1],points[e[2][1]]-points[i_max]]);
    eov:=List(eov,e->[e[1],e[2]/MyNorm(e[2])]);
    # find edge with smallest absolute value in x-coordinate
    j_min:=1;
    x_min:=AbsoluteValue(eov[1][2][1]);
    for j in [2..Size(eov)] do
        if AbsoluteValue(eov[j][2][1]) < x_min then
            j_min:=j;
            x_min:=AbsoluteValue(eov[j][2][1]);
        fi;
    od;
    # found the edge
    e:=eov[j_min][1];
    foe:=FacesOfEdge(t,e);
    foe:=List(foe,f->[f,VerticesOfFace(t,f)]);
    foe:=List(foe,f->[f[1],Crossproduct(points[f[2][2]]-points[f[2][1]],points[f[2][3]]-points[f[2][1]])]);
    foe:=List(foe,f->[f[1],f[2]/MyNorm(f[2])]);
    # find edge with smallest absolute value in x-coordinate
    j_max:=1;
    x_max:=AbsoluteValue(foe[1][2][1]);
    for j in [2..Size(foe)] do
        if AbsoluteValue(foe[j][2][1]) > x_max then
            j_max:=j;
            x_max:=AbsoluteValue(foe[j][2][1]);
        fi;
    od;
    # found the face
    f:=foe[j_max][1];
    n:=foe[j_max][2];
    if foe[j_max][2][1]<0. then
        n:=-n;
    fi;
    return [f,n];
end);

# Algorithm from Paper by M. Attene Title: As-exact-as-possible repair of unprintable stl files
# Input: a triangular complex s with coordinates vC and an outer face f with normal vector n
# Output: s, Restricted Complex to outer faces, Outer Faces, and correctly oriented normal vectors
InstallGlobalFunction(OuterHull,function(s,vC,f,n)
	local B, e, OuterTriangles, InnerTriangles, NormalVectors, b, edge, t, t_new, eps;
	eps := 1.0/(10^(12));
	B:=[];
	for e in EdgesOfFace(s,f) do
		Add(B,[f,e]);
	od;
	OuterTriangles:=[f];
	NormalVectors:=[];
	NormalVectors[f]:=n;
	InnerTriangles:=[];
	while not IsEmpty(B) do 
		t:=Remove(B);
		t_new:=UpwardContinuation(s,t[2],vC,t[1],NormalVectors[t[1]]);
		if not t_new[1] in OuterTriangles then

			Add(OuterTriangles,t_new[1]);
			NormalVectors[t_new[1]]:=MyRoundVector(t_new[2], eps);
			for edge in EdgesOfFace(s,t_new[1]) do
				if edge <> t[2] then
					Add(B,[t_new[1],edge]);
				fi;
			od;
		fi;
	od;
	return [s,SubcomplexByFaces(s,OuterTriangles),OuterTriangles,NormalVectors];
end);

InstallGlobalFunction(DrawSTLScratch,function(t,fileName,vC)
	local f,normals,normal,x,y,ccoords,vof;
	normals:=[];
	if IsOrientable(t) then
		vof:=VerticesOfFaces(t);
	else
		vof:=List(Orientation(t), f->VerticesAsList(f){[1,2,3]});
	fi;
	for f in Faces(t) do
		ccoords:=vC{vof[f]};
		x := ccoords[2]-ccoords[1];
        y := ccoords[3]-ccoords[1];
        normal := Crossproduct(x,y);
        normal := normal / Sqrt(normal*normal);
        normals[f] := normal;
    od;
    DrawSTLwithNormals(t,fileName,vC,normals,[]);
end);

#######################################################################################################################################################################
#Inputs
##
# This method takes a string and a list l in the coordinate format (l=[face1,face2,,face3,....]). The faces are also lists in the format 
# face1 = [vertex1,vertex2,vertex3,normal,vertexNumbs]. The vertices and normal then are lists as well: vertex1 = [x,y,z] where x,y,z are floats
# and vertexNumbs = [v1,v2,v3] where v1,v2,v3 are integers corresponding to the vertex indices in the simplicial face. Since those change quite a lot when fixing self-
# intersections, vertexNumbs is !depricated!
#
##
#   METHOD
##
# The method iself saves a STL file corresponding to the object with name specified by fileName
#
#######################################################################################################################################################################
InstallGlobalFunction(DrawSTLwithNormals,function(s,fileName, vC,normals,visualize_normal_list)
        local file, filesepr, name, output, x,y, i, j, k, r,l, coords, Copy_Coords, normal, eps, VoF, f,
        		middle, v2, v3, new_coords, new_normal, perms, perm,n;

        

        eps := 1.0/(10^(12));

        filesepr := SplitString(fileName, ".");
        name := filesepr[1];
        # test file name
        file := Filename( DirectoryCurrent(), Concatenation(name,".stl") );
        output := OutputTextFile( file, false ); # override other files
            
        if output = fail then
            Error(Concatenation("File ", String(file), " can't be opened.") );
        fi;        

        AppendTo(output, Concatenation("solid ", name, "\n"));
        
        for f in Faces(s) do 
            
                VoF:=VerticesOfFace(s,f);
                # get coords of vertices
                coords := [vC[VoF[1]],vC[VoF[2]],vC[VoF[3]]];
                normal := normals[f];
                
                # write normal

                AppendTo(output, "\tfacet normal ");
                for j in [1..3] do
                    AppendTo(output, Concatenation(String(normal[j])," "));
                od;
                AppendTo(output, "\n");
                    
                # write vertex coords according to the right hand rule
                perms:=[[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]];
                for perm in perms do
                	n:=Crossproduct(coords[perm[2]]-coords[perm[1]],coords[perm[3]]-coords[perm[1]]);
                	n:=-n/MyNorm(n);
                	if MyNorm(n-normal)<eps then

                		break;
                	fi;
                od;

                AppendTo(output, "\t\touter loop\n");
                for j in [1..3] do
                    AppendTo(output,"\t\t\tvertex ");
                       
                    for k in [1..3] do
                        AppendTo(output, Concatenation(String(coords[perm[j]][k])," "));
                    od;
                    AppendTo(output,"\n");
                od;
                AppendTo(output, "\t\tendloop\n");
                AppendTo(output,"\tendfacet\n");
                
                # visualize normal of face
                if f in visualize_normal_list then
                middle := 1./3 * coords[1] + 1./3 * coords[2] + 1./3 * coords[3];
                v2 := middle + 0.02 * (coords[1]-coords[2])/MyNorm(coords[1]-coords[2]);
                v3 := middle + 0.1 * normal;
                
                new_coords := [middle,v2,v3];
                new_normal := Crossproduct(middle-v2, middle - v3);
                new_normal := new_normal / MyNorm(new_normal);
                
                AppendTo(output, "\tfacet normal ");
                
                for j in [1..3] do
                    AppendTo(output, Concatenation(String(new_normal[j])," "));
                od;
                AppendTo(output, "\n");
                
                 # write vertex coords
                AppendTo(output, "\t\touter loop\n");
                for j in [1..3] do
                    AppendTo(output,"\t\t\tvertex ");
                       
                    for k in [1..3] do
                        AppendTo(output, Concatenation(String(new_coords[j][k])," "));
                    od;
                    AppendTo(output,"\n");
                od;
                AppendTo(output, "\t\tendloop\n");
                AppendTo(output,"\tendfacet\n");
                fi;
        od;
        AppendTo(output, Concatenation("endsolid ", name));
        Print("\n Saved file");
        CloseStream(output);
        return;
end);


InstallGlobalFunction(DrawSurfaceToObj,function(s,fileName, vC)
        local file, filesepr, name, output, x,y, i, j, k, r,l, coords, Copy_Coords, normal, eps, VoF, f,
        		middle, v2, v3, new_coords, new_normal,v;

        
        

        filesepr := SplitString(fileName, ".");
        name := filesepr[1];
        # test file name
        file := Filename( DirectoryCurrent(), Concatenation(name,".obj") );
        output := OutputTextFile( file, false ); # override other files
            
        if output = fail then
            Error(Concatenation("File ", String(file), " can't be opened.") );
        fi;        

        #AppendTo(output, Concatenation("solid ", name, "/n"));

        for v in vC do
        	AppendTo(output,"v ");
        	for i in [1..3] do
        		AppendTo(output,v[i]);
        		if i<3 then
        			AppendTo(output," ");
        		fi;
        	od;
        	AppendTo(output,"\n");
        od;
        
        for f in VerticesOfFaces(s) do 
        	AppendTo(output,"f ");
        	for i in [1..3] do
        		AppendTo(output,f[i]);
        		if i<3 then
        			AppendTo(output," ");
        		fi;
        	od;
        	AppendTo(output,"\n");            
        od;
        Print("\n Saved file");
        CloseStream(output);
        return;
end);




InstallGlobalFunction(DrawSTLScratchVOF,function(vof,fileName,vC)
	local f,normals,normal,x,y,ccoords,i;
	normals:=[];
	for i in [1..Size(vof)] do
		ccoords:=vC{vof[i]};
		x := ccoords[2]-ccoords[1];
        y := ccoords[3]-ccoords[1];
        normal := Crossproduct(x,y);
        normal := normal / Sqrt(normal*normal);
        normals[i] := normal;
    od;
    DrawSTLwithNormalsVOF(vof,fileName,vC,normals,[]);
end);


InstallGlobalFunction(DrawSTLwithNormalsVOF,function(vof,fileName, vC,normals,visualize_normal_list)
        local file, filesepr, name, output, x,y, i, j, k, r,l, coords, Copy_Coords, normal, eps, VoF, f,
        		middle, v2, v3, new_coords, new_normal, perms, perm,n;

        eps := 1.0/(10^(12));

        filesepr := SplitString(fileName, ".");
        name := filesepr[1];
        # test file name
        file := Filename( DirectoryCurrent(), Concatenation(name,".stl") );
        output := OutputTextFile( file, false ); # override other files
            
        if output = fail then
            Error(Concatenation("File ", String(file), " can't be opened.") );
        fi;        

        AppendTo(output, Concatenation("solid ", name, "\n"));
        
        for f in [1..Size(vof)] do 
            
                VoF:=vof[f];
                # get coords of vertices
                coords := [vC[VoF[1]],vC[VoF[2]],vC[VoF[3]]];
                normal := normals[f];
                
                # write normal

                AppendTo(output, "\tfacet normal ");
                for j in [1..3] do
                    AppendTo(output, Concatenation(String(normal[j])," "));
                od;
                AppendTo(output, "\n");
                    
                # write vertex coords according to the right hand rule
                perms:=[[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]];
                for perm in perms do
                	n:=Crossproduct(coords[perm[2]]-coords[perm[1]],coords[perm[3]]-coords[perm[1]]);
                	n:=-n/MyNorm(n);
                	if MyNorm(n-normal)<eps then

                		break;
                	fi;
                od;

                AppendTo(output, "\t\touter loop\n");
                for j in [1..3] do
                    AppendTo(output,"\t\t\tvertex ");
                       
                    for k in [1..3] do
                        AppendTo(output, Concatenation(String(coords[perm[j]][k])," "));
                    od;
                    AppendTo(output,"\n");
                od;
                AppendTo(output, "\t\tendloop\n");
                AppendTo(output,"\tendfacet\n");
                
                # visualize normal of face
                if f in visualize_normal_list then
                middle := 1./3 * coords[1] + 1./3 * coords[2] + 1./3 * coords[3];
                v2 := middle + 0.02 * (coords[1]-coords[2])/MyNorm(coords[1]-coords[2]);
                v3 := middle + 0.1 * normal;
                
                new_coords := [middle,v2,v3];
                new_normal := Crossproduct(middle-v2, middle - v3);
                new_normal := new_normal / MyNorm(new_normal);
                
                AppendTo(output, "\tfacet normal ");
                
                for j in [1..3] do
                    AppendTo(output, Concatenation(String(new_normal[j])," "));
                od;
                AppendTo(output, "\n");
                
                 # write vertex coords
                AppendTo(output, "\t\touter loop\n");
                for j in [1..3] do
                    AppendTo(output,"\t\t\tvertex ");
                       
                    for k in [1..3] do
                        AppendTo(output, Concatenation(String(new_coords[j][k])," "));
                    od;
                    AppendTo(output,"\n");
                od;
                AppendTo(output, "\t\tendloop\n");
                AppendTo(output,"\tendfacet\n");
                fi;
        od;
        AppendTo(output, Concatenation("endsolid ", name));
        Print("\n Saved file");
        CloseStream(output);
        return;
end);

BindGlobal("NormalsOuterHull",function(t,points)
    local data,f,n;
    data:=FindOuterTriangle(t,Sublist(points,Vertices(t)));
    f:=data[1];
    n:=data[2];
    return OuterHull(t,points,f,n)[4];
end);

# computing all the chambers of the given triangular complex t
# with vertex embedding in points
# returns two lists
# first list contains components and second list normal vectors
InstallGlobalFunction(ComputeChambers,function(t,points)
    local components,current_face,data,vof,n,c,count,faces;
    components:=[[],[]];
    faces:=ShallowCopy(Faces(t));
    for current_face in faces do
        count:=0;
        for c in components[1] do
            if current_face in Faces(c) then
                count:=count+1;
            fi;
        od;
        if not count=2 then
            vof:=VerticesOfFace(t,current_face);
            n:=Crossproduct(points[vof[2]]-points[vof[1]],points[vof[3]]-points[vof[1]]);
            n:=n/MyNorm(n);
            data:=OuterHull(t,points,current_face,n);
            if not data[2] in components[1] then
                Add(components[1],data[2]);
                Add(components[2],data[4]);
                faces:=Difference(faces,Faces(data[2]));
            fi;
            data:=OuterHull(t,points,current_face,-n);
            if not data[2] in components[1] then
                Add(components[1],data[2]);
                Add(components[2],data[4]);
                faces:=Difference(faces,Faces(data[2]));
            fi;
        fi;
    od;
    return components;
end);


BindGlobal("ChambersMergeSurfacesCoordinates",function(data)
    local vof, coordinates,i;
    vof:=[];
    coordinates:=[];
    for i in [1..Size(data)] do
        vof:=Concatenation(vof,VerticesOfFaces(data[i][1])+Size(coordinates));
        coordinates:=Concatenation(coordinates,data[i][2]);
    od;
    return [TriangularComplexByVerticesInFaces(vof),coordinates];
end);

#This function requires the GAPic Package
#Use the output of the PrintableOuterHull function and a name for the resulting html file
InstallGlobalFunction(AnimateChamberExplode,function(data_tri,name)
    local comps,InvertList,FixCoordinateFormat,c,data_big,s_assembly,coordinates_assembly,y,explode_vector_trans,points,new_points,i,j,v,verticesPositions,printRecord,params;
    comps:=ComputeChambers(data_tri[1],data_tri[3]);;
    InvertList:=function(list)
      local new_list,i;
      new_list:=[];
      for i in [1..Size(list)] do
        new_list[list[i]]:=i;
      od;
      return new_list;
    end;;
        
    FixCoordinateFormat:=function(s,coordinates)
      local map,new_faces, i, j, new_coordinates;
      map:=InvertList(Vertices(s));
      new_faces:=List(VerticesOfFaces(s),f->List(f,j->map[j]));
      new_coordinates:=List([1..Size(Vertices(s))],i->coordinates[Vertices(s)[i]]);
      return [TriangularComplexByVerticesInFaces(new_faces),new_coordinates];
    end;;
    comps:=List([1..Size(comps[1])],i->FixCoordinateFormat(comps[1][i],data_tri[3]));;
    for c in comps do
      if Size(Vertices(c[1]))=Size(Vertices(data_tri[2])) then
        Remove(comps,Position(comps,c));
      fi;
    od;
    data_big:=ChambersMergeSurfacesCoordinates(comps);;
    s_assembly:=data_big[1];;
    coordinates_assembly:=data_big[2];;
    LoadPackage("GapIC");
    explode_vector_trans:=function(s,points,strength)
      local center,new_points,c,v,center_comp,vector;
      center:=Sum(points)/Size(points);
      new_points:=[];
      vector:=[];
      for c in ConnectedComponents(s) do
        center_comp:=Sum(points{Vertices(c)})/NumberOfVertices(c);
        vector:=Concatenation(vector,center_comp-center);
        
      od;
      return vector;
    end;;

    y:=explode_vector_trans(s_assembly,coordinates_assembly,0.1);;

    points:=coordinates_assembly;;
    new_points:=[];;
    for i in [1..Size(ConnectedComponents(s_assembly))] do
      c:=ConnectedComponents(s_assembly)[i];
      for v in Vertices(c) do
        new_points[v]:=[];
        for j in [1..3] do
          new_points[v][j]:=Concatenation(String(points[v][j]),"+a*",String(y[j+(i-1)*3]));
        od;
      od;
    od;

    verticesPositions:=new_points;;
    printRecord := SetVertexCoordinatesParameterized(s_assembly, verticesPositions, rec());;
    params := [["a", 0, [0,10]]];;
    printRecord := SetVertexParameters(s_assembly, params, printRecord);;
    #printRecord.faceColours:=Concatenation(List([1..Size(assembly[1][1])],i->"0x006ab3"),List([1..Size(assembly[2][1])],i->"0xec7d05"));
    DrawComplexToJavaScript(s_assembly, name, printRecord);;
    return [s_assembly,coordinates_assembly,printRecord];
end);