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