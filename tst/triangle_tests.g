eps:=1.*10^-6;;
faces:=[[1,2,3],[4,5,6]];;
s:=SimplicialSurfaceByVerticesInFaces(faces);;
list:=[
[[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0.5, 2], [0, 0, 0], [1.5, 1.5, 1]], 
[[2, 2, 0], [2, 0, 2], [0, 2, 2], [1, 1, 0], [3, 1, 2], [1, 3, 2]], 
[[-1, 0, 0], [1, 0, 0], [0, 1, 1], [0, 0, 0], [1, 1, 2], [-1, -1, 2]], 
[[0, 0, 0], [2, 0, 0], [1, 2, 0], [1, 1, -1], [2, 1, -1], [1, 0, 1]],
  [[0, 0, 0], [1, 1, 0], [0, 1, 1], [1, 0, 1], [2, 2, 0], [-2, 0, 2]], 
  [[1, 0, 0], [0, 1, 0], [0, 0, 2], [-2, -1, 1], [-1, -2, 2], [2, 2, 1]], 
  [[0, 1, 2], [2, 0, 1], [1, 2, 0], [1, 1, 1], [0, 2, 1], [2, 1, 0]], 
  [[0, 1, 2], [2, 0, 1], [1, 2, 0], [-1, -1, 1], [0, 2, 1], [2, 1, 0]], 
  [[2, 0, 0], [0, 2, 0], [0, 0, 0], [2, 2, 0], [1, -1, 0], [-1, 1, 0]],
  [[2, 0, 0], [0, 2, 0], [-1, -1, 0], [2, 2, 0], [1, -1, 0], [-1, 1, 0]], 
  [[1, 0, 0], [0, 1, 0], [0, 0, 0], [0.5, 0.5, 0], [0.5, -0.5, 0], [1.5, 1.5, 0]],
   [[2, 2, 0], [2, 0, 0], [0, 2, 0], [1, 1, 0], [3, 1, 0], [1, 3, 0]],
    [[-1, 0, 0], [1, 0, 0], [0, 1, 0], [0, -1, 0], [1, 1, 0], [-1, -1, 0]],
     [[0, 0, 0], [2, 0, 0], [1, 2, 0], [1, 1, 0], [2, 1, 0], [1, 0, 0]],
      [[0, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0], [0.5, 0.5, 0], [2, 0, 0]],
      [[0, 0, 0], [2, 0, 0], [0, 2, 0], [0.5, 0.5, 0], [1.5, 0.5, 0], [0.5, 1.5, 0]],
      [[0, 0, 0], [2, 0, 0], [0, 2, 0], [0.5, 0.5, 0], [1., 0.5, 0], [0.5, 1., 0]],
      [[0, 0, 0], [2, 0, 0], [0, 2, 0], [-0.5, -0.5, 0], [1., 0.5, 0], [0.5, 1., 0]],
      [[0, 0, 0], [2, 0, 0], [0, 2, 0], [2.5, 2.5, 0], [1., 0.5, 0], [0.5, 1., 0]],
      [[0, 0, 0], [2, 0, 0], [1, 2, 0], [1, 1.5, 0], [0, 3, 0], [2, 3, 0]], 
      [[0, 0, 0], [2, 0, 0], [1, 2, 0], [-1, 0, 0], [1, 0, 0], [-2, -3, 0]]];;

data:=PrintableOuterHull(s,list[1]*1.,Concatenation("tst/Triangles/Triangle",String(1)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[2]*1.,Concatenation("tst/Triangles/Triangle",String(2)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[3]*1.,Concatenation("tst/Triangles/Triangle",String(3)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[4]*1.,Concatenation("tst/Triangles/Triangle",String(4)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[5]*1.,Concatenation("tst/Triangles/Triangle",String(5)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[6]*1.,Concatenation("tst/Triangles/Triangle",String(6)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[7]*1.,Concatenation("tst/Triangles/Triangle",String(7)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[8]*1.,Concatenation("tst/Triangles/Triangle",String(8)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[9]*1.,Concatenation("tst/Triangles/Triangle",String(9)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[10]*1.,Concatenation("tst/Triangles/Triangle",String(10)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[11]*1.,Concatenation("tst/Triangles/Triangle",String(11)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[12]*1.,Concatenation("tst/Triangles/Triangle",String(12)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[13]*1.,Concatenation("tst/Triangles/Triangle",String(13)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[14]*1.,Concatenation("tst/Triangles/Triangle",String(14)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[15]*1.,Concatenation("tst/Triangles/Triangle",String(15)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[16]*1.,Concatenation("tst/Triangles/Triangle",String(16)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[17]*1.,Concatenation("tst/Triangles/Triangle",String(17)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[18]*1.,Concatenation("tst/Triangles/Triangle",String(18)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[19]*1.,Concatenation("tst/Triangles/Triangle",String(19)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[20]*1.,Concatenation("tst/Triangles/Triangle",String(20)),eps,0.1,Group(()),true,false);
data:=PrintableOuterHull(s,list[21]*1.,Concatenation("tst/Triangles/Triangle",String(21)),eps,0.1,Group(()),true,false);
