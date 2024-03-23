gap> eps:=1.*10^-6;;
gap> faces:=[[1,2,3],[4,5,6]];;
gap> s:=SimplicialSurfaceByVerticesInFaces(faces);;
gap> list:=[
> [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0.5, 2], [0, 0, 0], [1.5, 1.5, 1]], 
> [[2, 2, 0], [2, 0, 2], [0, 2, 2], [1, 1, 0], [3, 1, 2], [1, 3, 2]], 
> [[-1, 0, 0], [1, 0, 0], [0, 1, 1], [0, 0, 0], [1, 1, 2], [-1, -1, 2]], 
> [[0, 0, 0], [2, 0, 0], [1, 2, 0], [1, 1, -1], [2, 1, -1], [1, 0, 1]],
>   [[0, 0, 0], [1, 1, 0], [0, 1, 1], [1, 0, 1], [2, 2, 0], [-2, 0, 2]], 
>   [[1, 0, 0], [0, 1, 0], [0, 0, 2], [-2, -1, 1], [-1, -2, 2], [2, 2, 1]], 
>   [[0, 1, 2], [2, 0, 1], [1, 2, 0], [1, 1, 1], [0, 2, 1], [2, 1, 0]], 
>   [[0, 1, 2], [2, 0, 1], [1, 2, 0], [-1, -1, 1], [0, 2, 1], [2, 1, 0]], 
>   [[2, 0, 0], [0, 2, 0], [0, 0, 0], [2, 2, 0], [1, -1, 0], [-1, 1, 0]],
>   [[2, 0, 0], [0, 2, 0], [-1, -1, 0], [2, 2, 0], [1, -1, 0], [-1, 1, 0]], 
>   [[1, 0, 0], [0, 1, 0], [0, 0, 0], [0.5, 0.5, 0], [0.5, -0.5, 0], [1.5, 1.5, 0]],
>    [[2, 2, 0], [2, 0, 0], [0, 2, 0], [1, 1, 0], [3, 1, 0], [1, 3, 0]],
>     [[-1, 0, 0], [1, 0, 0], [0, 1, 0], [0, -1, 0], [1, 1, 0], [-1, -1, 0]],
>      [[0, 0, 0], [2, 0, 0], [1, 2, 0], [1, 1, 0], [2, 1, 0], [1, 0, 0]],
>       [[0, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0], [0.5, 0.5, 0], [2, 0, 0]],
>       [[0, 0, 0], [2, 0, 0], [0, 2, 0], [0.5, 0.5, 0], [1.5, 0.5, 0], [0.5, 1.5, 0]],
>       [[0, 0, 0], [2, 0, 0], [0, 2, 0], [0.5, 0.5, 0], [1., 0.5, 0], [0.5, 1., 0]],
>       [[0, 0, 0], [2, 0, 0], [0, 2, 0], [-0.5, -0.5, 0], [1., 0.5, 0], [0.5, 1., 0]],
>       [[0, 0, 0], [2, 0, 0], [0, 2, 0], [2.5, 2.5, 0], [1., 0.5, 0], [0.5, 1., 0]],
>       [[0, 0, 0], [2, 0, 0], [1, 2, 0], [1, 1.5, 0], [0, 3, 0], [2, 3, 0]], 
>       [[0, 0, 0], [2, 0, 0], [1, 2, 0], [-1, 0, 0], [1, 0, 0], [-2, -3, 0]]];;
gap> 
gap> data:=PrintableOuterHull(s,list[1]*1.,Concatenation("tst/Triangles/Triangle",String(1)),eps,0.1,Group(()),true,false);

 Saved file[ triangular complex (8 vertices, 15 edges, and 8 faces)
    , triangular complex (7 vertices, 13 edges, and 7 faces)
    , 
  [ [ 1., 0., 0. ], [ 0., 1., 0. ], [ 0., 0., 1. ], 
      [ 0.166667, 0.166667, 0.666667 ], [ 0.375, 0.375, 0.25 ], 
      [ 0.5, 0.5, 2. ], [ 0., 0., 0. ], [ 1.5, 1.5, 1. ] ], 
  [ [ 0.57735, 0.57735, 0.57735 ], [ 0.57735, 0.57735, 0.57735 ], 
      [ 0.57735, 0.57735, 0.57735 ], [ 0.57735, 0.57735, 0.57735 ], 
      [ 0.57735, 0.57735, 0.57735 ], [ 0.707107, -0.707107, 0. ],, 
      [ 0.707107, -0.707107, -0. ] ] ]
gap> data:=PrintableOuterHull(s,list[2]*1.,Concatenation("tst/Triangles/Triangle",String(2)),eps,0.1,Group(()),true,false);

 Saved file[ triangular complex (8 vertices, 13 edges, and 6 faces)
    , simplicial surface (5 vertices, 7 edges, and 3 faces)
    , [ [ 2., 2., 0. ], [ 2., 0., 2. ], [ 0., 2., 2. ], [ 2., 1., 1. ], 
      [ 1., 2., 1. ], [ 1., 1., 0. ], [ 3., 1., 2. ], [ 1., 3., 2. ] ], 
  [ [ 0.57735, 0.57735, 0.57735 ],,,, [ 0.57735, 0.57735, -0.57735 ], 
      [ 0.57735, 0.57735, -0.57735 ] ] ]
gap> data:=PrintableOuterHull(s,list[3]*1.,Concatenation("tst/Triangles/Triangle",String(3)),eps,0.1,Group(()),true,false);

 Saved file[ triangular complex (6 vertices, 8 edges, and 3 faces)
    , simplicial surface (4 vertices, 5 edges, and 2 faces)
    , 
  [ [ -1., 0., 0. ], [ 1., 0., 0. ], [ 0., 1., 1. ], [ 0., 0., 0. ], 
      [ 1., 1., 2. ], [ -1., -1., 2. ] ], 
  [ [ 0., -0.707107, 0.707107 ], [ 0., -0.707107, 0.707107 ] ] ]
gap> data:=PrintableOuterHull(s,list[4]*1.,Concatenation("tst/Triangles/Triangle",String(4)),eps,0.1,Group(()),true,false);

 Saved file[ triangular complex (8 vertices, 15 edges, and 8 faces)
    , triangular complex (6 vertices, 11 edges, and 6 faces)
    , 
  [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 1., 2., 0. ], [ 1.5, 0.5, 0. ], 
      [ 1., 0.5, 0. ], [ 1., 1., -1. ], [ 2., 1., -1. ], [ 1., 0., 1. ] 
     ], 
  [ [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], 
      [ 0., 0., 1. ],, [ 0., 0.894427, 0.447214 ] ] ]
gap> data:=PrintableOuterHull(s,list[5]*1.,Concatenation("tst/Triangles/Triangle",String(5)),eps,0.1,Group(()),true,false);

 Saved file[ triangular complex (6 vertices, 8 edges, and 3 faces)
    , simplicial surface (4 vertices, 5 edges, and 2 faces)
    , [ [ 0., 0., 0. ], [ 1., 1., 0. ], [ 0., 1., 1. ], [ 1., 0., 1. ], 
      [ 2., 2., 0. ], [ -2., 0., 2. ] ], 
  [ , [ 0.301511, 0.301511, 0.904534 ], 
      [ 0.301511, 0.301511, 0.904534 ] ] ]
gap> data:=PrintableOuterHull(s,list[6]*1.,Concatenation("tst/Triangles/Triangle",String(6)),eps,0.1,Group(()),true,false);

 Saved file[ triangular complex (8 vertices, 14 edges, and 7 faces)
    , triangular complex (7 vertices, 12 edges, and 6 faces)
    , 
  [ [ 1., 0., 0. ], [ 0., 1., 0. ], [ 0., 0., 2. ], 
      [ -4.44089e-16, 0.5, 1. ], [ 0.294118, 0., 1.41176 ], 
      [ -2., -1., 1. ], [ -1., -2., 2. ], [ 2., 2., 1. ] ], 
  [ [ 0.666667, 0.666667, 0.333333 ], [ 0.666667, 0.666667, 0.333333 ],
      , [ 0.348743, -0.464991, -0.813733 ], 
      [ 0.348743, -0.464991, -0.813733 ], 
      [ 0.348743, -0.464991, -0.813733 ], 
      [ 0.348743, -0.464991, -0.813733 ] ] ]
gap> data:=PrintableOuterHull(s,list[7]*1.,Concatenation("tst/Triangles/Triangle",String(7)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (10 vertices, 18 edges, and 9 faces)
    , simplicial surface (10 vertices, 18 edges, and 9 faces)
    , 
  [ [ 0., 1., 2. ], [ 2., 0., 1. ], [ 1., 2., 0. ], [ 0.5, 1.5, 1. ], 
      [ 1., 1., 1. ], [ 1.5, 1., 0.5 ], [ 0.666667, 1.66667, 0.666667 ],
      [ 1.33333, 1.33333, 0.333333 ], [ 0., 2., 1. ], [ 2., 1., 0. ] ], 
  [ [ 0.57735, 0.57735, 0.57735 ], [ 0.57735, 0.57735, 0.57735 ], 
      [ 0.57735, 0.57735, 0.57735 ], [ 0.57735, 0.57735, 0.57735 ], 
      [ 0.57735, 0.57735, 0.57735 ], [ 0.57735, 0.57735, 0.57735 ], 
      [ 0.57735, 0.57735, 0.57735 ], [ 0.57735, 0.57735, 0.57735 ], 
      [ 0.57735, 0.57735, 0.57735 ] ] ]
gap> data:=PrintableOuterHull(s,list[8]*1.,Concatenation("tst/Triangles/Triangle",String(8)),eps,0.1,Group(()),true,false);

 Saved file[ triangular complex (8 vertices, 13 edges, and 6 faces)
    , simplicial surface (5 vertices, 7 edges, and 3 faces)
    , 
  [ [ 0., 1., 2. ], [ 2., 0., 1. ], [ 1., 2., 0. ], 
      [ 1.33333, 1.33333, 0.333333 ], [ 0.666667, 1.66667, 0.666667 ], 
      [ -1., -1., 1. ], [ 0., 2., 1. ], [ 2., 1., 0. ] ], 
  [ [ 0.57735, 0.57735, 0.57735 ], [ 0.57735, 0.57735, 0.57735 ], 
      [ 0.57735, 0.57735, 0.57735 ] ] ]
gap> data:=PrintableOuterHull(s,list[9]*1.,Concatenation("tst/Triangles/Triangle",String(9)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (10 vertices, 17 edges, and 8 faces)
    , simplicial surface (10 vertices, 17 edges, and 8 faces)
    , 
  [ [ 2., 0., 0. ], [ 0., 2., 0. ], [ 0., 0., 0. ], [ 1.5, 0.5, 0. ], 
      [ 1.33333, 0., 0. ], [ 0.5, 1.5, 0. ], [ 0., 1.33333, 0. ], 
      [ 2., 2., 0. ], [ 1., -1., 0. ], [ -1., 1., 0. ] ], 
  [ [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], 
      [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ] ] ]
gap> data:=PrintableOuterHull(s,list[10]*1.,Concatenation("tst/Triangles/Triangle",String(10)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (12 vertices, 21 edges, and 10 faces)
    , simplicial surface (12 vertices, 21 edges, and 10 faces)
    , 
  [ [ 2., 0., 0. ], [ 0., 2., 0. ], [ -1., -1., 0. ], [ 1.5, 0.5, 0. ], 
      [ 1.25, -0.25, 0. ], [ 0.5, 1.5, 0. ], [ -0.25, 1.25, 0. ], 
      [ 0.5, -0.5, 0. ], [ -0.5, 0.5, 0. ], [ 2., 2., 0. ], 
      [ 1., -1., 0. ], [ -1., 1., 0. ] ], 
  [ [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], 
      [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], 
      [ 0., 0., 1. ], [ 0., 0., 1. ] ] ]
gap> data:=PrintableOuterHull(s,list[11]*1.,Concatenation("tst/Triangles/Triangle",String(11)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (9 vertices, 15 edges, and 7 faces)
    , simplicial surface (9 vertices, 15 edges, and 7 faces)
    , 
  [ [ 1., 0., 0. ], [ 0., 1., 0. ], [ 0., 0., 0. ], [ 0.5, 0.5, 0. ], 
      [ 0.5, 0., 0. ], [ 0.833333, 0.166667, 0. ], [ 0.75, 0., 0. ], 
      [ 0.5, -0.5, 0. ], [ 1.5, 1.5, 0. ] ], 
  [ [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], 
      [ -0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ] ] ]
gap> data:=PrintableOuterHull(s,list[12]*1.,Concatenation("tst/Triangles/Triangle",String(12)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (8 vertices, 13 edges, and 6 faces)
    , simplicial surface (8 vertices, 13 edges, and 6 faces)
    , [ [ 2., 2., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 2., 1., 0. ], 
      [ 1., 1., 0. ], [ 1., 2., 0. ], [ 3., 1., 0. ], [ 1., 3., 0. ] ], 
  [ [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], 
      [ 0., 0., 1. ], [ 0., 0., 1. ] ] ]
gap> data:=PrintableOuterHull(s,list[13]*1.,Concatenation("tst/Triangles/Triangle",String(13)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (10 vertices, 17 edges, and 8 faces)
    , simplicial surface (10 vertices, 17 edges, and 8 faces)
    , 
  [ [ -1., 0., 0. ], [ 1., 0., 0. ], [ 0., 1., 0. ], [ 0.5, 0., 0. ], 
      [ 0.666667, 0.333333, 0. ], [ 0., 0., 0. ], [ 0.5, 0.5, 0. ], 
      [ 0., -1., 0. ], [ 1., 1., 0. ], [ -1., -1., 0. ] ], 
  [ [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], 
      [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], 
      [ 0., 0., -1. ] ] ]
gap> data:=PrintableOuterHull(s,list[14]*1.,Concatenation("tst/Triangles/Triangle",String(14)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (8 vertices, 14 edges, and 7 faces)
    , simplicial surface (8 vertices, 14 edges, and 7 faces)
    , [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 1., 2., 0. ], [ 1., 1., 0. ], 
      [ 1., 0., 0. ], [ 1.5, 1., 0. ], [ 1.66667, 0.666667, 0. ], 
      [ 2., 1., 0. ] ], 
  [ [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], 
      [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ] ] ]
gap> data:=PrintableOuterHull(s,list[15]*1.,Concatenation("tst/Triangles/Triangle",String(15)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (5 vertices, 7 edges, and 3 faces)
    , simplicial surface (5 vertices, 7 edges, and 3 faces)
    , 
  [ [ 0., 0., 0. ], [ 1., 1., 0. ], [ 0., 1., 0. ], [ 0.5, 0.5, 0. ], 
      [ 2., 0., 0. ] ], 
  [ [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ] ] ]
gap> data:=PrintableOuterHull(s,list[16]*1.,Concatenation("tst/Triangles/Triangle",String(16)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (6 vertices, 10 edges, and 5 faces)
    , simplicial surface (6 vertices, 10 edges, and 5 faces)
    , 
  [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 0.5, 0.5, 0. ], 
      [ 1.5, 0.5, 0. ], [ 0.5, 1.5, 0. ] ], 
  [ [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], 
      [ 0., 0., -1. ] ] ]
gap> data:=PrintableOuterHull(s,list[17]*1.,Concatenation("tst/Triangles/Triangle",String(17)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (6 vertices, 12 edges, and 7 faces)
    , simplicial surface (6 vertices, 12 edges, and 7 faces)
    , 
  [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 0.5, 0.5, 0. ], 
      [ 1., 0.5, 0. ], [ 0.5, 1., 0. ] ], 
  [ [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], 
      [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ] ] ]
gap> data:=PrintableOuterHull(s,list[18]*1.,Concatenation("tst/Triangles/Triangle",String(18)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (8 vertices, 16 edges, and 9 faces)
    , simplicial surface (8 vertices, 16 edges, and 9 faces)
    , 
  [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 1., 0.5, 0. ], 
      [ 0.5, 1., 0. ], [ 0.25, 0., 0. ], [ 0., 0.25, 0. ], 
      [ -0.5, -0.5, 0. ] ], 
  [ [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], 
      [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], 
      [ 0., 0., 1. ] ] ]
gap> data:=PrintableOuterHull(s,list[19]*1.,Concatenation("tst/Triangles/Triangle",String(19)),eps,0.1,Group(()),true,false);

 Saved file[ triangular complex (8 vertices, 16 edges, and 10 faces)
    , triangular complex (8 vertices, 16 edges, and 10 faces)
    , 
  [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 1., 0.5, 0. ], 
      [ 0.5, 1., 0. ], [ 1.21429, 0.785714, 0. ], 
      [ 0.785714, 1.21429, 0. ], [ 2.5, 2.5, 0. ] ], 
  [ [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], 
      [ 0., 0., -1. ], [ 0., 0., 1. ], [ 0., 0., -1. ], [ 0., 0., -1. ],
      [ 0., 0., 1. ], [ 0., 0., -1. ] ] ]
gap> data:=PrintableOuterHull(s,list[20]*1.,Concatenation("tst/Triangles/Triangle",String(20)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (8 vertices, 15 edges, and 8 faces)
    , simplicial surface (8 vertices, 15 edges, and 8 faces)
    , 
  [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 1., 2., 0. ], 
      [ 0.857143, 1.71429, 0. ], [ 1., 1.5, 0. ], 
      [ 1.14286, 1.71429, 0. ], [ 0., 3., 0. ], [ 2., 3., 0. ] ], 
  [ [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], 
      [ 0., 0., -1. ], [ 0., 0., -1. ], [ 0., 0., -1. ], 
      [ 0., 0., -1. ] ] ]
gap> data:=PrintableOuterHull(s,list[21]*1.,Concatenation("tst/Triangles/Triangle",String(21)),eps,0.1,Group(()),true,false);

 Saved file[ simplicial surface (6 vertices, 9 edges, and 4 faces)
    , simplicial surface (6 vertices, 9 edges, and 4 faces)
    , [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 1., 2., 0. ], [ 1., 0., 0. ], 
      [ -1., 0., 0. ], [ -2., -3., 0. ] ], 
  [ [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ], [ 0., 0., 1. ] ] ]
gap>