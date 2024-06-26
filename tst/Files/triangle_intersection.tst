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
gap> Retriangulation(s,list[1]*1.);
[ [ [ 1., 0., 0. ], [ 0., 1., 0. ], [ 0., 0., 1. ], 
      [ 0.166667, 0.166667, 0.666667 ], [ 0.375, 0.375, 0.25 ], 
      [ 0.5, 0.5, 2. ], [ 0., 0., 0. ], [ 1.5, 1.5, 1. ] ], 
  [ [ 1, 2, 5 ], [ 1, 3, 4 ], [ 1, 4, 5 ], [ 2, 3, 4 ], [ 2, 4, 5 ], 
      [ 4, 5, 6 ], [ 4, 5, 7 ], [ 5, 6, 8 ] ] ]
gap> Retriangulation(s,list[2]*1.);
[ [ [ 2., 2., 0. ], [ 2., 0., 2. ], [ 0., 2., 2. ], [ 2., 1., 1. ], 
      [ 1., 2., 1. ], [ 1., 1., 0. ], [ 3., 1., 2. ], [ 1., 3., 2. ] ], 
  [ [ 1, 4, 5 ], [ 2, 3, 4 ], [ 3, 4, 5 ], [ 4, 5, 6 ], [ 4, 5, 8 ], 
      [ 4, 7, 8 ] ] ]
gap> Retriangulation(s,list[3]*1.);
[ [ [ -1., 0., 0. ], [ 1., 0., 0. ], [ 0., 1., 1. ], [ 0., 0., 0. ], 
      [ 1., 1., 2. ], [ -1., -1., 2. ] ], 
  [ [ 1, 3, 4 ], [ 2, 3, 4 ], [ 4, 5, 6 ] ] ]
gap> Retriangulation(s,list[4]*1.);
[ [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 1., 2., 0. ], [ 1.5, 0.5, 0. ], 
      [ 1., 0.5, 0. ], [ 1., 1., -1. ], [ 2., 1., -1. ], [ 1., 0., 1. ] ], 
  [ [ 1, 2, 5 ], [ 1, 3, 5 ], [ 2, 3, 4 ], [ 2, 4, 5 ], [ 3, 4, 5 ], 
      [ 4, 5, 6 ], [ 4, 5, 8 ], [ 4, 6, 7 ] ] ]
gap> Retriangulation(s,list[5]*1.);
[ [ [ 0., 0., 0. ], [ 1., 1., 0. ], [ 0., 1., 1. ], [ 1., 0., 1. ], 
      [ 2., 2., 0. ], [ -2., 0., 2. ] ], 
  [ [ 1, 2, 3 ], [ 3, 4, 5 ], [ 3, 4, 6 ] ] ]
gap> Retriangulation(s,list[6]*1.);
[ [ [ 1., 0., 0. ], [ 0., 1., 0. ], [ 0., 0., 2. ], [ -4.44089e-16, 0.5, 1. ],
      [ 0.294118, 0., 1.41176 ], [ -2., -1., 1. ], [ -1., -2., 2. ], 
      [ 2., 2., 1. ] ], 
  [ [ 1, 2, 4 ], [ 1, 4, 5 ], [ 3, 4, 5 ], [ 4, 5, 6 ], [ 4, 5, 8 ], 
      [ 5, 6, 7 ], [ 5, 7, 8 ] ] ]
gap> Retriangulation(s,list[7]*1.);
[ [ [ 0., 1., 2. ], [ 2., 0., 1. ], [ 1., 2., 0. ], [ 0.5, 1.5, 1. ], 
      [ 1., 1., 1. ], [ 1.5, 1., 0.5 ], [ 0.666667, 1.66667, 0.666667 ], 
      [ 1.33333, 1.33333, 0.333333 ], [ 0., 2., 1. ], [ 2., 1., 0. ] ], 
  [ [ 1, 2, 5 ], [ 1, 4, 5 ], [ 2, 5, 6 ], [ 3, 7, 8 ], [ 4, 5, 7 ], 
      [ 4, 7, 9 ], [ 5, 6, 8 ], [ 5, 7, 8 ], [ 6, 8, 10 ] ] ]
gap> Retriangulation(s,list[8]*1.);
[ [ [ 0., 1., 2. ], [ 2., 0., 1. ], [ 1., 2., 0. ], 
      [ 1.33333, 1.33333, 0.333333 ], [ 0.666667, 1.66667, 0.666667 ], 
      [ -1., -1., 1. ], [ 0., 2., 1. ], [ 2., 1., 0. ] ], 
  [ [ 1, 2, 5 ], [ 2, 4, 5 ], [ 3, 4, 5 ], [ 4, 5, 6 ], [ 4, 6, 8 ], 
      [ 5, 6, 7 ] ] ]
gap> Retriangulation(s,list[9]*1.);
[ [ [ 2., 0., 0. ], [ 0., 2., 0. ], [ 0., 0., 0. ], [ 1.5, 0.5, 0. ], 
      [ 1.33333, 0., 0. ], [ 0.5, 1.5, 0. ], [ 0., 1.33333, 0. ], 
      [ 2., 2., 0. ], [ 1., -1., 0. ], [ -1., 1., 0. ] ], 
  [ [ 1, 4, 5 ], [ 2, 6, 7 ], [ 3, 4, 5 ], [ 3, 4, 6 ], [ 3, 5, 9 ], 
      [ 3, 6, 7 ], [ 3, 7, 10 ], [ 4, 6, 8 ] ] ]
gap> Retriangulation(s,list[10]*1.);
[ [ [ 2., 0., 0. ], [ 0., 2., 0. ], [ -1., -1., 0. ], [ 1.5, 0.5, 0. ], 
      [ 1.25, -0.25, 0. ], [ 0.5, 1.5, 0. ], [ -0.25, 1.25, 0. ], 
      [ 0.5, -0.5, 0. ], [ -0.5, 0.5, 0. ], [ 2., 2., 0. ], [ 1., -1., 0. ], 
      [ -1., 1., 0. ] ], 
  [ [ 1, 4, 5 ], [ 2, 6, 7 ], [ 3, 8, 9 ], [ 4, 5, 8 ], [ 4, 6, 8 ], 
      [ 4, 6, 10 ], [ 5, 8, 11 ], [ 6, 7, 9 ], [ 6, 8, 9 ], [ 7, 9, 12 ] ] ]
gap> Retriangulation(s,list[11]*1.);
[ [ [ 1., 0., 0. ], [ 0., 1., 0. ], [ 0., 0., 0. ], [ 0.5, 0.5, 0. ], 
      [ 0.5, 0., 0. ], [ 0.833333, 0.166667, 0. ], [ 0.75, 0., 0. ], 
      [ 0.5, -0.5, 0. ], [ 1.5, 1.5, 0. ] ], 
  [ [ 1, 6, 7 ], [ 2, 3, 4 ], [ 3, 4, 5 ], [ 4, 5, 6 ], [ 4, 6, 9 ], 
      [ 5, 6, 7 ], [ 5, 7, 8 ] ] ]
gap> Retriangulation(s,list[12]*1.);
[ [ [ 2., 2., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 2., 1., 0. ], 
      [ 1., 1., 0. ], [ 1., 2., 0. ], [ 3., 1., 0. ], [ 1., 3., 0. ] ], 
  [ [ 1, 4, 6 ], [ 1, 4, 7 ], [ 1, 6, 8 ], [ 2, 4, 5 ], [ 3, 5, 6 ], 
      [ 4, 5, 6 ] ] ]
gap> Retriangulation(s,list[13]*1.);
[ [ [ -1., 0., 0. ], [ 1., 0., 0. ], [ 0., 1., 0. ], [ 0.5, 0., 0. ], 
      [ 0.666667, 0.333333, 0. ], [ 0., 0., 0. ], [ 0.5, 0.5, 0. ], 
      [ 0., -1., 0. ], [ 1., 1., 0. ], [ -1., -1., 0. ] ], 
  [ [ 1, 3, 6 ], [ 2, 4, 5 ], [ 3, 6, 7 ], [ 4, 5, 7 ], [ 4, 6, 7 ], 
      [ 4, 6, 8 ], [ 5, 7, 9 ], [ 6, 8, 10 ] ] ]
gap> Retriangulation(s,list[14]*1.);
[ [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 1., 2., 0. ], [ 1., 1., 0. ], 
      [ 1., 0., 0. ], [ 1.5, 1., 0. ], [ 1.66667, 0.666667, 0. ], 
      [ 2., 1., 0. ] ], 
  [ [ 1, 3, 4 ], [ 1, 4, 5 ], [ 2, 5, 7 ], [ 3, 4, 6 ], [ 4, 5, 7 ], 
      [ 4, 6, 7 ], [ 6, 7, 8 ] ] ]
gap> Retriangulation(s,list[15]*1.);
[ [ [ 0., 0., 0. ], [ 1., 1., 0. ], [ 0., 1., 0. ], [ 0.5, 0.5, 0. ], 
      [ 2., 0., 0. ] ], [ [ 1, 3, 4 ], [ 1, 4, 5 ], [ 2, 3, 4 ] ] ]
gap> Retriangulation(s,list[16]*1.);
[ [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 0.5, 0.5, 0. ], 
      [ 1.5, 0.5, 0. ], [ 0.5, 1.5, 0. ] ], 
  [ [ 1, 2, 4 ], [ 1, 3, 4 ], [ 2, 4, 5 ], [ 3, 4, 6 ], [ 4, 5, 6 ] ] ]
gap> Retriangulation(s,list[17]*1.);
[ [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 0.5, 0.5, 0. ], 
      [ 1., 0.5, 0. ], [ 0.5, 1., 0. ] ], 
  [ [ 1, 2, 5 ], [ 1, 3, 6 ], [ 1, 4, 5 ], [ 1, 4, 6 ], [ 2, 3, 5 ], 
      [ 3, 5, 6 ], [ 4, 5, 6 ] ] ]
gap> Retriangulation(s,list[18]*1.);
[ [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 1., 0.5, 0. ], 
      [ 0.5, 1., 0. ], [ 0.25, 0., 0. ], [ 0., 0.25, 0. ], [ -0.5, -0.5, 0. ] 
     ], [ [ 1, 6, 7 ], [ 1, 6, 8 ], [ 1, 7, 8 ], [ 2, 3, 4 ], [ 2, 4, 6 ], 
      [ 3, 4, 5 ], [ 3, 5, 7 ], [ 4, 5, 6 ], [ 5, 6, 7 ] ] ]
gap> Retriangulation(s,list[19]*1.);
[ [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 0., 2., 0. ], [ 1., 0.5, 0. ], 
      [ 0.5, 1., 0. ], [ 1.21429, 0.785714, 0. ], [ 0.785714, 1.21429, 0. ], 
      [ 2.5, 2.5, 0. ] ], 
  [ [ 1, 2, 4 ], [ 1, 3, 5 ], [ 1, 4, 5 ], [ 2, 4, 6 ], [ 3, 5, 7 ], 
      [ 4, 5, 6 ], [ 4, 5, 7 ], [ 4, 6, 7 ], [ 5, 6, 7 ], [ 6, 7, 8 ] ] ]
gap> Retriangulation(s,list[20]*1.);
[ [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 1., 2., 0. ], [ 0.857143, 1.71429, 0. ],
      [ 1., 1.5, 0. ], [ 1.14286, 1.71429, 0. ], [ 0., 3., 0. ], 
      [ 2., 3., 0. ] ], 
  [ [ 1, 2, 5 ], [ 1, 4, 5 ], [ 2, 5, 6 ], [ 3, 4, 6 ], [ 3, 4, 7 ], 
      [ 3, 6, 8 ], [ 3, 7, 8 ], [ 4, 5, 6 ] ] ]
gap> Retriangulation(s,list[21]*1.);
[ [ [ 0., 0., 0. ], [ 2., 0., 0. ], [ 1., 2., 0. ], [ 1., 0., 0. ], 
      [ -1., 0., 0. ], [ -2., -3., 0. ] ], 
  [ [ 1, 3, 4 ], [ 1, 4, 6 ], [ 1, 5, 6 ], [ 2, 3, 4 ] ] ]
gap> 