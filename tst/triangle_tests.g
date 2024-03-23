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

calculate_intersections(faces,list[1]*1.,false,Group(()));
calculate_intersections(faces,list[2]*1.,false,Group(()));
calculate_intersections(faces,list[3]*1.,false,Group(()));
calculate_intersections(faces,list[4]*1.,false,Group(()));
calculate_intersections(faces,list[5]*1.,false,Group(()));
calculate_intersections(faces,list[6]*1.,false,Group(()));
calculate_intersections(faces,list[7]*1.,false,Group(()));
calculate_intersections(faces,list[8]*1.,false,Group(()));
calculate_intersections(faces,list[9]*1.,false,Group(()));
calculate_intersections(faces,list[10]*1.,false,Group(()));
calculate_intersections(faces,list[11]*1.,false,Group(()));
calculate_intersections(faces,list[12]*1.,false,Group(()));
calculate_intersections(faces,list[13]*1.,false,Group(()));
calculate_intersections(faces,list[14]*1.,false,Group(()));
calculate_intersections(faces,list[15]*1.,false,Group(()));
calculate_intersections(faces,list[16]*1.,false,Group(()));
calculate_intersections(faces,list[17]*1.,false,Group(()));
calculate_intersections(faces,list[18]*1.,false,Group(()));
calculate_intersections(faces,list[19]*1.,false,Group(()));
calculate_intersections(faces,list[20]*1.,false,Group(()));
calculate_intersections(faces,list[21]*1.,false,Group(()));
