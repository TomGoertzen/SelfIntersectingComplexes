#############################################################################
##  testall.g
#############################################################################
##
##  This file runs package tests.
##
#############################################################################
##
##  This file is part of the SelfIntersectingComplexes package.
##
##  This file's authors include Christian Amend and Tom Goertzen.
##
##  Please refer to the COPYRIGHT file for details.
##
##  SPDX-License-Identifier: GPL-2.0-or-later
##
#############################################################################

LoadPackage("SelfIntersectingComplexes");

TestDirectory(DirectoriesPackageLibrary( "SelfIntersectingComplexes", "tst/files" ),
  rec(exitGAP := true));


##############################################################################################
##########
########
######## test for rectification of non-manifold parts

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
#### can test other functionalities here as well
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
f := Remedy_NonManifold(datas,points_fix,shift_param);
DrawSTLwithNormals(f[1],"ico_3_1_path_repaired",f[2],datas[4],[]);

######## 
########
##########
##############################################################################################

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
