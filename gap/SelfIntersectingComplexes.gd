#############################################################################
##  SelfIntersectingComplexes.gd
#############################################################################
##
##  Declaration part of the package
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

## <#GAPDoc Label="SetNonManifoldVertexShift">
## <ManSection>
## <Func Name="SetNonManifoldVertexShift" Arg="shift"/>
## <Description>
##   Set the vertex shift parameter to <A>shift</A>.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("SetNonManifoldVertexShift");

## <#GAPDoc Label="SetNonManifoldEdgeShift">
## <ManSection>
## <Func Name="SetNonManifoldEdgeShift" Arg="shift"/>
## <Description>
##   Set the edge shift parameter to <A>shift</A>.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("SetNonManifoldEdgeShift");

## <#GAPDoc Label="SetPrecisionForIntersections">
## <ManSection>
## <Func Name="SetPrecisionForIntersections" Arg="eps"/>
## <Description>
##   Set the precision used for internal computations to <A>eps</A>.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("SetPrecisionForIntersections");

## <#GAPDoc Label="ExtractChamber">
## <ManSection>
## <Func Name="ExtractChamber" Arg="t,coordinates,f,n"/>
## <Description>
##   Extracts the chamber of a given triangular complex <A>t</A>, with embedding given by <A>coordinates</A>, containing face <A>f</A> with normal vector <A>n</n>.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("ExtractChamber");

## <#GAPDoc Label="OuterHull">
## <ManSection>
## <Func Name="OuterHull" Arg="t,coordinates,f,n"/>
## <Description>
##   Extracts the outer hull of a given triangular complex <A>t</A>, with embedding given by <A>coordinates</A>.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("OuterHull");

## <#GAPDoc Label="DrawSTLScratch">
## <ManSection>
## <Func Name="DrawSTLScratch" Arg="t,fileName,coordinates"/>
## <Description>
##   Generates STL-file with name given by the string <A>fileName</A> of the triangular complex <A>t</A> with embedding given by <A>coordinates</A>. As an STL-file requires normal vector for each face, random normal vectors will be used. If <A>t</A> is an oriented simplicial surface normal vectors are computed combinatorially.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("DrawSTLScratch");

## <#GAPDoc Label="DrawSTLwithNormals">
## <ManSection>
## <Func Name="DrawSTLScratch" Arg="t,fileName, coordinates,normals,visualize_normal_list"/>
## <Description>
##   Generates STL-file with name given by the string <A>fileName</A> of the triangular complex <A>t</A> with embedding given by <A>coordinates</A>. As an STL-file requires normal vector for each face, normal vectors are given in the list <A>normals</A> for each face. The list <A>visualize_normal_list</A> contains faces for which a normal vector is visualized using a triangle. Normally, this list will be kept empty and is primarly useful for debugging purposes.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("DrawSTLwithNormals");

## <#GAPDoc Label="DrawSurfaceToObj">
## <ManSection>
## <Func Name="DrawSurfaceToObj" Arg="t,fileName,coordinates"/>
## <Description>
##   Generates Obj-file with name given by the string <A>fileName</A> of the triangular complex <A>t</A> with embedding given by <A>coordinates</A>.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("DrawSurfaceToObj");

## <#GAPDoc Label="DrawSTLScratchVOF">
## <ManSection>
## <Func Name="DrawSTLScratchVOF" Arg="vof,fileName,coordinates"/>
## <Description>
##   Generates STL-file with name given by the string <A>fileName</A> of a triangular complex given by vertices of faces <A>vof</A> with embedding given by <A>coordinates</A>. As an STL-file requires normal vector for each face, random normal vectors will be used.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("DrawSTLScratchVOF");

## <#GAPDoc Label="DrawSTLwithNormalsVOF">
## <ManSection>
## <Func Name="DrawSTLwithNormalsVOF" Arg="vof,fileName, coordinates,normals,visualize_normal_list"/>
## <Description>
##   Generates STL-file with name given by the string <A>fileName</A> of a triangular complex given by its vertices of faces <A>vof</A> with embedding given by <A>coordinates</A>. As an STL-file requires normal vector for each face, normal vectors are given in the list <A>normals</A> for each face. The list <A>visualize_normal_list</A> contains faces for which a normal vector is visualized using a triangle. Normally, this list will be kept empty and is primarly useful for debugging purposes.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("DrawSTLwithNormalsVOF");

## <#GAPDoc Label="ComputeChambers">
## <ManSection>
## <Func Name="ComputeChambers" Arg="t,coordinates"/>
## <Description>
##   Computes all chambers of a given triangular complex <A>t</A> with embedding given by <A>coordinates</A>.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("ComputeChambers");

#DeclareGlobalFunction("AnimateChamberExplode");

## <#GAPDoc Label="RectifyDiscIntersections">
## <ManSection>
## <Func Name="RectifyDiscIntersections" Arg="l"/>
## <Description>
##   Given a list <A>l</A> that is given by vertices <A>l[1]</A> and edges connecting them <A>l[2]</A>. We compute all intersections of edges contained in <A>l[2]</A> and return a new list <C>l</D> which is geometrically equivalent but contains no self-intersecting edges.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("RectifyDiscIntersections");

## <#GAPDoc Label="DiscTriangulation">
## <ManSection>
## <Func Name="DiscTriangulation" Arg="l"/>
## <Description>
##   Given a list <A>l</A> that is given by vertices <A>l[1]</A> and edges connecting them <A>l[2]</A>. We assume that the first three vertices of <A>l[1]</A> span a triangle such that all other vertices are contained in this triangle. This function computes a triangulation and returns a disc <C>D</C> with vertices given by <A>l[1]</A> such that <C>D</C> contains all edges <A>l[2]</A>. This disc is obtained by iteratively adding the shortest possible edge without causing intersections.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("DiscTriangulation");

## <#GAPDoc Label="MergeOuterHull">
## <ManSection>
## <Func Name="MergeOuterHull" Arg="t,coordinates,name"/>
## <Description>
##   TODO
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("MergeOuterHull");

## <#GAPDoc Label="ComponentsOuterHull">
## <ManSection>
## <Func Name="ComponentsOuterHull" Arg="t,coordinates,name"/>
## <Description>
##   First compute a retriangulation and then fix degenerations of an embedded triangular complex <A>t</A> with coordinates <A>coordinates</A> and symmetry group acting on the faces together withs its orthogonoal embedding in <M>O(3)</M>. Moreover, an STL-file with name given by the string <A>name</A> is generated. 
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("ComponentsOuterHull");

## <#GAPDoc Label="JoinComponents">
## <ManSection>
## <Func Name="JoinComponents" Arg="t,coordinates,name"/>
## <Description>
##   TODO 
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("JoinComponents");

## <#GAPDoc Label="PrintableOuterHull">
## <ManSection>
## <Func Name="PrintableOuterHull" Arg="t,coordinates,name"/>
## <Description>
##   First compute a retriangulation and then fix degenerations of an embedded triangular complex <A>t</A> with coordinates <A>coordinates</A> and symmetry group acting on the faces together withs its orthogonoal embedding in <M>O(3)</M>. Moreover, an STL-file with name given by the string <A>name</A> is generated. 
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("PrintableSymmetricOuterHull");

## <#GAPDoc Label="PrintableOuterHull">
## <ManSection>
## <Func Name="PrintableOuterHull" Arg="t,coordinates,name"/>
## <Description>
##   First compute a retriangulation and then fix degenerations of an embedded triangular complex <A>t</A> with coordinates <A>coordinates</A>. Moreover, an STL-file with name given by the string <A>name</A> is generated. 
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("PrintableOuterHull");

## <#GAPDoc Label="ComputeSelfIntersections">
## <ManSection>
## <Func Name="ComputeSelfIntersections" Arg="t,coordinates"/>
## <Description>
##   Compute all self intersections of faces of an embedded triangular complex <A>t</A> with coordinates <A>coordinates</A>.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("ComputeSelfIntersections");

## <#GAPDoc Label="Retriangulation">
## <ManSection>
## <Func Name="Retriangulation" Arg="t,coordinates"/>
## <Description>
##   Compute retriangulation of an embedded triangular complex <A>t</A> with coordinates <A>coordinates</A>.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("Retriangulation");

## <#GAPDoc Label="ComponentsRetriangulation">
## <ManSection>
## <Func Name="ComponentsRetriangulation" Arg="t,coordinates"/>
## <Description>
##   Compute retriangulation of an embedded triangular complex <A>t</A> with coordinates <A>coordinates</A> for each component. Intersecting faces between components are retriangulated accordingly.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("ComponentsRetriangulation");

## <#GAPDoc Label="SymmetricRetriangulation">
## <ManSection>
## <Func Name="SymmetricRetriangulation" Arg="arg"/>
## <Description>
##   Compute retriangulation of an embedded triangular complex <A>t</A> with coordinates <A>coordinates</A>. Moreover, you can give either the symmetry group acting on the vertices or the symmetry group acting on the faces together withs its orthogonoal embedding in <M>O(3)</M>. This functions is faster than  if used correctly given correct input.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("SymmetricRetriangulation");

## <#GAPDoc Label="TwoTriangleIntersection">
## <ManSection>
## <Func Name="TwoTriangleIntersection" Arg="triangle1,triangle2"/>
## <Description>
##   Compute all intersections of two embedded triangles <A>triangle1,triangle2</A> and return it as a list of embedded edges and vertices describing the intersections.
## </Description>
## </ManSection>
## <#/GAPDoc>
DeclareGlobalFunction("TwoTriangleIntersection");



DeclareGlobalFunction("OrderNMEdges");

DeclareGlobalFunction("FixNMEdgePath");

DeclareGlobalFunction("FixNMVerts");

DeclareGlobalFunction("RemedyNonManifold");

# recursive function
DeclareGlobalFunction("_OrderPath");

# recursive function
DeclareGlobalFunction("_FixNMIntersection");

