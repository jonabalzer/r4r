/*////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013, Jonathan Balzer
//
// All rights reserved.
//
// This file is part of the R4R library.
//
// The R4R library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The R4R library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the R4R library. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////*/

#ifndef TRIMESH_H_
#define TRIMESH_H_

#define COTAN_EPSILON 1e-10

#include <set>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Subdivider/Adaptive/Composite/CompositeT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::Subdivider::Adaptive::CompositeTraits> TriangleMesh;

#include "bbox.h"
#include "types.h"

namespace R4R {

/*! \brief triangle mesh
 *
 *
 *
 */
class CTriangleMesh:public TriangleMesh {

public:

	//! Standard constructor.
	CTriangleMesh();

    //! Area.
    float FaceArea(FaceHandle fh);

     //! Computes the integration weight for a vertex.
    float VertexQuadratureWeight(VertexHandle vh);

    //! Computes the Voronoi area around a vertex.
    float VoronoiArea(VertexHandle vh);

	//! Returns barycenter of a mesh face.
    vec3f Barycenter(FaceHandle fh);

	//! Exterior normal at a boundary vertex.
    CTriangleMesh::Normal ExteriorNormal(VertexHandle vh);

	//! Refines the mesh uniformly with Loop's scheme.
    void UniformMeshRefinement(size_t n);

    //! Collapses edges with aspect ratio to high to do FE-analysis.
    void EdgeCollapse(double threshold);

    /*! Smoothes the boundary of mesh.
     *
     * \details Only works for surfaces with topology of a disk.
     *
     */
    bool SmoothBoundary(size_t n);

    /*! Finds a vertex handle that is on the boundary of the mesh.
     *
     * \details If the boundary is multiply-connected, the first hit will be returned.
     */
    OpenMesh::VertexHandle FindBoundaryVertex();

    //! Deletes a given set of vertices.
    void DeleteVertices(std::set<VertexHandle> vertices);

    //! Deletes a given set of faces.
    void DeleteFaces(std::set<FaceHandle> faces);

    //! Random normal displacement.
    void Disturb(float mu, float sigma);

    //! Simple mesh smoothing by vertex averaging.
    void SimpleSmooth(u_int n, bool boundary);

    //! Casts point into R4R data structure.
    vec3f Point(VertexHandle vh);

    //! Casts normal into R4R data structure.
    vec3f Normal(VertexHandle vh);

    //! Casts normal into R4R data structure.
    vec3f Normal(FaceHandle vh);

    //! Computes the bounding box of the mesh.
    void BoundingBox(vec3f& lower, vec3f& upper);

    //! Computes bounding box of the mesh.
    CBoundingBox<float> BoundingBox();

    //! Average edge length.
    float MeanEdgeLength();

    //! Compute dual half edge.
    TriangleMesh::Normal DualEdgeVector(HalfedgeHandle heh);

    //! Gets handle of the half edge opposite to a vertex on a face.
    OpenMesh::HalfedgeHandle OppositeHalfedgeHandle(FaceHandle fh, VertexHandle vh);

    //! Computes the discrete gradient operator for vertex-based functions.
    void ComputeGradientOperator(R4R::smatf& nabla);

    //! Cotangent of the angle opposite to a half edge.
    float CotanOppositeAngle(HalfedgeHandle heh);

private:

    //! Previous half edge on boundary from vertex.
    OpenMesh::HalfedgeHandle PrevBoundaryHalfedge(VertexHandle vh);

    //! Next half edge on boundary from vertex.
    OpenMesh::HalfedgeHandle NextBoundaryHalfedge(VertexHandle vh);

    //! Next vertex on boundary from vertex.
    OpenMesh::VertexHandle PreviousBoundaryVertex(VertexHandle vh);

    //! Next vertex on boundary from vertex.
    OpenMesh::VertexHandle NextBoundaryVertex(VertexHandle vh);

    //! Checks for an obtuse triangle.
    bool IsTriangleObtuse(FaceHandle fh);

};

}

#endif /* TRIMESH_H_ */
