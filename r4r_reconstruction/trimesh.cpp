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

#include <assert.h>
#include <limits>
#include <algorithm>

#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>

#include "trimesh.h"

using namespace std;
using namespace OpenMesh;

namespace R4R {

CTriangleMesh::CTriangleMesh():
        TriangleMesh() {

}

float CTriangleMesh::CotanOppositeAngle(HalfedgeHandle heh) {

    float sp, cotanAngle;

	cotanAngle = 0;

    VertexHandle vh1, vh2, vh3;

	vh1 = from_vertex_handle(heh);
	vh2 = to_vertex_handle(heh);
	vh3 = opposite_vh(heh);

    if (vh3!=TriangleMesh::InvalidVertexHandle) {

        TriangleMesh::Point v1, v2, v3;
		v1 = point(vh1);
		v2 = point(vh2);
		v3 = point(vh3);

        TriangleMesh::Normal v3v1, v3v2;
		v3v1 = v1 - v3;
		v3v2 = v2 - v3;

		v3v1 = v3v1.normalize();
		v3v2 = v3v2.normalize();
		sp = OpenMesh::dot(v3v1,v3v2);

		cotanAngle = sp / sqrt(1 - sp * sp);

        if(std::isnan(cotanAngle))
            cotanAngle = sp / sqrt(1 - sp * sp + COTAN_EPSILON);

	}

	return cotanAngle;

}

float CTriangleMesh::FaceArea(FaceHandle fh) {

    return calc_sector_area(fh_iter(fh));

}

float CTriangleMesh::VertexQuadratureWeight(VertexHandle vh) {

    VertexFaceIter vf_it;

    float w = 0;

    for (vf_it = this->vf_iter(vh); vf_it; ++vf_it)
        w += this->FaceArea(vf_it.handle());

    return w/3.0;

}

bool CTriangleMesh::IsTriangleObtuse(FaceHandle fh) {

    bool obtuse = false;

    FaceHalfedgeIter foh_it;

    for (foh_it=this->fh_iter(fh); foh_it; ++foh_it) {

        if(this->calc_sector_angle(foh_it)>0.5*M_PI) {

            obtuse = true;
            break;

        }

    }

    return obtuse;

}

float CTriangleMesh::VoronoiArea(VertexHandle vh) {

    VertexOHalfedgeIter voh_it;

    float voronoiArea = 0;
    float denom, angle, triArea, previousLength, currentLength, previousWeight, currentWeight;

    TriangleMesh::Normal previousEdge, currentEdge;

    for (voh_it = voh_iter(vh); voh_it; ++voh_it) {

        if(face_handle(voh_it.handle()) != TriangleMesh::InvalidFaceHandle) {

            calc_edge_vector(voh_it, currentEdge);
            calc_edge_vector(cw_rotated_halfedge_handle(voh_it), previousEdge);

            denom = currentEdge.norm()*previousEdge.norm();

             // degenerate triangle, no area gain
            if(denom > 0) {

                angle = acos(OpenMesh::dot(currentEdge, previousEdge)/denom);

                    // obtuse triangle: area
                    if(IsTriangleObtuse(this->face_handle(voh_it.handle()))) {
                    //if(false) {

                        triArea = this->FaceArea(this->face_handle(voh_it.handle()));

                        // obtuse angle at center vertex
                        if(angle>0.5*M_PI)
                            voronoiArea = voronoiArea + 0.5*triArea;
                        else
                            voronoiArea = voronoiArea + 0.25*triArea;

                    }
                    else {

                        previousLength = this->calc_edge_sqr_length(this->cw_rotated_halfedge_handle(voh_it));
                        currentLength = this->calc_edge_sqr_length(voh_it);
                        previousWeight = this->CotanOppositeAngle(voh_it.handle());
                        currentWeight = this->CotanOppositeAngle(this->opposite_halfedge_handle(voh_it));

                        triArea = (currentLength*currentWeight + previousLength*previousWeight)/8.0;
                        voronoiArea = voronoiArea + triArea;

                    }

            }

        }

    }

    return voronoiArea;
}

vec3f CTriangleMesh::Barycenter(FaceHandle fh) {

    vec3f cog;

    float valence = 0;

    TriangleMesh::FaceVertexIter fv_it;

	for (fv_it = fv_iter(fh); fv_it; ++fv_it) {

        cog += this->Point(fv_it);
    	++valence;

	}

	cog = cog/valence;

	return cog;

}

vec3f CTriangleMesh::Barycenter() {

    vec3f barycenter;
    size_t counter = 0;

    TriangleMesh::VertexIter v_it;

    for (v_it=vertices_begin(); v_it!=vertices_end(); ++v_it) {

        barycenter = barycenter + Point(v_it);
        counter++;

    }

    return (1.0/float(counter))*barycenter;

}

TriangleMesh::Normal CTriangleMesh::DualEdgeVector(HalfedgeHandle heh) {

    assert(face_handle(heh)!=TriangleMesh::InvalidFaceHandle);

    HalfedgeHandle a, b;
	a = prev_halfedge_handle(heh);
	b = prev_halfedge_handle(a);

    TriangleMesh::Normal aVec, bVec;
	calc_edge_vector(a, aVec);
	calc_edge_vector(b, bVec);
	bVec = -bVec;

    TriangleMesh::Normal result = -(aVec*CotanOppositeAngle(a) + bVec*CotanOppositeAngle(b));

	return result;

}

HalfedgeHandle CTriangleMesh::OppositeHalfedgeHandle(FaceHandle fh, VertexHandle vh) {

    HalfedgeHandle result = TriangleMesh::InvalidHalfedgeHandle;

    TriangleMesh::FaceHalfedgeIter fh_it;

	for (fh_it = fh_iter(fh); fh_it; ++fh_it) {

		if(opposite_vh(fh_it) == vh) {

			result = fh_it;
			break;

		}

	}

	return result;

}


HalfedgeHandle CTriangleMesh::PrevBoundaryHalfedge(VertexHandle vh) {

    HalfedgeHandle result;

	if(!is_boundary(vh))
        return TriangleMesh::InvalidHalfedgeHandle;

    TriangleMesh::VertexOHalfedgeIter voh_it;

	for (voh_it = voh_iter(vh); voh_it; ++voh_it) {

        if(face_handle(voh_it.handle())==TriangleMesh::InvalidFaceHandle) {

			result = voh_it.handle();
			break;

		}

	}

	return opposite_halfedge_handle(result);

}

HalfedgeHandle CTriangleMesh::NextBoundaryHalfedge(VertexHandle vh) {

    HalfedgeHandle result;

	if(!is_boundary(vh))
        return CTriangleMesh::InvalidHalfedgeHandle;

    TriangleMesh::VertexOHalfedgeIter voh_it;

	for (voh_it = voh_iter(vh); voh_it; ++voh_it) {

        if(face_handle(opposite_halfedge_handle(voh_it.handle()))==TriangleMesh::InvalidFaceHandle) {

			result = voh_it.handle();

			break;

		}

	}

	return result;

}

VertexHandle CTriangleMesh::PreviousBoundaryVertex(VertexHandle vh) {

    HalfedgeHandle res = PrevBoundaryHalfedge(vh);

    if(res != TriangleMesh::InvalidHalfedgeHandle)
		return from_vertex_handle(res);
	else
        return TriangleMesh::InvalidVertexHandle;

}

VertexHandle CTriangleMesh::NextBoundaryVertex(VertexHandle vh) {

    HalfedgeHandle res = NextBoundaryHalfedge(vh);

    if(res != TriangleMesh::InvalidHalfedgeHandle)
		return to_vertex_handle(res);
	else
        return TriangleMesh::InvalidVertexHandle;

}

TriangleMesh::Normal CTriangleMesh::ExteriorNormal(VertexHandle vh) {

    HalfedgeHandle hehNext = NextBoundaryHalfedge(vh);
    HalfedgeHandle hehPrev = PrevBoundaryHalfedge(vh);

	return -((DualEdgeVector(hehNext).normalize() + DualEdgeVector(hehPrev).normalize())*0.5).normalize();

}

void CTriangleMesh::UniformMeshRefinement(size_t n) {

    OpenMesh::Subdivider::Uniform::LoopT<TriangleMesh> sd;

    sd.attach(*this);
    sd(n);
    sd.detach();

    cout << "No of vertices: " << n_vertices() << endl;

}

void CTriangleMesh::EdgeCollapse(double threshold) {

	size_t before = n_vertices();

	if (!has_vertex_status())
		request_vertex_status ();

	if (!has_edge_status())
		request_edge_status ();

	if (!has_face_status())
		request_face_status ();

	if (!has_halfedge_status())
		request_halfedge_status ();

    TriangleMesh::EdgeIter e_it;
    vector<EdgeHandle> toDel;

	double weight;


	for (e_it = edges_begin(); e_it != edges_end(); ++e_it) {

		weight = CotanOppositeAngle(halfedge_handle(e_it,0)) + CotanOppositeAngle(halfedge_handle(e_it,1));

        if(weight>threshold || std::isnan(weight))
			toDel.push_back(e_it.handle());

	}

	for(size_t i=0; i<toDel.size(); i++) {

		if(is_collapse_ok(halfedge_handle(toDel[i],0)))
			collapse(halfedge_handle(toDel[i],0));

		if(is_collapse_ok(halfedge_handle(toDel[i],1)))
			collapse(halfedge_handle(toDel[i],1));

	}


	garbage_collection();

	release_vertex_status();
	release_edge_status();
	release_halfedge_status();
	release_face_status();

	size_t after = n_vertices();

	cout << "No of deleted vertices: " << before-after << std::endl;

}

bool CTriangleMesh::SmoothBoundary(size_t n) {

    VertexHandle init, actual, prev, next;
    TriangleMesh::Point newPt;

	init = FindBoundaryVertex();

	for(size_t k=0; k<n; k++) {

		cout << "Smoothing step k=" << k << endl;

		actual = init;

		do {

            if (actual == TriangleMesh::InvalidVertexHandle) {

				cout << "WARNING: Mesh has no boundary..." << endl;
				return 1;

			}

			prev = PreviousBoundaryVertex(actual);
			next = NextBoundaryVertex(actual);

			newPt = (point(prev) + point(next))*0.5;

			point(actual) = newPt;

			actual = next; //NextBoundaryVertex(actual);

		}
		while (actual!=init);

	}

	return 0;

}

VertexHandle CTriangleMesh::FindBoundaryVertex() {

    VertexHandle result = CTriangleMesh::InvalidVertexHandle;

    TriangleMesh::VertexIter v_it, v_end(vertices_end());

	for (v_it = vertices_begin(); v_it != v_end; ++v_it) {

		if(is_boundary(v_it)) {

			result = v_it.handle();
			break;

		}

	}

	return result;

}

void CTriangleMesh::DeleteVertices(set<VertexHandle> vertices) {

	if (!has_vertex_status())
		request_vertex_status ();

	if (!has_edge_status())
		request_edge_status ();

	if (!has_face_status())
		request_face_status ();

	if (!has_halfedge_status())
		request_halfedge_status ();

    TriangleMesh::VertexIter v_it, v_end(vertices_end());

	for (v_it = vertices_begin(); v_it != v_end; ++v_it) {

		//if(property(toDelete, v_it))
		if(vertices.find(v_it) != vertices.end())
			delete_vertex(v_it, true);

	}

	garbage_collection();

	release_vertex_status();
	release_edge_status();
	release_halfedge_status();
	release_face_status();
	//remove_property(toDelete);

}

void CTriangleMesh::DeleteFaces(set<FaceHandle> faces) {

	if (!has_vertex_status())
		request_vertex_status ();

	if (!has_edge_status())
		request_edge_status ();

	if (!has_face_status())
		request_face_status ();

	if (!has_halfedge_status())
		request_halfedge_status ();

    TriangleMesh::FaceIter f_it, f_end(faces_end());

	for (f_it = faces_begin(); f_it != f_end; ++f_it) {

		if(faces.find(f_it) != faces.end())
			delete_face(f_it, true);

	}

	garbage_collection();

	release_vertex_status();
	release_edge_status();
	release_halfedge_status();
	release_face_status();

}

void CTriangleMesh::Disturb(float mu, float sigma) {

    // generate random normal displacement
    vecf noise(this->n_vertices());
    noise.RandN(mu,sigma);

    TriangleMesh::VertexIter v_it;

    this->update_normals();

    for (v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
        this->point(v_it) = this->point(v_it) + this->normal(v_it)*noise.Get(v_it.handle().idx());

    this->update_normals();

}

void CTriangleMesh::SimpleSmooth(u_int n, bool boundary) {

    OpenMesh::VPropHandleT<TriangleMesh::Point> newPoints;
    this->add_property(newPoints);

    TriangleMesh::VertexIter v_it, v_end(this->vertices_end());
    TriangleMesh::VertexVertexIter vv_it;

    TriangleMesh::Point sumPoints;

    for (int i=0; i<n; ++i) {

        for (v_it = this->vertices_begin(); v_it != v_end; ++v_it) {

            sumPoints = TriangleMesh::Point(0.0,0.0,0.0);

            for (vv_it = this->vv_iter(v_it); vv_it; ++vv_it)
                sumPoints = sumPoints + this->point(vv_it);

            this->property(newPoints, v_it) = sumPoints/float(this->valence(v_it));

        }

        for (v_it = this->vertices_begin(); v_it != v_end; ++v_it) {

            if(!boundary) {

                if(!this->is_boundary(v_it))  // shrinking?
                    this->point(v_it) = this->property(newPoints, v_it);

            }
            else
                this->point(v_it) = this->property(newPoints, v_it);

        }

    }

    this->remove_property(newPoints);

}

void CTriangleMesh::ComputeGradientOperator(smatf& nabla) {

    nabla = smatf(3*this->n_faces(),this->n_vertices());

    // row pointer
    size_t row = 0;

    TriangleMesh::FaceIter f_it;
    TriangleMesh::FaceVertexIter fv_it;

    for(f_it=this->faces_begin(); f_it!=this->faces_end(); ++f_it) {

        // face area
        float A =  this->calc_sector_area(this->fh_iter(f_it));
        A = sqrt(A);

        if(A>0) {

            for (fv_it = this->fv_iter(f_it); fv_it; ++fv_it) {

                TriangleMesh::Normal grad = this->DualEdgeVector(this->OppositeHalfedgeHandle(f_it,fv_it));

                // gradient is (0.5/A) times sqrt(A) for integration = (0.5/A)*sqrt(A) = 1/2*sqrt(A)
                grad *= (0.5/A);

                // find col index
                size_t col = fv_it.handle().idx();

                for(uint i=0; i<3; i++)
                    nabla.Set(row+i,col,grad[i]);

            }

        }

        row += 3;

    }

}



vec3f CTriangleMesh::Point(VertexHandle vh) {

    vec3f point = { this->point(vh)[0], this->point(vh)[1], this->point(vh)[2] };

    return point;

}

vec3f CTriangleMesh::Normal(VertexHandle vh) {

    vec3f normal = { this->normal(vh)[0], this->normal(vh)[1], this->normal(vh)[2] };

    return normal;

}

vec3f CTriangleMesh::Normal(FaceHandle fh) {

    vec3f normal = { this->normal(fh)[0], this->normal(fh)[1], this->normal(fh)[2] };

    return normal;

}

void CTriangleMesh::BoundingBox(vec3f& lower, vec3f& upper) {

    for(u_int i=0; i<3; i++) {

        lower(i) = std::numeric_limits<float>::max();
        upper(i) = -std::numeric_limits<float>::max();

    }

    TriangleMesh::VertexIter v_it;

    for (v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it) {

        vec3f point = this->Point(v_it);

        for(u_int i=0; i<3; i++) {

            if(point.Get(i)>=upper.Get(i))
                upper(i) = point.Get(i);

            if(point.Get(i)<=lower.Get(i))
                lower(i) = point.Get(i);

        }

    }

}


CBoundingBox<float> CTriangleMesh::BoundingBox() {

    vec3f lower, upper;

    this->BoundingBox(lower,upper);

    return CBoundingBox<float>(lower,upper);

}

float CTriangleMesh::MeanEdgeLength() {

    TriangleMesh::EdgeIter e_it;

    float length = 0;

    for (e_it = this->edges_begin(); e_it != this->edges_end(); ++e_it)
        length += this->calc_edge_sqr_length(e_it);

    return length/float(this->n_edges());

}


void CTriangleMesh::Deform(VertexHandle vh, size_t n, float h) {

    // first get neighborhood
    // then distances in neighborhood
    // then exponential function of -d (sigma \propto n)

}

void CTriangleMesh::GetNRingNeighborhood(VertexHandle vh, size_t h, set<VertexHandle>& neighbors) {

    neighbors.insert(vh);

    // base case of recursion
    if(h==0)
        return;

    TriangleMesh::VertexVertexIter vv_it;
    for (vv_it = this->vv_iter(vh); vv_it; ++vv_it) {

        if(neighbors.find(vv_it.handle()) == neighbors.end())
            GetNRingNeighborhood(vv_it.handle(),h-1,neighbors);

    }

}

map<VertexHandle,float> CTriangleMesh::Dijkstra(VertexHandle vh, const set<VertexHandle>& vertices) {

    map<VertexHandle,float> result;     // distance container
    map<VertexHandle,bool> visited;     // visited vertices
    vector<vd> vds;                     // heap storage

    // if the seed is not even in the subset of vertices, do nothing
    if(vertices.find(vh)==vertices.end())
        return result;

    set<VertexHandle>::const_iterator it;
    for(it=vertices.begin(); it!=vertices.end(); ++it) {

        if((*it)==vh) {

            result.insert(pair<VertexHandle,float>(vh,0));
            visited.insert(pair<VertexHandle,bool>(vh,true));
            vds.push_back(pair<VertexHandle,float>(vh,0));

        }
        else {

            result.insert(pair<VertexHandle,float>(*it,std::numeric_limits<float>::max()));
            visited.insert(pair<VertexHandle,bool>(*it,false));

        }

    }

    // move the initial value to the top of the heap
    std::make_heap(vds.begin(),vds.end(),CVertexDistanceComparator());

    while(!vds.empty()) {

        // get the current smallest element
        vd current = vds.front();

        // remove it from the priority queue
        std::pop_heap(vds.begin(),vds.end(),CVertexDistanceComparator());
        vds.pop_back();

        TriangleMesh::VertexOHalfedgeIter voh_it;
        for (voh_it=this->voh_iter(current.first); voh_it; ++voh_it) {

            VertexHandle tv = this->to_vertex_handle(voh_it.handle());

            if(vertices.find(tv)!=vertices.end() && !visited[tv]) {

                float tentdist = this->calc_edge_length(voh_it) + result[current.first];

                if(tentdist<result[tv]) {

                    // set trsult
                    result[tv]  = tentdist;

                    // update the heap
                    vds.push_back(pair<VertexHandle,float>(tv,tentdist));
                    std::make_heap(vds.begin(),vds.end(),CVertexDistanceComparator());

                }

            }

        }

    }

    return result;

}




/*void CTriangleMesh::Diffuse(VertexHandle vh, size_t h, std::map<int,float>& funcin, std::map<int,float>& funcout) {

    TriangleMesh::VertexVertexIter vv_it;
    float sum = 0;
    for (vv_it = this->vv_iter(vh); vv_it; ++vv_it) {
        sum += funcin[vv_it.handle().idx()];
    }
    funcout[vh.idx()] = funcin[vh.idx()] + sum;
    funcout[vh.idx()] /= float(this->valence(vh));

    // base case of recursion
    if(h==0)
        return;

    for (vv_it = this->vv_iter(vh); vv_it; ++vv_it) {

        if(funcout.find(vv_it.handle().idx())==funcout.end())
            Diffuse(vv_it.handle(),h-1,funcin,funcout);

    }

}*/



}
