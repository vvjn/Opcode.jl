#include <vector>
#include <iostream>
#include <cmath>
#include "Opcode.h"
using namespace Opcode;

#include "jlcxx/array.hpp"
#include "jlcxx/jlcxx.hpp"

IceMaths::Point* create_ice_points(jlcxx::ArrayRef<double,2> verts, int64_t nverts) {
  IceMaths::Point *iceverts = new IceMaths::Point [nverts];
  for (int64_t ki=0; ki < nverts; ki++) {
    iceverts[ki] = IceMaths::Point(verts[3*ki], verts[3*ki+1], verts[3*ki+2]);
  }
  return iceverts;
}
IceMaths::IndexedTriangle* create_ice_triangles(jlcxx::ArrayRef<int64_t,2> faces, int64_t nfaces) {
  IceMaths::IndexedTriangle* icefaces = new IceMaths::IndexedTriangle [nfaces];
  for (int64_t ki=0; ki < nfaces; ki++) {
    icefaces[ki] = IceMaths::IndexedTriangle(faces[3*ki]-1, faces[3*ki+1]-1, faces[3*ki+2]-1);
  }
  return icefaces;
}

Model* aabb_tree_create(jlcxx::ArrayRef<double,2> verts, int64_t nverts, 
			jlcxx::ArrayRef<int64_t,2> faces, int64_t nfaces) {

  IceMaths::Point *cVertices = create_ice_points(verts, nverts);
  IceMaths::IndexedTriangle *cIndices = create_ice_triangles(faces, nfaces);

  Model *opcmodel = new Model;
  MeshInterface *cMesh = new MeshInterface;
  cMesh->SetInterfaceType(MESH_TRIANGLE);
  cMesh->SetNbTriangles(nfaces);
  cMesh->SetNbVertices(nverts);
  cMesh->SetPointers(cIndices, cVertices);
  OPCODECREATE OPCC;
  OPCC.mIMesh = cMesh;
  BuildSettings cBS;
  cBS.mRules = SPLIT_SPLATTER_POINTS | SPLIT_GEOM_CENTER;
  OPCC.mSettings = cBS;
  OPCC.mNoLeaf = true;
  OPCC.mQuantized = false;
  OPCC.mKeepOriginal = false; // For debug
  bool status = opcmodel->Build(OPCC);
  if (!status) {
    jl_error("Failed to build AABB tree");
  }
  return(opcmodel);
}

void aabb_tree_intersect_mut(Model *opcmodel, jlcxx::ArrayRef<double,2> origins,  // origin point of ray
			     jlcxx::ArrayRef<double,2> dirs,  // direction of ray
			     jlcxx::ArrayRef<int64_t> hits, // whether there is an intersection
			     jlcxx::ArrayRef<double> dists, // distance to intersection
			     jlcxx::ArrayRef<int64_t> faceinds, // index of the face that intersected the ray
			     jlcxx::ArrayRef<double,2> barys, // barycentric coordinates of intersection
			     int64_t nrays) {
  RayCollider RC;
  RC.SetFirstContact(false);
  RC.SetClosestHit(true);
  RC.SetCulling(false);
  float inf = MAX_FLOAT; //1.0f/0.0f;
  RC.SetMaxDist(inf);
  CollisionFaces CF;
  RC.SetDestination(&CF);

  const MeshInterface *mIMesh = opcmodel->GetMeshInterface();

  VertexPointers VP;
  IceMaths::Point v1,v2,v3, cStart, cDir;
  IceMaths::Ray cRay;
  double b1,b2,b3;
  const CollisionFace* colFaces;
  bool hit,status;
  static udword Cache;

  for(int64_t ki=0; ki<nrays; ki++) {
    cStart = IceMaths::Point(origins[3*ki], origins[3*ki+1], origins[3*ki+2]);
    cDir = IceMaths::Point(dirs[3*ki], dirs[3*ki+1], dirs[3*ki+2]);
    cDir.Normalize();
    cRay = IceMaths::Ray(cStart,cDir);
    
    status = RC.Collide(cRay, *opcmodel, NULL, &Cache);
    if (!status) jl_error("Error when hitting.");
    hit = RC.GetContactStatus();
    colFaces = CF.GetFaces();
    
    if (hit) {
      hits[ki] = ki+1;
      dists[ki] = colFaces[0].mDistance;
      faceinds[ki] = colFaces[0].mFaceID + 1;
      barys[2*ki] = colFaces[0].mU;
      barys[2*ki+1] = colFaces[0].mV;
    }
    else {
      hits[ki] = 0;
    }
  }
}
 

void aabb_tree_intersect(Model *opcmodel, jlcxx::ArrayRef<double,2> origins,  // origin point of ray
			 jlcxx::ArrayRef<double,2> dirs,  // direction of ray
			 jlcxx::ArrayRef<int64_t> hits, // whether there is an intersection
			 jlcxx::ArrayRef<double> dists, // distance to intersection
			 jlcxx::ArrayRef<int64_t> faceinds, // index of the face that intersected the ray
			 jlcxx::ArrayRef<double> barys,  // barycentric coordinates of intersection
			 int64_t nrays) {

  RayCollider RC;
  RC.SetFirstContact(false);
  RC.SetClosestHit(true);
  RC.SetCulling(false);
  float inf = MAX_FLOAT; //1.0f/0.0f;
  RC.SetMaxDist(inf);
  CollisionFaces CF;
  RC.SetDestination(&CF);

  const MeshInterface *mIMesh = opcmodel->GetMeshInterface();

  VertexPointers VP;
  IceMaths::Point v1,v2,v3, cStart, cDir;
  IceMaths::Ray cRay;
  double b1,b2,b3;
  const CollisionFace* colFaces;
  bool hit,status;
  static udword Cache;

  for(int64_t ki=0; ki<nrays; ki++) {
    cStart = IceMaths::Point(origins[3*ki], origins[3*ki+1], origins[3*ki+2]);
    cDir = IceMaths::Point(dirs[3*ki], dirs[3*ki+1], dirs[3*ki+2]);
    cDir.Normalize();
    cRay = IceMaths::Ray(cStart,cDir);
    
    status = RC.Collide(cRay, *opcmodel, NULL, &Cache);
    if (!status) jl_error("Error when hitting.");
    hit = RC.GetContactStatus();
    colFaces = CF.GetFaces();
    
    if (hit) {
      hits.push_back(ki+1);
      dists.push_back(colFaces[0].mDistance);
      faceinds.push_back(colFaces[0].mFaceID + 1);
      barys.push_back(colFaces[0].mU);
      barys.push_back(colFaces[0].mV);
    }
  }
}

void aabb_tree_deform(Model *opcmodel, jlcxx::ArrayRef<double,2> verts, int64_t nverts) {
  const MeshInterface *mIMesh = opcmodel->GetMeshInterface();
  if (nverts != mIMesh->GetNbVertices()) {
    jl_error("Input vertices need to match current mesh.");
  }
  IceMaths::Point* iceverts = (IceMaths::Point*)mIMesh->GetVerts();
  for (int ki=0; ki<nverts; ki++) {
    iceverts[ki] = IceMaths::Point(verts[ki*3],verts[ki*3+1],verts[ki*3+2]);
  }
  opcmodel->Refit();
}

void aabb_tree_delete(Model *opcmodel) {
  free(opcmodel);
}

JULIA_CPP_MODULE_BEGIN(registry)
jlcxx::Module& opcode = registry.create_module("Opcode");

opcode.add_type<Model>("Model");

opcode.method("aabb_tree_create", &aabb_tree_create);
opcode.method("aabb_tree_intersect", &aabb_tree_intersect);
opcode.method("aabb_tree_intersect_mut", &aabb_tree_intersect_mut);
opcode.method("aabb_tree_deform", &aabb_tree_deform);
opcode.method("aabb_tree_delete", &aabb_tree_delete);

JULIA_CPP_MODULE_END

