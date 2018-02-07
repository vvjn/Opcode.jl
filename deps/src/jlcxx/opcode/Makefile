SRCS = OPC_AABBCollider.cpp  OPC_Common.cpp         OPC_OBBCollider.cpp     OPC_SphereCollider.cpp  Opcode.cpp OPC_AABBTree.cpp      OPC_HybridModel.cpp    OPC_OptimizedTree.cpp   OPC_SweepAndPrune.cpp   StdAfx.cpp OPC_BaseModel.cpp     OPC_LSSCollider.cpp    OPC_Picking.cpp         OPC_TreeBuilders.cpp OPC_BoxPruning.cpp    OPC_MeshInterface.cpp  OPC_PlanesCollider.cpp  OPC_TreeCollider.cpp OPC_Collider.cpp      OPC_Model.cpp          OPC_RayCollider.cpp     OPC_VolumeCollider.cpp
OBJS = $(SRCS:.cpp=.o)
INCS = OPC_AABBCollider.h   OPC_IceHook.h         OPC_Picking.h            OPC_SphereAABBOverlap.h  OPC_TriTriOverlapGilvan.h OPC_AABBTree.h       OPC_LSSAABBOverlap.h  OPC_PlanesAABBOverlap.h  OPC_SphereCollider.h     OPC_VolumeCollider.h OPC_BaseModel.h      OPC_LSSCollider.h     OPC_PlanesCollider.h     OPC_SphereTriOverlap.h   Opcode.h OPC_BoxBoxOverlap.h  OPC_LSSTriOverlap.h   OPC_PlanesTriOverlap.h   OPC_SweepAndPrune.h      Stdafx.h OPC_BoxPruning.h     OPC_MeshInterface.h   OPC_RayAABBOverlap.h     OPC_TreeBuilders.h OPC_Collider.h       OPC_Model.h           OPC_RayCollider.h        OPC_TreeCollider.h OPC_Common.h         OPC_OBBCollider.h     OPC_RayTriOverlap.h      OPC_TriBoxOverlap.h OPC_HybridModel.h    OPC_OptimizedTree.h   OPC_Settings.h           OPC_TriTriOverlap.h

ICESRCS = Ice/IceAABB.cpp       Ice/IceIndexedTriangle.cpp  Ice/IceOBB.cpp    Ice/IceRandom.cpp          Ice/IceSegment.cpp Ice/IceContainer.cpp  Ice/IceMatrix3x3.cpp        Ice/IcePlane.cpp  Ice/IceRay.cpp             Ice/IceTriangle.cpp Ice/IceHPoint.cpp     Ice/IceMatrix4x4.cpp        Ice/IcePoint.cpp  Ice/IceRevisitedRadix.cpp  Ice/IceUtils.cpp
ICEOBJS = $(ICESRCS:.cpp=.o)
ICEINCS = Ice/IceAABB.h            Ice/IceFPU.h              Ice/IceMatrix3x3.h     Ice/IcePairs.h         Ice/IceRandom.h          Ice/IceTriList.h Ice/IceAxes.h            Ice/IceHPoint.h           Ice/IceMatrix4x4.h     Ice/IcePlane.h         Ice/IceRay.h             Ice/IceTriangle.h Ice/IceBoundingSphere.h  Ice/IceIndexedTriangle.h  Ice/IceMemoryMacros.h  Ice/IcePoint.h         Ice/IceRevisitedRadix.h  Ice/IceTypes.h Ice/IceContainer.h       Ice/IceLSS.h              Ice/IceOBB.h           Ice/IcePreprocessor.h  Ice/IceSegment.h         Ice/IceUtils.h


%.o: %.cpp %.h Stdafx.h
	g++ -I . -fPIC -c $< -o $@

all: $(OBJS) $(ICEOBJS)
	g++ -shared -fPIC $(OBJS) $(ICEOBJS) -o libOpcode.so


clean:
	rm -f *.o Ice/*.o libOpcode.so

