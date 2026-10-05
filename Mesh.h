#ifndef _MESH_H_
#define _MESH_H_

#define _CRT_SECURE_NO_WARNINGS //i won't use fopen_s, sprintf_s, etc; so, don't warn me

//#define GRAPHICS_STUFF //define this to see/update the results on 3D viewer using Open Inventor library; undefine for a portable compilation

#define TINY						1e-6 //a small number bigger than DBL_MIN=2.2250738585072014e-308
#define INF							1e38 //largest float is 3.4x10^38; 1e+38 (or 1e38) also works (1e39 makes overflow; prints 1.#INF); DBL_MAX is 1.79769e+308 (any addition to DBL_MAX overflows so use INF instead which is still very big)
#define PI							3.141592653
#define EXP							2.718281828
#define USE_FIB_HEAP //much faster (70x faster, e.g., 55secs -> 0.8secs for 10 dijkstras on 52K-vertex mesh) than the array-based minimum search during dijkstraShortestPaths()

#include <time.h>
#include <iostream>
#include <vector>
#include <set>
#include <queue>
using namespace std;

//acceleration structures
#ifdef USE_FIB_HEAP
#include "FibHeap.h"
#endif

#ifdef GRAPHICS_STUFF
#define HAVE_INT8_T //similarly Eigen & Inventor/Inventor/system/inttypes.h has a redefinition problem: int8_t; so prevent them with this in files using include Inventor
#endif

struct Edge
{
	//edges[idx] is this edge
	int v1i, v2i;
	//distance between v1i & v2i
	double length, curvature;

	Edge(int i1, int i2, double l) : v1i(i1), v2i(i2), length(l) {};
};

struct Triangle
{
	//tris[idx] is this tri
	int v1i, v2i, v3i; //idx to verts forming this tri
	//normal and unit normal and unsigned area
	double* normal = NULL, * unormal, area;

	Triangle(int i1, int i2, int i3) : v1i(i1), v2i(i2), v3i(i3) {};

	bool setNormal(double* v1, double* v2, double* v3, bool print1)
	{
		if (! normal) //don't malloc if already malloced
		{
			normal = new double[3];
			unormal = new double[3];
		}

		//normal is the A, B, C triplet of Ax + By + Cz + D = 0; and unormal is the unitized version (direction)
		normal[0] = v1[1]*(v2[2]-v3[2]) + v2[1]*(v3[2]-v1[2]) + v3[1]*(v1[2]-v2[2]);
		normal[1] = v1[2]*(v2[0]-v3[0]) + v2[2]*(v3[0]-v1[0]) + v3[2]*(v1[0]-v2[0]);
		normal[2] = v1[0]*(v2[1]-v3[1]) + v2[0]*(v3[1]-v1[1]) + v3[0]*(v1[1]-v2[1]);
		//D = -v1[0]*(v2[1]*v3[2] - v3[1]*v2[2]) - v2[0]*(v3[1]*v1[2] - v1[1]*v3[2]) - v3[0]*(v1[1]*v2[2] - v2[1]*v1[2]);
		//D is needed only for distance from a point to tri computation; not needed in snake 3d

		//length of this vector
		double length = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

		//note that, contractToV2 or middle may produce temporary 0-area (due to collinearity) tris
		bool printed = false;
		if (length > 0)
		{
			unormal[0] = normal[0] / length;
			unormal[1] = normal[1] / length;
			unormal[2] = normal[2] / length; //frequent case
		}
		else
		{
			unormal[0] = unormal[1] = unormal[2] = 0.57735026918962576450914878050196; //rare case (this makes a unit normal but a constant one; better choice is to copy neighbor's unormal)
			if (print1)
			{
				cout << "WARNING: triangle w/ bad normal (there may be others; printing stopped)\n\n";
				printed = true;
			}
		}
		return printed;
	}

	void setArea(double* v1, double* v2, double* v3)
	{
		//if a triangle is defined by 3 points, say p, q and r, then its area is 0.5 * length of ((p - r) cross (q - r)) (see Real-Time Rendering book Appendix A)

		double* p_r = new double[3];
		for (int i = 0; i < 3; i++)
			p_r[i] = v2[i] - v1[i];

		double* q_r = new double[3];
		for (int i = 0; i < 3; i++)
			q_r[i] = v3[i] - v1[i];

		double* cross = new double[3];
		cross[0] = p_r[1]*q_r[2] - p_r[2]*q_r[1];
		cross[1] = p_r[2]*q_r[0] - p_r[0]*q_r[2];
		cross[2] = p_r[0]*q_r[1] - p_r[1]*q_r[0];

		//length of this vector
		double lengthSquared = cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2];
		if (lengthSquared >= 0.0) //seems redundant but observed -nan(ind) once in a deforming mesh's non-manifold triangle so prevent sqrt(-ve)
			area = 0.5 * sqrt(lengthSquared);
		else
			area = 0.0;

		//delete a, b, c; ignores b and c
		delete[] p_r;
		delete[] q_r;
		delete[] cross; //arrays allocated with new[] must be deallocated with delete[]
	}
};

struct Vertex
{
	int prev, binID, gtIdx, tmpIdx, matchIdx; //correspondence on the other mesh
	//spatial and embed, e.g., MDS, coordinates
	double* coords, * normal, * unormal, area, gc, patchGC, agd, tmp, dSpanned, distortion, distortionCpy, gtDistortion, geoToRef, geoToPath;
	vector< double > areas, areasRing, geoAll, geodesic; //double* geodesic does not release memory easily w/ classical delete[] geodesic (due to shallow copy etc) so use dynamic array instead; observed memory leaks w/ low-level double* geodesic approach so use vector instead
	float* color;
	bool sample, marked, pathVert, border, roi, roiBorder, coarseMapMatch, locMaxGC, locMaxGCrelaxed, geoNormalized; //geodesic from this vertex to all other are normalized (so don't normalize/divide again into tiny values during denseMapLargeN())
	//idx to vert neighbors of this vert (use set for efficient duplicate entry prevention); this is 1-ring neighborhood
	typedef set< int, std::less< int > > vertNeighborsIS; vertNeighborsIS vertNeighborsSet;
	vector< int > vertList1, triNeighbors, edgeList, sedgeList, patch, roiList; //1-ring neighbs, tri neighbs, edges/support edges incident to this vertex, etc
	vector< vector < int > > pathSamples, pathVerts, samplesToCoarse; //samples/verts on the path going from this sample to the jth sample are given by pathSamples/pathVerts[j]; similarly store samples from this sample to each coarseMap sample
#ifdef USE_FIB_HEAP
	FibHeapNode* heapNode; //corresponding heapNode for this vert to be used in DecreaseKey(heapNode, newKey)
#endif
	Vertex(double* c) : coords(c), color(NULL), gtIdx(-1), tmpIdx(-1), matchIdx(-1), coarseMapMatch(false), geoNormalized(false) {};

	bool addVertNeighbor(int vertInd)
	{
		//add vertInd as a vert neighb of this vert if it won't cause a duplication

		pair< vertNeighborsIS::const_iterator, bool > pa;
		//O(log N) search for duplicates by the use of red-black tree (balanced binary search tree)
		pa = vertNeighborsSet.insert(vertInd);
		//pa.second is false if item already exists, and true if insertion succeeds
		return pa.second;
	}

	void setNormal(double* nnormal)
	{
		normal = new double[3];
		unormal = new double[3];
		for (int i = 0; i < 3; i++)
			normal[i] = nnormal[i];
		//length of this vector
		double length = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
		//if (length >= ROUNDING_ERR) //0.00001f thrshld is bad 'cos 0.0000854606 and, undesiredly, 0.0000854606 >= 0.00001f holds
		if (length > 0)
		{
			unormal[0] = normal[0] / length;
			unormal[1] = normal[1] / length;
			unormal[2] = normal[2] / length;
		}
		else
		{
			unormal[0] = unormal[1] = unormal[2] = 0.57735026918962576450914878050196; //rare case (this makes a unit normal but a constant one; better choice is to copy neighbor's unormal)
			cout << "WARNING: vertex has bad normal!\t" << length << endl;
		}
	}
};

class Mesh
{
public:
	vector< Triangle* > tris;
	vector< Vertex* > verts;
	vector< Edge* > edges, sedges;
	vector< int > samples, samplesC2F, sphereIdxs, coarseMap; //store coarseMap to avoid a lot of coarseMap parameter passing from Correspondence.denseMap()
	size_t id, m1, m2, m3, nROIverts, nRingAreas, borderTreat, nonuniformPathSamplingHeuristic;
	float* color;
	double minEdgeLen, maxEdgeLen, avgEdgeLen, edgeLenTotal, maxGeoDist, avgGC, stdDevGC, totalArea, avgPatchArea, avgRingArea, avgPatchGC, radius, ** Q2D, *** Q, ** Q2Ddense, ***Qdense;
	bool printHeavy, mesh2, fpsFromGT, weightByGC, eucPatching, extraSamples; //extra path samples from a sample towards coarseMap samples (for the dense matching only where we have a non-empty coarseMap)
	Mesh(size_t i, size_t border, bool m2, size_t nonunifPath) : id(i), borderTreat(border), fpsFromGT(false), mesh2(m2), avgRingArea(0.0), nonuniformPathSamplingHeuristic(nonunifPath) { color = new float[3], color[0] = color[1] = color[2] = 0.7f, extraSamples = true, weightByGC = false, printHeavy = false; };
	void loadOff(char* meshFile, bool eucPatch = false);
	void loadGT(char* meshFile, Mesh* mesh1);
	void additionalRotations(int axis, double degree);
	void vertColors(int patchColoring, bool roiReady = true);
	void FPS(size_t N, bool curvatureBased);
	void FPSroi(size_t N, bool curvatureBased);
	void FPSroi(size_t k, const vector< int >& prevSamples);
	int FPSpatch(size_t n, double r, int v, bool curvatureBased, bool untouchedsOnly = false);
	void FPSfile(size_t N, bool competitorResult = false);
	void unitGeoScaling(double scaleFactor, bool updateMaxGeo = true);
	void curvaturesAll();
	void avgGeoDists();
	void areas();
	void areasForDisplay(int i);
	void doubleBorderAreas();
	void avgAreas(bool print = true);
	void locMaxSampling(int N, bool denseSampling = true);
	void dijkstraShortestPaths(int sourceVert);	
	void pathSampling(size_t n);
	void pathSamplingCoarse(size_t n, bool print = false);
	double ROI(const vector< int >& coarseMapFull, double upperBound = INF);
	void fillQ(size_t pathSize, double areaWeight, bool descriptorAvailable);
	void fillQdense(size_t pathSize, double areaWeight);
	void fillQdenseC2F(size_t pathSize, double areaWeight);
	void resultFromFile(Mesh* mesh2);
	void resultFromFile2(Mesh* mesh2);
	void resultFromFile3(Mesh* mesh2, int dl);
	void geoAll(vector< int > anchors, int nCloseAnchors);
private:
	void addVertexND(double* c);
	void addTriangle(int v1i, int v2i, int v3i);
	void addEdge(int v1i, int v2i);
	void addSupportEdge(int v1i, int v2i);
	void supportEdges();
	int borderVerts();	
	void unitAreaScaling();	
	void scaleEdges(double scaleFactor);
	void normalsAndAreas();
	void trianglesSharedBy(Edge* bigEdge, int* shTris);
	int edgeSharedBy(int vi, int vj);		
	int extremeVert(bool curvatureBased);
	void gaussianCurvatures();
	void patchAndDescriptors(int sourceVert);	
	void ROIij(int vi, int vj);
	vector< int > regionGrowing(int sourceVert, double mainArea);
	void clearSamples(size_t offset);
	bool replaceSample(int s);
	int localMaxGC(const vector< double >& curvs);
	void pathsFromVertex(int sourceVert);
	bool duplicateAhead(size_t k, size_t i, size_t j, int candidate);
};

#endif
