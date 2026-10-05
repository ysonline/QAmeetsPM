#include "Mesh.h"

void Mesh::loadOff(char* meshFile, bool eucPatch)
{
	cout << "Mesh initializing (file: " << meshFile << ")...\n";

	FILE* fPtr;
	if (! (fPtr = fopen(meshFile, "r")))
		cout << "cannot read " << meshFile << endl, exit(0);
	char off[25];
	fscanf(fPtr, "%s\n", &off); //cout << off << " type file\n";
	float a, b, c, d;	//for face lines and the 2nd line (that gives # of verts, faces, and edges) casting these to int will suffice
	fscanf(fPtr, "%f %f %f\n", &a, &b, &c);
	size_t nVerts = (size_t) a, v = 0;

	minEdgeLen = INF;
	maxEdgeLen = -INF;
	edgeLenTotal = 0.0;
	double upscale = 1;//(rescale ? 1.2 : 1); //rescale pfaust models
	while (v++ < nVerts) //go until the end of coords section
	{
		fscanf(fPtr, "%f %f %f\n", &a, &b, &c);
		double* coords = new double[3];
		coords[0] = a * upscale;
		coords[1] = b * upscale;
		coords[2] = c * upscale;
		addVertexND(coords); //ND: no duplicate check
	}
	//verts ready, time to fill triangles
	while (fscanf(fPtr, "%f %f %f %f\n", &d, &a, &b, &c) != EOF) //go until the end of file
		addTriangle((int) a, (int) b, (int) c); //no -1 'cos idxs start from 0 for off files
	avgEdgeLen = edgeLenTotal / edges.size();
	fclose(fPtr);
//	if (! pfaust) //shrec16 mode requires x-axis 90 degree rotation for visualization only
		additionalRotations(0, -90);

	//vertNeighborsSet used for efficient duplicate handling in addVertNeighbor(); now transfer this info to vector for simpler/shorter access via v.vertList1
	for (v = 0; v < verts.size(); v++)
	{
		set< int, std::less< int > >::const_iterator neighbVertIdx;
		for (neighbVertIdx = verts[v]->vertNeighborsSet.begin(); neighbVertIdx != verts[v]->vertNeighborsSet.end(); neighbVertIdx++)
			verts[v]->vertList1.push_back(*neighbVertIdx); //1-ring neighbors in v.vertList1
	}
	//normal and area info for vertices and triangles
	normalsAndAreas();
//	unitAreaScaling(); //may be good for HKS but bad for geodesics on partial matching (no problem for complete/nonpartial cases)
	int nb = borderVerts(); //areasRing[] descriptor should be doubled for border vertices as the surface around them are missing/partial
	cout << "Mesh " << id << " has " << tris.size() << " tris, " << verts.size() << " verts, " << edges.size() << " edges";
	if (printHeavy)
		cout << "\nmin <= avg <= maxEdgeLen: " << minEdgeLen << " <= " << avgEdgeLen << " <= " << maxEdgeLen << ", totalArea: " << totalArea << ", nBorderVerts: " << nb << "\n";
	else
		cout << "; avgEdgeLen: " << avgEdgeLen << ", totalArea: " << totalArea << endl;
	
//	if (id == 10)
//		eucPatch = true;
//	eucPatching = (verts.size() > 20000);//19000); //geodesic-based patching, despite the early break in patchAndDescriptors(), is significantly slow (13secs on 52K-vertex mesh w/ nSamples=nPathSamples=10 (eucPatching takes just 0.3secs), 20secs nSamples=30,nPathSamples=3) so use euclidean distance during patching for a 40x speedup
//	eucPatching = false;
	eucPatching = eucPatch;
	if (eucPatching)
		cout << "WARNING: path patches by euclidean distance for efficiency; geodesic is more accurate but slow\n\n";

	//additional support edges for more accurate geodesics
	supportEdges();
	if (printHeavy)
		cout << endl;
}

void Mesh::loadGT(char* gtFile, Mesh* mesh1)
{
	//loads ground-truth matches available for TOSCA-partial for evaluation only

	cout << "ground-truths for eval.. ";
	int nGTs = 0;
	FILE* fPtr;
	if (! (fPtr = fopen(gtFile, "r")))
	{
		if (mesh1->verts.size() == verts.size()) //assume i-i is the ground-truth for all i
			for (size_t v = 0; v < verts.size(); v++)
				mesh1->verts[v]->gtIdx = v, verts[v]->gtIdx = v, nGTs++;
		else
			cout << "WARNING: no gt set for evaluation\n\n";
		cout << nGTs << " gt matches available\n";
		//return;//exit(0); don't exit 'cos gt evaluation is optional
	}
	else
	{
		int v2 = 0, gtMatchIn1; //mesh1 is the null pose (w/ id=0) so gtMatchIn1-v2 means mesh2.v2 has the ground-truth match of gtMatchIn1
		while (fscanf(fPtr, "%d\n", &gtMatchIn1) != EOF) //go until the end of file
			gtMatchIn1--, mesh1->verts[gtMatchIn1]->gtIdx = v2, verts[v2++]->gtIdx = gtMatchIn1, nGTs++; //indices in file start from 1; gtIdx=-1 defaulted in constructor so some mesh1 verts will remain as gtIdx=-1 (no gt match for them)
		cout << nGTs << " gt matches available through " << gtFile << endl;
		fclose(fPtr);
	}	
}

inline double distanceBetween(double* v1, double* v2)
{
	return sqrt( (v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]));
}
inline double distanceBetween2(double* v1, double* v2)
{
	//returns squared distance which skips the sqrt part of distanceBetween() for efficiency
	return (v2[0]-v1[0])*(v2[0]-v1[0]) +  (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]);
}
inline double dotProduct(double* v1, double* v2)
{
	//returns the dot product b/w v1 & v2
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}
inline void normalize(double* v1)
{
	double len = sqrt( dotProduct(v1, v1) );
	if (len > TINY)
		v1[0] /= len, v1[1] /= len, v1[2] /= len;
}
/*inline double angleBetween(double* v1, double* v2)
{
	//returns the smallest radian angle b/w normalized/unit vectors v1 & v2 via v1 DOT v2 = v1Len * v2Len * cosTheta (len = 1 for normalized, hence acos(v1 DOT v2) is the desired theta)

	double cosTheta = dotProduct(v1, v2);
	if (cosTheta < -1.0)	cosTheta = -1.0;
	if (cosTheta > 1.0)		cosTheta = 1.0;
	return acos(cosTheta); //* (180/PI); no radian to float to save time (besides gaussianCurvature() expects radian output)
}angleFromCotan more robust: Discrete Differential Ops for Triangulated 2-Manifolds paper Table 1 caption*/
inline double angleFromCotan(double* v1, double* v2)
{
	//returns the smallest radian angle b/w unnormalized v1 & v2 via safer atan2 (angle from cotan tactic here is supposed to be more robust according to Discrete Differential Ops for Triangulated 2-Manifolds paper Table 1 caption that prefers atan2 over acos)
	//also v1 & v2 do not have to be normalized here (sqrt(-ve) may arise due to precision when normalized v1/2 used 'cos x in sqrt(x-y) below is always 1 (small) when normalized v1/2 in use)

	double v1DotV2 = dotProduct(v1, v2), len1sqLen2sq = dotProduct(v1, v1) * dotProduct(v2, v2), v1DotV2sq = v1DotV2 * v1DotV2;
	if (len1sqLen2sq <= v1DotV2sq)
		//cout << "WARNING: " << len1sqLen2sq - v1DotV2sq << " inside sqrt\n\n",
		len1sqLen2sq = v1DotV2sq + TINY; //prevent negative param to sqrt() and also that sqrt() term goes to denominator in atan2 so prevent it from being 0 too
	//double denom = sqrt(len1sqLen2sq - v1DotV2sq); //cotan b/w v1 & v2 copied from cotan()
	return abs(atan2(sqrt(len1sqLen2sq - v1DotV2sq), v1DotV2)); //tan = denom / v1DotV2 = y/x*/
}

void Mesh::resultFromFile(Mesh* mesh2)
{
	//loads matches of the Partial 3-D Correspondence from Shape Extremities paper

	char fName[250];	
//	id = 10, mesh2->id = 19;
	sprintf(fName, "p3d-results\\%d-%dp1.txt", id, mesh2->id); //pfaust partialID=1
//	id = 30, mesh2->id = 11;
	sprintf(fName, "p3d-results\\%d-%dp1set3.txt", id, mesh2->id); //pfaust partialID=1
	sprintf(fName, "p3d-results\\%d-%dp1set2.txt", id, mesh2->id); //pfaust partialID=1
	sprintf(fName, "p3d-results\\%d-%dp1.txt", id, mesh2->id); //pfaust partialID=1
	sprintf(fName, "p3d-results\\%d-%dp3.txt", id, mesh2->id); //partialID=3
	sprintf(fName, "p3d-results\\%d-%dp3set3.txt", id, mesh2->id); //partialID=3
	sprintf(fName, "p3d-results\\%d-%dp1.txt", id, mesh2->id); //partialID=1
	sprintf(fName, "p3d-results\\%d-%dp1set2.txt", id, mesh2->id); //partialID=1
	sprintf(fName, "p3d-results\\%d-%dp4set2.txt", id, mesh2->id); //partialID=4
	sprintf(fName, "p3d-results\\%d-%dp4.txt", id, mesh2->id); //partialID=4
	sprintf(fName, "p3d-results\\%d-%dp2.txt", id, mesh2->id); //partialID=2
	sprintf(fName, "p3d-results\\%d-%dp2set3.txt", id, mesh2->id); //partialID=2
	sprintf(fName, "p3d-results\\%d-%dp2.txt", id, mesh2->id); //partialID=2
	sprintf(fName, "p3d-results\\%d-%dp2set2.txt", id, mesh2->id); //partialID=2	
	sprintf(fName, "p3d-results\\%d-%dp2set3.txt", id, mesh2->id); //partialID=2
	sprintf(fName, "p3d-results\\%d-%dp3.txt", id, mesh2->id); //partialID=3
	sprintf(fName, "p3d-results\\%d-%dp3set3.txt", id, mesh2->id); //partialID=3
	sprintf(fName, "p3d-results\\%d-%dp3set2.txt", id, mesh2->id); //partialID=3
	FILE* fPtr;
	if (! (fPtr = fopen(fName, "r")))
	{
		cout << "WARNING: not found: " << fName << "\n\n";
		return; //exit(0); don't exit 'cos gt evaluation is optional
	}
	int v1, v2;
	while (fscanf(fPtr, "%d\t%d\n", &v1, &v2) != EOF) //go until the end of file
		samples.push_back(v1), verts[v1]->sample = true, verts[v1]->matchIdx = v2, mesh2->samples.push_back(v2), mesh2->verts[v2]->sample = true, mesh2->verts[v2]->matchIdx = v1;
	fclose(fPtr);
}

void Mesh::resultFromFile2(Mesh* mesh2)
{
	//loads matches of the qmatch-style matcher

	char fName[250];	
	sprintf(fName, "qmatchStyle-results\\%d-%d n50 k5 solqmatch.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n50 k5 solqucoop.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n50 k5 solqucm.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n50 k20 solqmatch.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n50 k20 solqucoop.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n50 k20 solqucm.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n50 k40 solqmatch.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n250 k40 solqmatch.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n250 k40 solqucoop.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n250 k40 solqucm.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n250 k20 solqmatch.txt", id, mesh2->id);
	sprintf(fName, "qmatchStyle-results\\%d-%d n250 k40 solqucm2.txt", id, mesh2->id);
	FILE* fPtr;
	if (! (fPtr = fopen(fName, "r")))
	{
		cout << "WARNING: not found: " << fName << "\n\n";
		return; //exit(0); don't exit 'cos gt evaluation is optional
	}
	int v1 = 0, v2;
	while (fscanf(fPtr, "%d\n", &v2) != EOF) //go until the end of file
		verts[ samples[v1] ]->matchIdx = mesh2->samples[v2], mesh2->verts[ mesh2->samples[v2++] ]->matchIdx = samples[v1++];
	fclose(fPtr);
	//adapted from mapQualityGT()
	double gtDistortion = 0.0, nAdds = 0.0, g;
	g = maxGeoDist/*max(maxGeoDist, mesh2->maxGeoDist)*/, unitGeoScaling(g), mesh2->unitGeoScaling(g); //to get the gtDistortion (mean geodesic error) in [0, 1]
	size_t n = samples.size(), nMatches = 0;	
	for (size_t i = 0; i < n; i++)
		if (verts[ samples[i] ]->gtIdx != -1 && verts[ samples[i] ]->matchIdx != -1)
		{
			int v2 = verts[ samples[i] ]->matchIdx; //definitely a sample whose geodesics to all other verts are known for sure
			verts[ samples[i] ]->gtDistortion = (mesh2->verts[v2]->geodesic[ verts[ samples[i] ]->gtIdx ] / 1.0);//mesh2->avgEdgeLen);			
			gtDistortion += verts[ samples[i] ]->gtDistortion;
			nAdds++;
		}
	cout << "gt-distortion of this loaded map: " << gtDistortion / nAdds << " (over " << nAdds << " adds)\n";
}

void Mesh::resultFromFile3(Mesh* mesh2, int dl)
{
	//loads matches of the Hybrid Functional Maps, EchoMatch, and Unsupervised Learning of Spectral Shape Matching paper

	char fName[250], fOutName[250];
	if (dl == 1) //deep learning paper # 1: HFM
		sprintf(fName, "hybridfmaps-results-v2\\outputs\\shrec16\\converted_mats\\scape200_decimated_david\\%d-%d_p2p_nocomma.txt", id, mesh2->id), sprintf(fOutName, "hfm%d-%d.csv", id, mesh2->id);
		//sprintf(fName, "toscahorse_quadric_05-toscahorse_quadric_02_p2pnocomma.txt"), sprintf(fOutName, "hfm5-2.csv");
	else if (dl == 2) //deep learning paper # 2: ULRSS
		sprintf(fName, "ulrssm-shrec16\\converted_mats\\ulrssm\\shrec16_null_to_partial\\%d-%d_p2pnocomma.txt", id, mesh2->id), sprintf(fOutName, "ulrss%d-%d.csv", id, mesh2->id);
	else              //deep learning paper # 3: EchoMatch
		sprintf(fName, "EchoMatchResults\\shrec16\\echo_match_shrec16_benchmark_dino50_\\predicted_maps\\david%d_david%d_partial3.map", id, mesh2->id), sprintf(fOutName, "echomatch%d-%d.csv", id, mesh2->id);
	FILE* fPtr;
	if (! (fPtr = fopen(fName, "r")))
	{
		cout << "WARNING: not found: " << fName << "\n\n";
		return; //exit(0); don't exit 'cos gt evaluation is optional
	}
	cout << "freading map from " << fName << "..\n";
	int v1, v2 = 0;
	while (fscanf(fPtr, "%d\n", &v1) != EOF) //go until the end of file
		mesh2->verts[(int) v2]->matchIdx = (int) v1 - 1, verts[(int) v1 - 1]->matchIdx = (int) v2++; //dual operation and v2 increment; -1 'cos hfmaps output indices are 1-based
	fclose(fPtr);
	///////////////// manyToOneMesh1v2() to improve coloring on mesh1 enable manyTo1 map which matches the unmatched vertices of mesh1 via proximity (disable to see their original result w/ many unmatched vertices) /////////////////
	for (size_t v1 = 0; v1 < verts.size(); v1++) //for each unmatched mesh1 vert v1, find the closest v2 below
    {
        if (verts[v1]->matchIdx != -1) //match only the unmatched verts (this already knows its match)
            continue;
        double minDistDiff = INF;
        for (size_t vm = 0; vm < verts.size(); vm++)
            if (verts[vm]->matchIdx != -1)
            {
                double currDistDiff = distanceBetween2(verts[v1]->coords, verts[vm]->coords);
                if (currDistDiff < minDistDiff)// && angleBetween(verts[v1]->unormal, verts[v2]->unormal) < normalCompatibilityAngle) //2nd condition ensures that normals are compatible b/w this and mesh2 vert normals
                {
                    minDistDiff = currDistDiff;
                    verts[v1]->matchIdx = verts[vm]->matchIdx;
                }
            }
    }
	///////////////// manyToOneMesh1v2() to improve coloring on mesh1 enable manyTo1 map which matches the unmatched vertices of mesh1 via proximity (disable to see their original result w/ many unmatched vertices) ends /////////////////*/
	
	///////////////// adapted from transferColors() /////////////////	
	for (size_t j = 0; j < mesh2->verts.size(); j++) //geodesic distances to the reference vertex for visualization only (to transfer smooth colors from mesh2 to mesh1)
	{
		size_t v = 3643;//0; //reference geo/seo source vertex
		v = (v < mesh2->verts.size() ? v : 0); //stay in bounds
		mesh2->dijkstraShortestPaths(v); //geodesic[] ready (returns immediately if v.geodesic is already not empty)
		mesh2->verts[j]->geoToRef = mesh2->verts[v]->geodesic[j]; //make sure vertColors() below uses geoToRef too
	}
	mesh2->vertColors(0, false); //predefined mesh2 colors are set
	bool print1 = true;
	int nMatched = 0;
	for (size_t i = 0; i < verts.size(); i++) //all or some (some in partial matching case) mesh1 vert has a match and it'll use that matchIdx to learn its color from the fully colored mesh2
	{
		verts[i]->color = new float[3]; //alloc memo for v.color and init it to mesh color (gray)
		verts[i]->color[0] = color[0]; //unmatched mesh2 verts will have this default color (gray set in Mesh's constructor)
		verts[i]->color[1] = color[1];
		verts[i]->color[2] = color[2];
		if (verts[i]->matchIdx != -1)
			verts[i]->color = mesh2->verts[ verts[i]->matchIdx ]->color, nMatched++;
		else if (print1)
			cout << "WARNING: unmatched mesh1 vert (there may be other unmatched ones; printing stopped)\n\n", print1 = false;
	}
	cout << nMatched << " colors pulled back from mesh2\n";
	///////////////// adapted from transferColors() ends /////////////////

	cout << "fprinting map to " << fOutName << "..\n";
	FILE* fPtrCsv = fopen(fOutName, "w");
	//mesh1 to mesh2
	for (size_t i = 0; i < verts.size(); i++)
		if (i == verts.size() - 1) //no comma after the last match
			fprintf(fPtrCsv, "%d\n", (verts[i]->matchIdx != -1 ? verts[i]->matchIdx : 0));
		else //comma after each match
			fprintf(fPtrCsv, "%d,", (verts[i]->matchIdx != -1 ? verts[i]->matchIdx : 0));
	fclose(fPtrCsv);
}

void Mesh::addVertexND(double* c)
{
	//quickly add verts w/ coords c to mesh without duplication check
	verts.push_back(new Vertex(c));
}

void Mesh::addTriangle(int v1i, int v2i, int v3i)
{
	//adds tri v1i-v2i-v3i to the triInd'th element of tris

	tris.push_back( new Triangle(v1i, v2i, v3i) );

	//update each vertex with its new tri and vert neighbors
	verts[v1i]->triNeighbors.push_back( tris.size() - 1 );
	//addVertNeighbor returns true if verts[v1i] is not a neighbor of [v2i]; so they can form a new edge
	if (verts[v1i]->addVertNeighbor(v2i))
		addEdge(v1i, v2i);
	if (verts[v1i]->addVertNeighbor(v3i))
		addEdge(v1i, v3i);

	verts[v2i]->triNeighbors.push_back( tris.size() - 1 );
	verts[v2i]->addVertNeighbor(v1i); //v1i v2i edge is (probably) added above; do not even try to add v2i v1i edge here
	if (verts[v2i]->addVertNeighbor(v3i))
		addEdge(v2i, v3i);	//v3i v2i is not added above, so add this if ok

	verts[v3i]->triNeighbors.push_back( tris.size() - 1 );
	verts[v3i]->addVertNeighbor(v1i);	//no addEdge risk, so no if control
	verts[v3i]->addVertNeighbor(v2i);
}

int Mesh::borderVerts()
{
	//computes border vertices, i.e., endpoints of the border edges (that are incident to 1 triangle only)

	for (size_t i = 0; i < verts.size(); i++)
		verts[i]->border = false;
	for (size_t i = 0; i < edges.size(); i++)
	{
		int* shrdTris = new int[50]; //int[2] is expected but for nonmanifold meshes, e.g., edge incident to 4 triangles, [2] causes segmentation fault so be generous
		trianglesSharedBy(edges[i], shrdTris);
		if (shrdTris[1] == -1)
			verts[ edges[i]->v1i ]->border = verts[ edges[i]->v2i ]->border = true;
		delete [] shrdTris;
	}
	int nBorders = 0;
	for (size_t i = 0; i < verts.size(); i++)
		nBorders += verts[i]->border;
	return nBorders;
}

void Mesh::supportEdges()
{
	//generates support edges that improve dijkstra accuracy by letting it passing through faces; for good dijkstra, edges should be isotropic, i.e. uniform in all directions, which is achieved by these support edges

	//given a triangle t and its 3 neighbor triangles, if they are ~coplanar, then support pair t.v1 and tn1.vo1, t.v2 and tn2.vo2, t.v3 and tn3.vo3, i.e. all valid vertex pairs where vo means other vertex not included in t and tn is neighbor triangle
	int* shrdTris = new int[20]; //int[2] is expected but for nonmanifold meshes, e.g., edge incident to 4 triangles, [2] causes segmentation fault so be generous
	for (size_t t = 0; t < tris.size(); t++)
	{
		//get the 3 triangles adjacent to t; if t + these three are almost coplanar, then add the support edges
		int v1i = tris[t]->v1i, v2i = tris[t]->v2i, v3i = tris[t]->v3i, t1 = -1, t2 = -1, t3 = -1,
			e1 = edgeSharedBy(v1i, v3i), e2 = edgeSharedBy(v1i, v2i), e3 = edgeSharedBy(v2i, v3i);
		trianglesSharedBy(edges[e1], shrdTris); //call-by-ref output to shrdTris
		if (shrdTris[1] > -1) //e1 has exactly 2 triangle neighbors so t1 is defined properly
			t1 = (shrdTris[0] == t ? shrdTris[1] : shrdTris[0]);
		trianglesSharedBy(edges[e2], shrdTris); //call-by-ref output to shrdTris
		if (shrdTris[1] > -1)
			t2 = (shrdTris[0] == t ? shrdTris[1] : shrdTris[0]);
		trianglesSharedBy(edges[e3], shrdTris); //call-by-ref output to shrdTris
		if (shrdTris[1] > -1)
			t3 = (shrdTris[0] == t ? shrdTris[1] : shrdTris[0]);	

		//coplanarity tests
		if (t1 != -1 && dotProduct(tris[t]->unormal, tris[t1]->unormal) < 0.95) //dot product = ~1 means that t & t1 are ~coplanar
			continue;
		if (t2 != -1 && dotProduct(tris[t]->unormal, tris[t2]->unormal) < 0.95) //dot product = ~1 means that t & t2 are ~coplanar
			continue;
		if (t3 != -1 && dotProduct(tris[t]->unormal, tris[t3]->unormal) < 0.95) //dot product = ~1 means that t & t3 are ~coplanar
			continue;

		//decide other vertices, i.e., the neighboring ones not in t
		int u1 = -1, u2 = -1, u3 = -1;
		if (t1 != -1) //u1 is the other vert of the 1st adj tri (t1)
		{
			if (! (tris[t1]->v1i == v1i || tris[t1]->v1i == v3i)) //v2i ignored 'cos t1 is created via v1i-v3i edge
				u1 = tris[t1]->v1i;
			else if (! (tris[t1]->v2i == v1i || tris[t1]->v2i == v3i))
				u1 = tris[t1]->v2i;
			else
				u1 = tris[t1]->v3i;
		}
		if (t3 != -1) //u2 is the other vert of the 3rd adj tri (t3)
		{
			if (! (tris[t3]->v1i == v3i || tris[t3]->v1i == v2i))
				u2 = tris[t3]->v1i;
			else if (! (tris[t3]->v2i == v3i || tris[t3]->v2i == v2i))
				u2 = tris[t3]->v2i;
			else
				u2 = tris[t3]->v3i;
		}
		if (t2 != -1) //u3 is the other vert of the 2nd adj tri (t2)
		{
			if (! (tris[t2]->v1i == v2i || tris[t2]->v1i == v1i))
				u3 = tris[t2]->v1i;
			else if (! (tris[t2]->v2i == v2i || tris[t2]->v2i == v1i))
				u3 = tris[t2]->v2i;
			else
				u3 = tris[t2]->v3i;
		}

		//add all possible supports if not added already (that test in addSupportEdge())
		addSupportEdge(v1i, u2);
		addSupportEdge(v2i, u1);
		addSupportEdge(v3i, u3);
		addSupportEdge(u1, u3);
		addSupportEdge(u1, u2);
		addSupportEdge(u2, u3);
	}
	delete [] shrdTris;
	if (printHeavy)
		cout << sedges.size() << " support edges generated\n";
}

void Mesh::addSupportEdge(int v1i, int v2i)
{
	//adds the support edge b/w v1i & v2i w/ length len if it has not been previously generated

	if (v1i == -1 || v2i == -1)
		return; //may validly be -1 if u1/2/3 is -1 in the caller
	//check if virtual support edge connection already exists
	for (size_t se = 0; se < verts[v1i]->sedgeList.size(); se++)
	{
		Edge* sedge = sedges[ verts[v1i]->sedgeList[se] ];
		int otherVert = (sedge->v1i == v1i ? sedge->v2i : sedge->v1i);
		if (otherVert == v2i) //virtual connection already exists
			return;
	}
	//add the support edge for the first time b/w these 2 verts
	verts[v1i]->sedgeList.push_back( sedges.size() );
	verts[v2i]->sedgeList.push_back( sedges.size() ); //do same thing for [v2i]
	sedges.push_back(new Edge(v1i, v2i, distanceBetween(verts[v1i]->coords, verts[v2i]->coords)));
}

int Mesh::edgeSharedBy(int vi, int vj)
{
	//returns idx of edge[vi, vj]

	for (size_t vie = 0; vie < verts[vi]->edgeList.size(); vie++)
		for (size_t vje = 0; vje < verts[vj]->edgeList.size(); vje++)
			if (verts[vj]->edgeList[vje] == verts[vi]->edgeList[vie])
				return verts[vj]->edgeList[vje];

	cout << "WARNING: missing edge b/w " << vi << " and " << vj << "\n\n";

	//edge not exist
	return -1;
}

void Mesh::addEdge(int v1i, int v2i)
{
	//adds edge v1i-v2i to the edgeInd'th element of edges

	Edge* newcomer = new Edge(v1i, v2i, distanceBetween(verts[v1i]->coords, verts[v2i]->coords));

	//update the list of edges of which verts[v1i] is a member
	verts[v1i]->edgeList.push_back( edges.size() );
	verts[v2i]->edgeList.push_back( edges.size() ); //do same thing for [v2i]
	edges.push_back(newcomer);

	//information
	if (newcomer->length < minEdgeLen)
		minEdgeLen = newcomer->length;
	if (newcomer->length > maxEdgeLen)
		maxEdgeLen = newcomer->length;
	edgeLenTotal += newcomer->length;
}

void Mesh::trianglesSharedBy(Edge* e, int* shTris) //result written to shTris via call-by-ref
{
	//tri t is shared by bigEdge iff it exists in the triNeigbors of verts[v1i] & [v2i]

	int v1i = e->v1i, v2i = e->v2i, a = 0, tr1, tr2;
	for (size_t t1 = 0; t1 < verts[v1i]->triNeighbors.size(); t1++)
	{
		tr1 = verts[v1i]->triNeighbors[t1];
		for (size_t t2 = 0; t2 < verts[v2i]->triNeighbors.size(); t2++)
		{
			tr2 = verts[v2i]->triNeighbors[t2];
			if  (tris[tr2] != 0 && tr1 == tr2)
				shTris[a++] = tr1; //or tr2
		}
	}
	//if e is a border edge (cannot happen in watertight meshes), take precautions by using -1 flag
	if (a == 1)
		shTris[1] = -1; //border indicator
	else if (a == 0)
	{
		shTris[0] = shTris[1] = -3;
		cout << "WARNING: an edge w/ no triangle neighbor; must not happen\n\n";
	}
	else if (a > 2)
		shTris[1] = -2; //non-manifold indicator
}

void Mesh::scaleEdges(double scaleFactor)
{
	//scales edge lengths, e.g., after unitAreaScaling()

	edgeLenTotal = 0.0;
	minEdgeLen = INF;
	maxEdgeLen = -INF;
	for (size_t e = 0; e < edges.size(); e++)
	{
		//edges[e]->length = distanceBetween(verts[ edges[e]->v1i ]->coords, verts[ edges[e]->v2i ]->coords);
		edges[e]->length *= scaleFactor;
		if (edges[e]->length < minEdgeLen)
			minEdgeLen = edges[e]->length;
		if (edges[e]->length > maxEdgeLen)
			maxEdgeLen = edges[e]->length;
		edgeLenTotal += edges[e]->length;
	}
	avgEdgeLen = edgeLenTotal / edges.size();
	for (size_t e = 0; e < sedges.size(); e++)
		sedges[e]->length *= scaleFactor;
}

bool printOnce0 = true;
void Mesh::unitGeoScaling(double scaleFactor, bool updateMaxGeo)
{
	//scales only the geodesic[] entries to [0 ,1] (tanh/cubic function that weights a given map distortion needs this normalization for its [0, 1] input requirement)

	if (scaleFactor == 1.0 && printOnce0)
		cout << "WARNING: normalizer=1 makes normalization ineffective; is this the intent? (there may be others; printing stopped)\n\n", printOnce0 = false;
	double maxG = -INF;
	for (size_t i = 0; i < verts.size(); i++)
		if (! verts[i]->geodesic.empty() && ! verts[i]->geoNormalized) //don't normalize already-normalized verts (may happen only in the evaluation part of denseMapLargeN())
		{
			for (size_t j = 0; j < verts.size(); j++)
				verts[i]->geodesic[j] /= scaleFactor, maxG = (verts[j]->sample && verts[i]->geodesic[j] > maxG ? verts[i]->geodesic[j] : maxG); //scaleFactor = max(mesh1.maxGeoDist, mesh2.maxGeoDist)
			verts[i]->geoNormalized = true; //currently normalized
		}
	//radius /= scaleFactor; Correspondence constructor calls here so keep radius intact by not doing this
	if (updateMaxGeo)
		maxGeoDist = maxG, cout << "max geo dist after geo normalization: " << maxGeoDist << endl;
}

void Mesh::unitAreaScaling()
{
	//scales this mesh such that its total surface area is 1.0 (eigen stuff, e.g., hks, is supposed to work better w/ this scaling)

	double scaleFactor = sqrt(totalArea), totAreaBef = totalArea;
	for (int i = 0; i < (int) verts.size(); i++)
		for (int c = 0; c < 3; c++)
			verts[i]->coords[c] /= scaleFactor;
	//boundingBox(), centerOfMass();
	normalsAndAreas(); //t.area totalArea and so on updated	
	scaleEdges(1.0 / scaleFactor); //so that upcoming geodesics use the scaled edge lengths	
	cout << "total mesh area before/after scaling: " << totAreaBef << "/" << totalArea << "\n\n";
}

void Mesh::normalsAndAreas()
{
	//sets normals and areas of tris and verts in mesh

	//triangles
	totalArea = 0.0;
	bool printOnce = true;
	for (size_t t = 0; t < tris.size(); t++)
	{
		int v1i = tris[t]->v1i, v2i = tris[t]->v2i, v3i = tris[t]->v3i;
		bool printed = tris[t]->setNormal(verts[v1i]->coords, verts[v2i]->coords, verts[v3i]->coords, printOnce);
		tris[t]->setArea(verts[v1i]->coords, verts[v2i]->coords, verts[v3i]->coords);
		totalArea += tris[t]->area;
		if (printed)
			printOnce = false;
	}
	//vertices
	bool printOnce1 = true, areaWeightedNormal = true;
	double normalTotal[3], nl1[3], nl2[3];
	int weightType = 3; //0: uniform, 1: area, 2: angle, 3: area*angle
	for (size_t v = 0; v < verts.size(); v++)
	{
		normalTotal[0] = normalTotal[1] = normalTotal[2] = 0.0;
		verts[v]->area = 0.0; //init patchGC in case patchAndDescriptors() not called, i.e., no path sampling (for printing only)
		double howManyTriNeighb = 0.0;
		for (size_t tneighb = 0; tneighb < verts[v]->triNeighbors.size(); tneighb++)		
			if (tris[ verts[v]->triNeighbors[tneighb] ]->area >= DBL_MIN)
			{
				int t = verts[v]->triNeighbors[tneighb];
				double area = tris[t]->area;
				for (int i = 0; i < 3; i++)
					normalTotal[i] += (areaWeightedNormal ? area : 1.0) * tris[t]->normal[i];
				verts[v]->area += area; //local area around v is the sum of areas of adjacent triangles (1-ring)
				howManyTriNeighb++;
				if (weightType < 2) //uniform or area weights
					for (size_t i = 0; i < 3; i++)
						normalTotal[i] += (weightType == 1 ? area : 1.0) * tris[t]->unormal[i];
				else				//angle-based weights
				{
					double angle;
					int v1 = (tris[t]->v1i == v ? tris[t]->v2i : tris[t]->v1i), v2 = tris[t]->v1i;
					if (v2 == v1 || v2 == v)
						v2 = ((tris[t]->v2i != v1 && tris[t]->v2i != v) ? tris[t]->v2i : tris[t]->v3i);
					for (size_t c = 0; c < 3; c++)
					{
						nl1[c] = verts[v1]->coords[c] - verts[v]->coords[c]; //from v to neighbor vert v1
						nl2[c] = verts[v2]->coords[c] - verts[v]->coords[c]; //from v to neighbor vert v2
					}
					angle = angleFromCotan(nl1, nl2); //angle of t at vertex v (radian or degree does not matter 'cos we normalize the end result hence only the ratios contribute as weights)
					for (size_t i = 0; i < 3; i++)
						normalTotal[i] += (weightType == 2 ? angle : area*angle) * tris[t]->unormal[i];
				}
			}
			else if (printOnce1)
			{
				cout << "WARNING: there exists a tiny triangle " << verts[v]->triNeighbors[tneighb] << " w/ area: " << tris[ verts[v]->triNeighbors[tneighb] ]->area << " (there may be others; printing stopped)\n\n"; //identical vertex indices in the same face may cause this
				printOnce1 = false;
			}
		if (howManyTriNeighb != 0.0) //if 0, go on w/ undefined normal (should never happen)
		{
			for (size_t i = 0; i < 3; i++)
				normalTotal[i] /= howManyTriNeighb;
			verts[v]->setNormal(normalTotal);
		}
		else
		{
			verts[v]->setNormal(normalTotal); //normalTotal=0 here so warning will be printed during setNormal(); setNormal() defines/mallocs the normal which is necessary so this call is mandatory
			cout << "WARNING: undefined normal for vertex " << v << " (relax area threshold)\n\n"; //happens if there exist isolated vertices w/o any connection to the mesh (howManyTriNeighb=0)
		}
	}
}

void Mesh::additionalRotations(int axis, double degree)
{
	//applies additional x- y- and/or z-axis rotations to mesh for better visualization (makes more sense when 2 meshes on screen)

//	cout << "additional rotations for visual convenience..\n";
	double x = 0, y = 0, z = 0, theta = degree * (PI / 180.0); //radian rotation angle	
	if (axis == 2 || axis == 3) //3: z-rot followed by y-rot
		//z-axis rotation by theta degree
		for (int i = 0; i < (int) verts.size(); i++)
		{
			double x = verts[i]->coords[0]*cos(theta) - verts[i]->coords[1]*sin(theta);
			double y = verts[i]->coords[0]*sin(theta) + verts[i]->coords[1]*cos(theta);
			verts[i]->coords[0] = x;
			verts[i]->coords[1] = y;
		}
	if (axis == 1 || axis == 3)
		//y-axis rotation by theta degree
		for (int i = 0; i < (int) verts.size(); i++)
		{
			x = verts[i]->coords[2]*sin(theta) + verts[i]->coords[0]*cos(theta);
			z = verts[i]->coords[2]*cos(theta) - verts[i]->coords[0]*sin(theta);
			verts[i]->coords[0] = x;
			verts[i]->coords[2] = z;
		}
	if (axis == 0)
		//x-axis rotation by theta degree
		for (int i = 0; i < (int) verts.size(); i++)
		{
			y = verts[i]->coords[1]*cos(theta) - verts[i]->coords[2]*sin(theta);
			z = verts[i]->coords[1]*sin(theta) + verts[i]->coords[2]*cos(theta);
			verts[i]->coords[1] = y;
			verts[i]->coords[2] = z;
		}
}

void Mesh::dijkstraShortestPaths(int sourceVert)
{
	//computes shortest paths (geodesics) from sourceVert to all other verts
	
	if (! verts[sourceVert]->geodesic.empty()) //don't recompute if already computed (if(verts[sourceVert]->geodesic) would have been used if double* geodesic; was in use)
		return;
	//verts[sourceVert]->geodesic = new double[ verts.size() ]; //length of geodesic from sourceVert to any mesh vertex
	verts[sourceVert]->geodesic.resize( verts.size() ); //preallocations to save tiny time (requesting more space than needed is also OK but i request the exact amount now); resize() is less efficient as it default-constructs/initializes all n elements (to 0 here)
														//but it's mandatory 'cos push_back's required after reserve() to increase the size: writing to uninitialized memory via [] after reserve() can corrupt memory, crash your program, or appear to work silently while hiding serious bugs (so reserve must work with push_backs)
#ifdef USE_FIB_HEAP
	//use Fibonacci heap to get O(VlgV + E); minHeap takes O(VlgV + ElgV), simple array takes O(V^2)
	FibHeap* h = new FibHeap(); //fibheap 0.008secs vs array 0.8secs on an 8K-vertex component; 0.001 vs. 0.03secs for 268-vertex component; so fibheap is preferred
#endif
	//initialization
	for (size_t v = 0; v < verts.size(); v++)
	{
		verts[v]->dSpanned = (v == sourceVert ? 0.0 : INF);
#ifdef USE_FIB_HEAP
		FibHeapNode* hn = new FibHeapNode(verts[v]->dSpanned, v); //v is a handle from headNode to verts,
		verts[v]->heapNode = hn;								  //and hn is an handle from verts to heapNode
		h->Insert(hn);
#else
		verts[v]->tmp = 0; //1 means shortest path distance is finalized
#endif
		verts[v]->prev = -1; //path info via prev values (sourceVert.prev remains -1 (0-length path to itself))	
	}
	while (true)
	{
		//extract vert w/ min dSpanned value
		//recall that Dijkstra overall cost: O(extractMin) + O(aggregate relax) = O(VlgV) + O(E) w/ FibHeap, O(VlgV) + O(ElgV) w/ minHeap, O(V^2) w/ array; sparse graphs (our meshes) E=O(V), dense graphs E = O(V^2)
#ifdef USE_FIB_HEAP
		FibHeapNode* minNode = h->ExtractMin(); //ExtractMin O(lgV) for both minHeap and FibHeap but DecreaseKey amortized cost is O(1) in FibHeap (no such feature in minHeap, always O(lgV) DecreaseKey)
		if (minNode == NULL)
			break; //no more verts w/ unknown shortest dist; so, we're done
		int minDidx = minNode->element, e, va;
		double minD = minNode->key;
#else
		int minDidx = -1, e, va;
		double minD = INF;
		for (size_t v = 0; v < verts.size(); v++)
			if (verts[v]->tmp == 0 && verts[v]->dSpanned < minD)
				minD = verts[v]->dSpanned, minDidx = v;
		if (minDidx == -1)
			break;
		verts[minDidx]->tmp = 1;
#endif
		verts[sourceVert]->geodesic[minDidx] = minD; //original unnormalized value
		delete minNode; //new hn above is deleted here to save memo
		//minD is an ascending value; so, when i reach far enough verts, i can safely skip the rest since all verts visited so far have finalized their shortest paths,
		//which means any vertex that i visit after this pnt will be have dist >= minD to sourceVert
		//if (minD > far)  //far enough vert reached
		//	break;

		//relax each edge incident to minDidx
		for (size_t ve = 0; ve < verts[minDidx]->edgeList.size(); ve++)
		{
			e = verts[minDidx]->edgeList[ve];
			va = (edges[e]->v1i == minDidx ? edges[e]->v2i : edges[e]->v1i);
			if (verts[minDidx]->dSpanned + edges[e]->length < verts[va]->dSpanned)
			{
				verts[va]->dSpanned = verts[minDidx]->dSpanned + edges[e]->length; //relaxation (for length)
#ifdef USE_FIB_HEAP
				h->DecreaseKey(verts[va]->heapNode, verts[va]->dSpanned); //worst O(logV), amortized O(1)
#endif
				verts[va]->prev = minDidx; //relaxation (for path)
			}
		}
		//relax each support edge, if any, incident to minDidx
		for (size_t ve = 0; ve < verts[minDidx]->sedgeList.size(); ve++)
		{
			e = verts[minDidx]->sedgeList[ve];
			va = (sedges[e]->v1i == minDidx ? sedges[e]->v2i : sedges[e]->v1i);
			if (verts[minDidx]->dSpanned + sedges[e]->length < verts[va]->dSpanned)
			{
				verts[va]->dSpanned = verts[minDidx]->dSpanned + sedges[e]->length; //relaxation (for length)
#ifdef USE_FIB_HEAP
				h->DecreaseKey(verts[va]->heapNode, verts[va]->dSpanned); //worst: O(logV), amortized O(1)
#endif
				verts[va]->prev = minDidx; //relaxation (for path)
			}
		}
	} //end of while (true)
#ifdef USE_FIB_HEAP
	delete h;
#endif
	verts[sourceVert]->geoNormalized = false; //currently unnormalized
	//store paths from this sourceVert to each target vertex to avoid dijkstra recomputations later (prev values must be updated for each new sourceVert so doing this later causes dijkstra recomputations)
//	pathsFromVertex(sourceVert); //based on the current set of v.prev values (this causes overflow memory for n=250 samples on 52K-vertex mesh so defer it until i find all n=250 samples and then fill it from samples to samples only, not samples to vertices which causes memory overflow)
}

bool comparator(pair<double, int> a, pair<double, int> b) { return a.first > b.first; } //by default sort() sorts the vector of pairs in ascending order of the 1st member of pairs; provide this function for descending order
void Mesh::FPS(size_t N, bool curvatureBased)
{
	//samples N evenly-spaced FPS (Farthest Point Sampling) points based on the paper The Farthest Point Strategy for Progressive Image Sampling (default curvatureBased=false)
	//curvatureBased=true is same default FPS except instead of taking the farthest point we first create a pool of 1% farthest points and from this pool we select the one w/ max curvature

	cout << (curvatureBased ? "Curv. e" : "E") << "venly-spaced sampling [" << N << "]..";
	clock_t t = clock();
	int nGTs = 0; //number of vertices with valid ground-truth match labels
	for (size_t v = 0; v < verts.size(); v++)
	{
/*		if (verts[v]->geodesic) //in case already deleted
			delete [] verts[v]->geodesic; //arrays allocated with new[] must be deallocated with delete[]
		verts[v]->geodesic = NULL; still caused memory leaks (due to shallow copy etc) so switched to dynamic array vector*/
		verts[v]->geodesic.clear(), verts[v]->sample = verts[v]->marked = false; //marked=true means v is in the vicinity of an existing FPS sample (for replaceSample() only; FPS() is already evenly-spaced (curvatureBased slightly distracts it but anyway) so this marked is for replaceSample() that has no evenly-spacing respect)
		nGTs += (verts[v]->gtIdx != -1);
	}
	fpsFromGT = (nGTs < verts.size()*0.05 ? false : fpsFromGT); //too few gt-equipped vertices (less than 5%) so switch to typical fpsFromGT=false mode
	samples.clear();

	vector< pair<double, int> > pool;
	//fixed number (N) of fps (farthest point sampling) samples goes to samples[] first of which is an extremity vertex
	int nReplacements = 0, extVert = extremeVert(curvatureBased); //as i want the 1st ones to be compatible, e.g. extremities, use extVert instead of a random startVert (might pick the max curvature vertex here in curvatureBased mode but selecting the extremity is also a good idea)
	//verts[extVert]->sample = true; samples.push_back(extVert); dijkstraShortestPaths(extVert); already made inside extremeVert()
	//fill the remaining N-1 items via fps for uniform distribution that lets nice coverage on models; before that add locMaxGC vertices (low-curvature ones and close ones trimmed by locMaxGCs()) to samples
//	locMaxSampling(N, false); not needed initially 'cos main().qapMap() works on very low sample resolutions, e.g., 10 vs 3, 10 vs 5, 10 vs 10 (complete matching); locMaxSampling() currently processes only the ROI vertices so this would be ineffective anyway
	//now fill the remaining items via fps for uniform distribution that lets nice coverage on models
	while (samples.size() < N)
	{
		//in fps, the next sample candidate x_i \in all vertices finds its closest existing sample x_clo and remembers d_i = d(x_i, x_clo); if d_i is the max among all other candidates closest remembered dists, then x_i is selected as the next sample;
		//here my candidates are restricted to nothing, i.e. i use all verts as candidates; so, find the closest existing sample e2 \in samples to s \in verts and add s to samples as nextSample if that closest dist is very big, i.e. the biggest among all closest dists for different s's
		double maxDist = -INF;
		int nextSample = -1;
		pool.clear();
		for (size_t s = 0; s < verts.size(); s++) //fps selects from amongst all vertices
		{
			if (verts[s]->edgeList.empty() || //isolated verts should never be samples (even if they'd have dSapnning=INF values); they occur on microholes/holes/viewing classes of shrec11 dataset
				verts[s]->sample) //pure FPS not need this but coming here from another sampler may insert duplicates that are already made sample by that sampler
				continue;
			if (fpsFromGT && verts[s]->gtIdx == -1) //for fair evaluation select amongs gt-equipped vertices
				continue;
			int e2 = -1;
			double minDist = INF;
			for (size_t e = 0; e < samples.size(); e++) //find the existing sample e \in samples that is closest to s \in verts
				if (verts[ samples[e] ]->geodesic[s] <= minDist)
				{
					minDist = verts[ samples[e] ]->geodesic[s];
					e2 = samples[e];
				}
			//instead of just taking the farthest point, we first create a pool of the top 1% farthest points, and from this pool we then select the one with highest curvature
			if (curvatureBased) //instead of making nextSample=s if verts[e2]->geodesic[s] is max (else if part below), insert verts[e2]->geodesic[s] value to the pool and later use the 1% of this pool (to pick the max curvature item as the nextSample)
				pool.push_back({verts[e2]->geodesic[s], s});
			else if (verts[e2]->geodesic[s] > maxDist)
			{
				maxDist = verts[e2]->geodesic[s];
				nextSample = s;
			}
		} //end of for s
		if (curvatureBased && ! pool.empty())
		{
			//now sort the pool (in descending order) to get the top 1% farthest points
			sort(pool.begin(), pool.end(), comparator); //if comparator not provided then sorts in ascending order
			double maxGC = -INF; //nextSample will be the max curvature item from the top 1% of the pool
			for (size_t i = 0; i < pool.size()*0.01; i++) //too large pool hurts evenly spacing, e.g., puts samples to all fingers makign small hand region excessively dense sampled
				//if (abs(verts[ pool[i].second ]->gc) > maxGC) //still meaningful as very negative gc means deep/radical saddle points which are also important
				if (verts[ pool[i].second ]->gc > maxGC) //locMaxGC=true made based on non-abs values so be consistent
					//maxGC = abs(verts[ pool[i].second ]->gc), nextSample = pool[i].second; use this if abs in if is enabled
					maxGC = verts[ pool[i].second ]->gc, nextSample = pool[i].second;
			if (nextSample == -1 || verts[nextSample]->border)//&& ! verts[ pool[0].second ]->border) //border.gc not reliable so replace it w/ the farthest point, i.e., no curvatureBased decision, at index 0 (descending order)
				nextSample = pool[0].second; //pool[0] allowed to be border (2nd condition disabled) but if it was not allowed then i'd have checked its 1-ring for a non-border
		}
		if (nextSample == -1)
		{
			cout << "WARNING: FPS could not find a valid sample; quiting early (pool[] may be empty)\n\n";
			break;
		}
		verts[nextSample]->sample = true;
		samples.push_back(nextSample);
		dijkstraShortestPaths(nextSample); //upcoming iteration needs samples.geodesic
		if (curvatureBased)
		{
			nReplacements += replaceSample(nextSample); //replace nextSample with an appropriate max-curvature locMaxGC vertex
			for (size_t v = 0; v < verts.size(); v++)
				if (verts[ samples[ samples.size() - 1 ] ]->geodesic[v] < radius) //last sample (either nextSample or its replacement that is local to replaceSample())
					verts[v]->marked = true;
		}
	} //end of while N
/*	//shuffle sample points (makes qucoop more interesting, e.g., w/o shuffling q_free is ~0 requiring ~no changes in identity permutation matrix)
	if (id > 0 && samples.size() > 4)
	{
		//swap
		int tmp = samples[1];
		samples[1] = samples[4];
		samples[4] = tmp;
		//swap
		tmp = samples[0];
		samples[0] = samples[2];
		samples[2] = tmp;
	}//*/
/*	if (id == 0 && samples.size() > 4)
	{
		//swap
		int tmp = samples[1];
		samples[1] = samples[8];
		samples[8] = tmp;
		//swap
		tmp = samples[0];
		samples[0] = samples[6];
		samples[6] = tmp;
	}//*/
	//set maxGeoDist properly which in turn gives me a chance to normalize geodesics if necessary; stores path vertices b/w pairs of samples too
	maxGeoDist = -INF; //maxGeoDist can be from sample to non-sample vertex, so after normalization geodesic b/w most-separated sample pair may remain < 1 (not the case anymore; see j-loop below)
	cout << "path verts..";
	for (size_t i = 0; i < samples.size(); i++)
	{
		for (size_t j = 0; j < samples.size(); j++)//for (size_t v = 0; v < verts.size(); v++) matching samples to samples so let maxGeoDist b/w sample to sample (curvatureBased=false selects samples[0] and [1] so v-loop also ok but =true may pick non-sample w/ v-loop)
			if (verts[ samples[i] ]->geodesic[ samples[j] ] > maxGeoDist)
				maxGeoDist = verts[ samples[i] ]->geodesic[ samples[j] ];//, max1=samples[i], max2=v;
		//store paths from samples[i] to each sample (to each vertex causes memory overflow) to avoid dijkstra recomputations later (prev values must be updated for each new sourceVert so doing this later causes dijkstra recomputations)
		pathsFromVertex(samples[i]);
	}
	if (nReplacements > 0)
		cout << "max geo dist: " << maxGeoDist << " (nReplace: " << nReplacements << ") in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs\n";//" from " << max1 << " to " << max2 << endl;	
	else
		cout << "max geo dist: " << maxGeoDist << " in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs\n";
	if (printHeavy)
	{
		for (size_t i = 0; i < samples.size(); i++)
		{
			if (i < 6 || (int) i >= samples.size() - 4) //print first and last sample ids
				cout << samples[i] << " ";
			else if (i == 6)
				cout << " .. ";
		}
		cout << endl;
	}
}

void Mesh::FPSroi(size_t N, bool curvatureBased)
{
	//same as FPS() except sampling is done on ROI here, i.e., on verts w/ v.roi=true

	cout << (curvatureBased ? "Curv. e" : "E") << "venly-spaced sampling [" << N << "] on ROI..";
	clock_t t = clock();
	size_t offset = (mesh2 ? 1 : 0);
	clearSamples(offset); //marked values set too
	//fill the rest of the samples (clearSamples() keeps the coarseMap[], not coarseMapFull[], samples)
	vector< pair<double, int> > pool;
	N = (N > nROIverts ? nROIverts : N); //N cannot exceed number of ROI vertices 'cos only they are eligible to be samples in this function
	//first add locMaxGC vertices (low-curvature ones and close ones trimmed by locMaxGCs()) to samples
	locMaxSampling(N); //marked values set too
	//now fill the remaining items via fps for uniform distribution that lets nice coverage on models
//	bool eucSamples = false;//(verts.size() > 50000 && N >= 250) || (verts.size() <= 50000 && N >= 1000); //avoid too many dijkstra (20secs euc vs 66secs geo for N=500 on V=52K mesh) by using fast extrinsic euclidean distance during sampling (instead of intrinsic/better geodesic distance)
	int nReplacements = 0;
	while (samples.size() < N)
	{
		//in fps, the next sample candidate x_i \in all vertices finds its closest existing sample x_clo and remembers d_i = d(x_i, x_clo); if d_i is the max among all other candidates closest remembered dists, then x_i is selected as the next sample;
		//here my candidates are restricted to nothing, i.e. i use all verts as candidates; so, find the closest existing sample e2 \in samples to s \in verts and add s to samples as nextSample if that closest dist is very big, i.e. the biggest among all closest dists for different s's
		double maxDist = -INF;
		int nextSample = -1;
		pool.clear();
		for (size_t s = 0; s < verts.size(); s++) //fps selects from amongst all vertices
		{
			if (! verts[s]->roi || //process only the ROI vertices during this sampling
				verts[s]->edgeList.empty() || //isolated verts should never be samples (even if they'd have dSapnning=INF values); they occur on microholes/holes/viewing classes of shrec11 dataset
				verts[s]->sample) //pure FPS not need this but coming here from another sampler may insert duplicates that are already made sample by that sampler
				continue;
			if (fpsFromGT && verts[s]->gtIdx == -1) //for fair evaluation select amongs gt-equipped vertices
				continue;
			int e2 = -1;
			double minDist = INF;
			for (size_t e = 0; e < samples.size(); e++) //find the existing sample e \in samples that is closest to s \in verts
				if (/*! eucSamples &&*/ verts[ samples[e] ]->geodesic[s] <= minDist) //regular geodesic version
					minDist = verts[ samples[e] ]->geodesic[s], e2 = samples[e];
//				else if (eucSamples && distanceBetween2(verts[ samples[e] ]->coords, verts[s]->coords) <= minDist) //fast euclidean version
//					minDist = distanceBetween2(verts[ samples[e] ]->coords, verts[s]->coords), e2 = samples[e];
//			if (! eucSamples)
			{
				//instead of just taking the farthest point, we first create a pool of the top 1% farthest points, and from this pool we then select the one with highest curvature
				if (curvatureBased) //instead of making nextSample=s if verts[e2]->geodesic[s] is max (else if part below), insert verts[e2]->geodesic[s] value to the pool and later use the 1% of this pool (to pick the max curvature item as the nextSample)
					pool.push_back({verts[e2]->geodesic[s], s});
				else if (verts[e2]->geodesic[s] > maxDist)
					maxDist = verts[e2]->geodesic[s], nextSample = s;
			}
/*			else //geodesic not available
			{
				if (curvatureBased) //geodesic above is replaced w/ distanceBetween2()
					pool.push_back({distanceBetween2(verts[e2]->coords, verts[s]->coords), s});
				else if (distanceBetween2(verts[e2]->coords, verts[s]->coords) > maxDist)
					maxDist = distanceBetween2(verts[e2]->coords, verts[s]->coords), nextSample = s;	
			}*/
		} //end of for s
		if (curvatureBased && ! pool.empty())
		{
			//now sort the pool (in descending order) to get the top 1% farthest points
			sort(pool.begin(), pool.end(), comparator); //if comparator not provided then sorts in ascending order
//			if (verts[ pool[0].second ]->roiBorder && ! verts[ pool[0].second ]->border) //not a regular border that may confuse areasRing (borderTreat=1) but an informative roiBorder so select it regardless of its gc
//				nextSample = pool[0].second; not that critical so just make it consistent with FPS()
//			else
			{
				double maxGC = -INF; //nextSample will be the max curvature item from the top 1% of the pool
				for (size_t i = 0; i < pool.size()*0.01; i++) //too large pool hurts evenly spacing, e.g., puts samples to all fingers makign small hand region excessively dense sampled
					//if (abs(verts[ pool[i].second ]->gc) > maxGC) //still meaningful as very negative gc means deep/radical saddle points which are also important
					if (verts[ pool[i].second ]->gc > maxGC) //locMaxGC=true made based on non-abs values so be consistent
						//maxGC = abs(verts[ pool[i].second ]->gc), nextSample = pool[i].second; use this if abs in if is enabled
						maxGC = verts[ pool[i].second ]->gc, nextSample = pool[i].second;
				if (nextSample == -1 || verts[nextSample]->border)//&& ! verts[ pool[0].second ]->border) //border.gc not reliable so replace it w/ the farthest point, i.e., no curvatureBased decision, at index 0 (descending order)
					nextSample = pool[0].second; //pool[0] allowed to be border (2nd condition disabled) but if it was not allowed then i'd have checked its 1-ring for a non-border
			}
		}
		if (nextSample == -1)
		{
			cout << "WARNING: FPS could not find a valid sample; quiting early (pool[] may be empty)\n\n";
			break;
		}
		verts[nextSample]->sample = true;
		samples.push_back(nextSample);
//		if (! eucSamples)
			dijkstraShortestPaths(nextSample); //upcoming iteration needs samples.geodesic
		if (curvatureBased)
		{
			nReplacements += replaceSample(nextSample); //replace nextSample with an appropriate max-curvature locMaxGC vertex
			for (size_t v = 0; v < verts.size(); v++)
				if (verts[ samples[ samples.size() - 1 ] ]->geodesic[v] < radius) //last sample (either nextSample or its replacement that is local to replaceSample())
					verts[v]->marked = true;
		}
//if(id==10)cout << pool.size()*0.01 << " " << radius << " " << nextSample << "\t" << verts[2220]->geodesic[2189] << " // " << verts[2220]->geodesic[nextSample] << "\n"; //2144 --> 2189   -59.6066, 1.23088sssssssssssssssssss
//if (nextSample == 2144) for (int i = 0; i < 110; i++) cout << pool.size() << "\t" << pool[i].first << " &&&&&& " << pool[i].second << "\tgc: " << verts[ pool[i].second ]->gc << " " << verts[ pool[i].second ]->locMaxGC << endl;
	} //end of while N
	//set maxGeoDist properly which in turn gives me a chance to normalize geodesics if necessary; stores path vertices b/w pairs of samples too
	maxGeoDist = -INF; //maxGeoDist can be from sample to non-sample vertex, so after normalization geodesic b/w most-separated sample pair may remain < 1 (not the case anymore; see j-loop below)
	cout << "path verts..";
	for (size_t i = 0; i < samples.size(); i++)
	{
		for (size_t j = 0; j < samples.size(); j++)//for (size_t v = 0; v < verts.size(); v++) matching samples to samples so let maxGeoDist b/w sample to sample (curvatureBased=false selects samples[0] and [1] so v-loop also ok but =true may pick non-sample w/ v-loop)
//			if (eucSamples && distanceBetween(verts[ samples[i] ]->coords, verts[ samples[j] ]->coords) > maxGeoDist) //this is max euclidean distance so naming sucks
//				maxGeoDist = distanceBetween(verts[ samples[i] ]->coords, verts[ samples[j] ]->coords);
//			else if (! eucSamples && verts[ samples[i] ]->geodesic[ samples[j] ] > maxGeoDist)
			if (verts[ samples[i] ]->geodesic[ samples[j] ] > maxGeoDist)
				maxGeoDist = verts[ samples[i] ]->geodesic[ samples[j] ];
		//store paths from samples[i] to each sample (to each vertex causes memory overflow) to avoid dijkstra recomputations later (prev values must be updated for each new sourceVert so doing this later causes dijkstra recomputations)
		pathsFromVertex(samples[i]); //verts[ samples[i] ]->tmp = (i / 30); //for visualization only
	}
	if (nReplacements > 0)
		cout << "max geo dist on this ROI submesh: " << maxGeoDist << " (nReplace: " << nReplacements << ") in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs\n"; //(eucSamples ? "Euclidean" : "geodesic") << "
	else
		cout << "max geo dist on this ROI submesh: " << maxGeoDist << " in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs\n"; //(eucSamples ? "Euclidean" : "geodesic") << "
/*	cout << "pairwise structure:\n";
	for (size_t i = 0; i < samples.size(); i++)
	{
		for (size_t j = 0; j < samples.size(); j++)
			cout << verts[ samples[i] ]->geodesic[ samples[j] ] << "\t";
		cout << endl;
	}//*/
//	return eucSamples;
}

int Mesh::FPSpatch(size_t n, double r, int v, bool curvatureBased, bool untouchedsOnly)
{
	//returns n-size subset of samples[] representing the r-radius patch centered at v
	
	if (n > samples.size())
		cout << "WARNING: nSamplesPerSubregion is too big for the available number of samples\n\n";
	for (size_t s = 0; s < samples.size(); s++)
		verts[ samples[s] ]->marked = false; //true means added to samplesC2F below (needed by rare nextSample==-1 case below); not updating bool sample 'cos that variable is true for all ROI samples (FPSroi())
	samplesC2F.clear(), samplesC2F.push_back(v), verts[v]->marked = true; //vector< int > fpsLocal = {v};
	vector< pair<double, int> > pool;
	int nOffRadius = 0; //number of samples added by nextSample==-1 case below which does not respect the radius constraint
	while (samplesC2F.size() < n)
	{
		double maxDist = -INF;
		int nextSample = -1;
		pool.clear();
		for (size_t i = 0; i < samples.size(); i++) //fps selects from amongst all samples ('cos all the unary/pairwise descriptors are already precomputed for the samples)
		{
			size_t s = samples[i];
			if (verts[v]->geodesic[s] > r || //normalized geo distance vs. normalized radius (c2f radius not the patch radius)
				(untouchedsOnly && verts[v]->geoToRef != -1.0)) //process untoucheds only but v is touched so skip it (geoToRef!=-1 means touched, naming sucks)
				continue; //cannot be in the current patch 'cos s is too far away from the seed v
			bool exists = false; //exists in samplesC2F
			for (size_t j = 0; j < samplesC2F.size(); j++)
				if (samplesC2F[j] == s)
					exists = true;
			if (exists)
				continue;
			int e2 = -1;
			double minDist = INF;
			for (size_t e = 0; e < samplesC2F.size(); e++) //find the existing sample e \in samplesC2F that is closest to s \in samples
				if (verts[ samplesC2F[e] ]->geodesic[s] <= minDist)
					minDist = verts[ samplesC2F[e] ]->geodesic[s], e2 = samplesC2F[e];
//cout << "existing " << e2 << " is closest to candid " << s << " w/ geo: " << minDist << "\tcandid's geoToCenter: " << verts[v]->geodesic[s] << endl;
			//instead of just taking the farthest point, we first create a pool of the top 1% (see 10 candidates below) farthest points, and from this pool we then select the one with highest curvature
			if (curvatureBased) //instead of making nextSample=s if verts[e2]->geodesic[s] is max (else if part below), insert verts[e2]->geodesic[s] value to the pool and later use the 1% of this pool (to pick the max curvature item as the nextSample)
				pool.push_back({verts[e2]->geodesic[s], s});
			else if (verts[e2]->geodesic[s] > maxDist)
				maxDist = verts[e2]->geodesic[s], nextSample = s;
		} //end of for i
		if (curvatureBased && ! pool.empty()) //ineffective here as pool.size() = nDenseSamples = 30 (typically) which makes 30*0.01=0 below
		{
			//now sort the pool (in descending order) to get the top 1% farthest points
			sort(pool.begin(), pool.end(), comparator); //if comparator not provided then sorts in ascending order
//			if (verts[ pool[0].second ]->roiBorder && ! verts[ pool[0].second ]->border) //not a regular border that may confuse areasRing (borderTreat=1) but an informative roiBorder so select it regardless of its gc
//				nextSample = pool[0].second; not that critical so just make it consistent with FPS() and FPSroi()
//			else
//			{
				double maxGC = -INF; //nextSample will be the max curvature item from the top 1% of the pool
				//for (size_t i = 0; i < pool.size()*0.01; i++)
				for (size_t i = 0; i < 5/*10*/ && i < pool.size(); i++) //pool.size() = nDenseSamples = 30 (typically) which makes 30*0.01=0 so use a fixed number of 10 candidates here
					//if (abs(verts[ pool[i].second ]->gc) > maxGC) //still meaningful as very negative gc means deep/radical saddle points which are also important
					if (verts[ pool[i].second ]->gc > maxGC) //locMaxGC=true made based on non-abs values so be consistent
						//maxGC = abs(verts[ pool[i].second ]->gc), nextSample = pool[i].second; use this if abs in if is enabled
						maxGC = verts[ pool[i].second ]->gc, nextSample = pool[i].second;
				if (nextSample == -1 || verts[nextSample]->border)//&& ! verts[ pool[0].second ]->border) //border.gc not reliable so replace it w/ the farthest point, i.e., no curvatureBased decision, at index 0 (descending order)
					nextSample = pool[0].second; //pool[0] allowed to be border (2nd condition disabled) but if it was not allowed then i'd have checked its 1-ring for a non-border
//			}
		}
		if (nextSample == -1)
		{
			//fill the rest by relaxing r-distance constraint, i.e., insert the next closest(s) even if geodesic from v to that is longer than r
			while (samplesC2F.size() < n)
			{
				double minDist = INF;
				for (size_t s = 0; s < samples.size(); s++)
					if (! verts[ samples[s] ]->marked && verts[v]->geodesic[ samples[s] ] < minDist)
					{
						if (untouchedsOnly && verts[ samples[s] ]->geoToRef != -1.0) //!=-1 means samples[s] is touched so cannot be the nextSample in untouchedsOnly mode
							continue;
						minDist = verts[v]->geodesic[ samples[s] ], nextSample = samples[s];
					}
				if (nextSample == -1) //no valid nextSample due to the untouchedsOnly condition so fill the rest below
					break; //samplesC2F.size() < n for sure
				verts[nextSample]->marked = true, samplesC2F.push_back(nextSample), nOffRadius++;
//cout << v << ": r-relaxed and added " << nextSample << "\t" << minDist << "\t" << verts[nextSample]->roiBorder << " & " << verts[nextSample]->border << endl;
			}
			if (nextSample == -1) //always false in !untouchedsOnly mode, i.e., no extra loop below is needed 'cos samplesC2F.size() is definitely n when nextSample!=-1
			{
				while (samplesC2F.size() < n) //fill the rest after removing the geoToRef/untouched condition
				{
					double minDist = INF;
					for (size_t s = 0; s < samples.size(); s++)
						if (! verts[ samples[s] ]->marked && verts[v]->geodesic[ samples[s] ] < minDist)
							minDist = verts[v]->geodesic[ samples[s] ], nextSample = samples[s];
					if (nextSample == -1)
						cout << "WARNING: cannot find a valid nextSample (due to n > samples.size)\n\n";
					verts[nextSample]->marked = true, samplesC2F.push_back(nextSample), nOffRadius++;
//cout << v << ": r-relaxed and added " << nextSample << "\t" << minDist << "\t" << verts[nextSample]->roiBorder << " & " << verts[nextSample]->border << endl;
				}				
			}
			break;
		}
		verts[nextSample]->marked = true, samplesC2F.push_back(nextSample); //dijkstraShortestPaths(nextSample); //upcoming iteration needs fpsLocal.geodesic which is already computed by fpsROI() so disable here
//cout << "sample # " << samplesC2F.size() << ": " << nextSample << "\t" << verts[v]->geodesic[nextSample] << "\t" << verts[nextSample]->roiBorder << " & " << verts[nextSample]->border << endl; //first sample added above before the loop
	}
	return nOffRadius;
}

void Mesh::FPSfile(size_t N, bool competitorResult)
{
	//loads N FPS samples from file fName s that i can compare with qmatch results using the same sample set

	char fName[250];
	sprintf(fName, "qmatchStyle-results\\%d-n%d solqmatch.txt", id, N);
	sprintf(fName, "qmatchStyle-results\\%d-n%d solqucoop.txt", id, N);
	sprintf(fName, "qmatchStyle-results\\%d-n%d solqucm2.txt", id, N);
	if (! competitorResult) //our result so our proram will be executed on these file-based samples (same as competitor's samples for fairness)
		sprintf(fName, "qmatchStyle-results\\%d-n%d solqucm.txt", id, N);
	FILE* fPtr;
	if (! (fPtr = fopen(fName, "r")))
		cout << "cannot read " << fName << endl, exit(0);
	cout << "File-based evenly-spaced sampling [" << N << "].. ";
	samples.clear();
	for (size_t v = 0; v < verts.size(); v++)
		verts[v]->geodesic.clear(), verts[v]->sample = false;
	int a;
	while (fscanf(fPtr, "%d\n", &a) != EOF) //go until the end of file
		verts[a]->sample = true, samples.push_back(a), dijkstraShortestPaths(a);
		
	//set maxGeoDist properly which in turn gives me a chance to normalize geodesics if necessary
	maxGeoDist = -INF; //maxGeoDist can be from sample to non-sample vertex, so after normalization geodesic b/w most-separated sample pair may remain < 1 (not the case anymore; see j-loop below)
	for (size_t i = 0; i < samples.size(); i++)
	{
		for (size_t j = 0; j < samples.size(); j++)//for (size_t v = 0; v < verts.size(); v++) matching samples to samples so let maxGeoDist b/w sample to sample (curvatureBased=false selects samples[0] and [1] so v-loop also ok but =true may pick non-sample w/ v-loop)
			if (verts[ samples[i] ]->geodesic[ samples[j] ] > maxGeoDist)
				maxGeoDist = verts[ samples[i] ]->geodesic[ samples[j] ];//, max1=samples[i], max2=v;
		pathsFromVertex(samples[i]);
	}
	cout << "max geo dist: " << maxGeoDist << endl;//" from " << max1 << " to " << max2 << endl;	
	avgGeoDists(); //to see the agd values of each of the loaded samples	
	for (size_t i = 0; i < samples.size(); i++)
		cout << samples[i] << "\t" << verts[ samples[i] ]->agd << "\t" << verts[ samples[i] ]->coords[0] << " " << verts[ samples[i] ]->coords[1] << " " << verts[ samples[i] ]->coords[2] << endl;
}

void Mesh::FPSroi(size_t k, const vector< int >& prevSamples)
{
	//same as FPSroi() except prevSamples are not added as new samples here although they affect the locations of the new k samples

	cout << "Curv. evenly-spaced sampling [" << k << "] on ROI..";
	clock_t t = clock();
	size_t offset = (mesh2 ? 1 : 0);
	vector< pair<double, int> > pool;
	k = (k > nROIverts ? nROIverts : k); //k cannot exceed number of ROI vertices 'cos only they are eligible to be samples in this function
	if (prevSamples.empty()) //very first call to this function so do clearSamples() and locMaxSampling()
	{
		clearSamples(offset); //marked values set too
		//fill the rest of the samples (clearSamples() keeps the coarseMap[], not coarseMapFull[], samples), first w/ locMaxSampling() and then typical FPS below
		//first add locMaxGC vertices (low-curvature ones and close ones trimmed by locMaxGCs()) to samples
		locMaxSampling(k); //marked values set too
	}
	else
	{
		//fill all the samples w/ typical FPS below which will see that previous samples are already inside samples[] so adds the new/current samples farthest away from them
		for (size_t v = 0; v < verts.size(); v++)
			//verts[v]->sample = verts[v]->marked = false; //no v.geodesic.clear() to protect previous samples' geodesics; since protected geodesics are normalized they're not compatible w/ unnormalized nextSample.geodesic (after dijkstra) below so do v.geodesic.clear()
			verts[v]->geodesic.clear(), verts[v]->sample = verts[v]->marked = false; //do dijkstraShortestPaths(prevSamples[i]) below since geodesics are cleared now (inefficient recomputaions but accurate/safe)
		samples.clear();
		for (size_t i = 0; i < prevSamples.size(); i++)
		{
			verts[ prevSamples[i] ]->sample = true, samples.push_back(prevSamples[i]), dijkstraShortestPaths(prevSamples[i]);
			if (verts[ prevSamples[i] ]->geodesic.empty())
				cout << "WARNING: empty geodesics not expected for previous samples\n\n";
			for (size_t v = 0; v < verts.size(); v++)
				if (verts[ prevSamples[i] ]->geodesic[v] < radius) //radius*1.5 in replaceSample() is different as it's allowed to search for a good replacement in a larger vicinity; this 1*radius must match w/ the 1*radius in FPS()
					verts[v]->marked = true; //more marked=true will be made during FPSroi() below as it proceeds
		}
	}
	//now fill the remaining items via fps for uniform distribution that lets nice coverage on models
	int nReplacements = 0, nAdds = 0;
	do
	{
		//exit condition based on the very first call (prevSamples.empty() and other calls)
		if (prevSamples.empty() && samples.size() >= k) //==k also good but >=k handles push_back(coarseMap[i + offset])s in clearSamples() (which should never exceed k=30 but anyway)
			break;
		if (! prevSamples.empty() && nAdds++ == k) //not the very first call
			break;
		double maxDist = -INF;
		int nextSample = -1;
		pool.clear();
		for (size_t s = 0; s < verts.size(); s++) //fps selects from amongst all vertices
		{
			if (! verts[s]->roi || //process only the ROI vertices during this sampling
				verts[s]->edgeList.empty() || //isolated verts should never be samples (even if they'd have dSapnning=INF values); they occur on microholes/holes/viewing classes of shrec11 dataset
				verts[s]->sample) //pure FPS not need this but coming here from another sampler may insert duplicates that are already made sample by that sampler
				continue;
			if (fpsFromGT && verts[s]->gtIdx == -1) //for fair evaluation select amongs gt-equipped vertices
				continue;
			int e2 = -1;
			double minDist = INF;
			for (size_t e = 0; e < samples.size(); e++) //find the existing sample e \in samples that is closest to s \in verts
				if (verts[ samples[e] ]->geodesic[s] <= minDist)
					minDist = verts[ samples[e] ]->geodesic[s], e2 = samples[e];			
			pool.push_back({verts[e2]->geodesic[s], s}); //instead of making nextSample=s if verts[e2]->geodesic[s] is max, insert verts[e2]->geodesic[s] value to the pool and later use the 1% of this pool (to pick the max curvature item as the nextSample)
		} //end of for s
		if (! pool.empty())
		{
			//now sort the pool (in descending order) to get the top 1% farthest points
			sort(pool.begin(), pool.end(), comparator); //if comparator not provided then sorts in ascending order
			double maxGC = -INF; //nextSample will be the max curvature item from the top 1% of the pool
			for (size_t i = 0; i < pool.size()*0.01; i++) //too large pool hurts evenly spacing, e.g., puts samples to all fingers makign small hand region excessively dense sampled
				//if (abs(verts[ pool[i].second ]->gc) > maxGC) //still meaningful as very negative gc means deep/radical saddle points which are also important
				if (verts[ pool[i].second ]->gc > maxGC) //locMaxGC=true made based on non-abs values so be consistent
					//maxGC = abs(verts[ pool[i].second ]->gc), nextSample = pool[i].second; use this if abs in if is enabled
					maxGC = verts[ pool[i].second ]->gc, nextSample = pool[i].second;
			if (nextSample == -1 || verts[nextSample]->border)//&& ! verts[ pool[0].second ]->border) //border.gc not reliable so replace it w/ the farthest point, i.e., no curvatureBased decision, at index 0 (descending order)
				nextSample = pool[0].second; //pool[0] allowed to be border (2nd condition disabled) but if it was not allowed then i'd have checked its 1-ring for a non-border
		}
		if (nextSample == -1)
		{
			cout << "WARNING: FPS could not find a valid sample; quiting early (pool[] may be empty)\n\n";
			break;
		}
		verts[nextSample]->sample = true;
		samples.push_back(nextSample);
		dijkstraShortestPaths(nextSample); //upcoming iteration needs samples.geodesic
		if (prevSamples.empty()) //almost all k=30 are replaced for late iterations as those samples are not as carefully selected, i.e., no locMaxSampling(), as the very first iteration so skip replaceSample() to save time (regardless of skipping n=250 samples reached by denseMap() and denseMapLargeN() will be different)
			nReplacements += replaceSample(nextSample); //replace nextSample with an appropriate max-curvature locMaxGC vertex
		for (size_t v = 0; v < verts.size(); v++)
			if (verts[ samples[ samples.size() - 1 ] ]->geodesic[v] < radius) //last sample (either nextSample or its replacement that is local to replaceSample())
				verts[v]->marked = true;
	} while(true);
	//erase the prevSamples which occupies the first prevSamples.size() entries in samples (they were necessary to put the new/current samples farthest away from them)
    if (prevSamples.size() <= samples.size()) //safety check to ensure prevSamples.size() does not exceed vector size
        samples.erase(samples.begin(), samples.begin() + prevSamples.size());
	else
		cout << "WARNING: too many prevSamples\n\n";
	cout << "path verts..";
	if (prevSamples.empty()) //very first call to this function so set maxGeoDist
	{
		//set maxGeoDist properly which in turn gives me a chance to normalize geodesics if necessary; stores path vertices b/w pairs of samples too
		maxGeoDist = -INF; //maxGeoDist can be from sample to non-sample vertex, so after normalization geodesic b/w most-separated sample pair may remain < 1 (not the case anymore; see j-loop below)	
		for (size_t i = 0; i < samples.size(); i++)
		{
			for (size_t j = 0; j < samples.size(); j++)//for (size_t v = 0; v < verts.size(); v++) matching samples to samples so let maxGeoDist b/w sample to sample (curvatureBased=false selects samples[0] and [1] so v-loop also ok but =true may pick non-sample w/ v-loop)
				if (verts[ samples[i] ]->geodesic[ samples[j] ] > maxGeoDist)
					maxGeoDist = verts[ samples[i] ]->geodesic[ samples[j] ];
			//store paths from samples[i] to each sample (to each vertex causes memory overflow) to avoid dijkstra recomputations later (prev values must be updated for each new sourceVert so doing this later causes dijkstra recomputations)
			pathsFromVertex(samples[i]); //normalization of the samples[i].dijkstra's made inside pathsFromVertex()'ll be done by denseMapLargeN() on return
		}
	}
	else //just paths no maxGeoDist for this intermediate function call
	{
		for (size_t i = 0; i < samples.size(); i++)
			//store paths from samples[i] to each sample (to each vertex causes memory overflow) to avoid dijkstra recomputations later (prev values must be updated for each new sourceVert so doing this later causes dijkstra recomputations)
			pathsFromVertex(samples[i]); //verts[ samples[i] ]->tmp = (i / 30); //for visualization only
		//coarseMap[] was definitely included in samples[] (due to clearSamples()) for the prevSamples.empty() case above (very first call); now samples[] changed so coarseMap[] samples must recompute their paths to those new samples as done below
		for (size_t i = 0; i < coarseMap.size(); i += 2) //each 2 consecutive items make 1 match (even indices belong to mesh1, odd indices mesh2)
			pathsFromVertex(coarseMap[i + offset]); //normalization of the coarseMap[i].dijkstra's made inside pathsFromVertex()'ll be done by denseMapLargeN() on return
	}
	if (nReplacements > 0)
		cout << "max geo dist on this ROI submesh: " << maxGeoDist << " (nReplace: " << nReplacements << ") in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs\n"; //(eucSamples ? "Euclidean" : "geodesic") << "
	else
		cout << "max geo dist on this ROI submesh: " << maxGeoDist << " in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs\n"; //(eucSamples ? "Euclidean" : "geodesic") << "
}

void Mesh::pathsFromVertex(int sourceVert)
{
	//stores path vertices from sourceVert to each target sample t in sourceVert.pathVerts[idx of t in samples] (stored such that local pathVerts[0] is t and pathVerts[last] is sourceVert)

//if (! verts[sourceVert]->pathVerts.empty())return; cannot be done 'cos once samples[] change, e.g., from 10x10 initial coarse mode to 250x250 dense mode, pathVerts[] needs to be refilled; resize() deallocs the old memory and creates a new memory of new capacity items
	verts[sourceVert]->pathVerts.clear(); //all vertices on the path going from sourceVert to the ith sample will be stored in sourceVert.pathVerts[i] as one vector
	//get the path that goes from sourceVert to target t using the prev's set during the following dijkstra (includes start (sourceVert) and end (t) vertices)
	verts[sourceVert]->geodesic.clear(), dijkstraShortestPaths(sourceVert); //the 1st resize() by FPS/FPSroi() already allocated memory of size V and this clear sets size=0 but not capacity=0 so this 2nd call reuses the existing capacity, i.e., zero memory allocations for the 2nd time
	for (size_t i = 0; i < samples.size(); i++) //for (size_t t = 0; t < verts.size(); t++) is replaced by sample-based traversal to prevent memory overflow; all i need is sample-to-sample paths anyway so store N=250 entries instead of V=52565
	{
		vector< int > pathVerts; //path vertices local to this loop
		int t = samples[i], a = t; //from a to a.prev, and then a.prev becomes the new a
		pathVerts.push_back(a); //includes end vertex (a=t)
		while (a != sourceVert)
		{
			pathVerts.push_back(a); //path vertex stored for nonuniform sampling purposes only
			if (verts[a]->prev == -1)
				cout << "WARNING: prev of " << a << " is not defined for " << sourceVert << " to " << t << "\n\n";
			a = verts[a]->prev;				
		}
		if (a != t) //a = t happens when path is from t to itself in which case it must only contain t (duplication occurs w/o this condition 'cos t already added above)
			pathVerts.push_back(a); //a = sourceVert for sure 
		verts[sourceVert]->pathVerts.push_back(pathVerts); //vertices on the path going from sourceVert to t=samples[i] are stored in sourceVert.pathVerts[i] now
	}
}

void Mesh::clearSamples(size_t offset)
{
	//clears all samples except the coarseMap[] samples

	for (size_t v = 0; v < verts.size(); v++)
		verts[v]->geodesic.clear(), verts[v]->sample = verts[v]->marked = false; //marked=true means v is in the vicinity of an existing FPSroi sample (for replaceSample() only; FPSroi() is already evenly-spaced (curvatureBased slightly distracts it but anyway) so this marked is for replaceSample() that has no evenly-spacing respect)
	samples.clear();
	for (size_t i = 0; i < coarseMap.size(); i += 2) //each 2 consecutive items make 1 match (even indices belong to mesh1, odd indices mesh2)
	{
		verts[ coarseMap[i + offset] ]->sample = true, samples.push_back(coarseMap[i + offset]), dijkstraShortestPaths(coarseMap[i + offset]);
		if (! verts[ coarseMap[i + offset] ]->roi)
			cout << "WARNING: matched sample must have been included in the ROI\n\n";
		for (size_t v = 0; v < verts.size(); v++)
			if (verts[ coarseMap[i + offset] ]->geodesic[v] < radius) //radius*1.5 in replaceSample() is different as it's allowed to search for a good replacement in a larger vicinity; this 1*radius must match w/ the 1*radius in FPS()
				verts[v]->marked = true; //more marked=true will be made during FPSroi() as it proceeds
	}
}

bool Mesh::replaceSample(int s)
{
	//replaces sample verts[s] with the locMaxGC vertex that is closest to it (heuristic=1) or with the max-curvature locMaxGC vertex that is within radius distance to it (heuristic=2) or no replacement (0)

	int heuristic = 2/*0 1 2*/, newReplacement1 = -1, newReplacement2 = -1, replaced = 0;
	if (heuristic == 0)
		return false;
	double minDist = INF, maxGC = -INF;
	for (size_t v = 0; v < verts.size(); v++)
	{
		if (! verts[v]->roi || //process only the ROI vertices during this sampling
			verts[v]->edgeList.empty() || //isolated verts should never be samples (even if they'd have dSapnning=INF values); they occur on microholes/holes/viewing classes of shrec11 dataset
			verts[v]->sample || verts[v]->marked || //only nonsamples and nonmarked can be the replacements ('cos marked=true means v is in the vicinity of an existing FPS sample so can't be a new replacement sample)
			//! verts[v]->locMaxGC) //only local maxima of gaussian curvature can be the replacements
			! verts[v]->locMaxGCrelaxed) //use relaxed local maxima to increase chances of selecting a good replacement within radius (also marked condition decreases chances even for non-radius heuristic1 so increse those chances too by using a larger/relaxed local maxima set)
			continue;
		if (fpsFromGT && verts[v]->gtIdx == -1) //for fair evaluation select amongs gt-equipped vertices
			continue;
		if (verts[s]->geodesic[v] < minDist) //s to v is 0 distance when v=s which is ok 'cos that v passes the locMaxGC test above so deserves to be a sample
			minDist = verts[s]->geodesic[v], newReplacement1 = v; //cannot pick both closest and max-curvature because the closest one w/ tiny-curvature will prevent far-away max-curvature
		if (verts[s]->geodesic[v] < radius * 1.5 && verts[v]->gc > maxGC)
			maxGC = verts[v]->gc, newReplacement2 = v;
	}
	if (heuristic == 1 && newReplacement1 != -1 && newReplacement1 != s) //-1 means no replacement happened
		verts[s]->sample = false, verts[s]->geodesic.clear(), samples[ samples.size() - 1 ] = newReplacement1, verts[newReplacement1]->sample = true, dijkstraShortestPaths(newReplacement1), replaced = 1;
	if (heuristic == 2 && newReplacement2 != -1 && newReplacement2 != s) //-1 means no replacement happened
		verts[s]->sample = false, verts[s]->geodesic.clear(), samples[ samples.size() - 1 ] = newReplacement2, verts[newReplacement2]->sample = true, dijkstraShortestPaths(newReplacement2), replaced = 1;
//cout << (replaced == 1 ? newReplacement2 : s) << " is the next sample\t" << replaced << endl;
//if (newReplacement2 != -1)cout << s << " --> " << newReplacement2 << " || " << verts[newReplacement2]->geodesic[s] << " < " << radius*1.5 << "\n"; //2144 --> 2189 || 14.2789 < 15.0135
	return replaced == 1;
}

void Mesh::locMaxSampling(int N, bool denseSampling)
{
	//makes local maxima of gc values samples (low-curvature ones and close ones are already trimmed by locMaxGCs() but spaced=true prevents still-close, e.g., 2 fingertips, samples from being selected below)

	bool spaced = true; //prevent still-close samples
	int nLocMaxSamples = 0;
//	for(size_t v=0;v<verts.size();v++)verts[v]->marked=verts[v]->sample; //true means ineligible to be a new sample below (clearSamples() (or extremeVert() if !denseSampling) already marked samples from coarseMap[] and their vicinity so don't update that)
	while (true)
	{
		double maxGC = -INF;
		int winner = -1;
		for (size_t v = 0; v < verts.size(); v++)
			if (! verts[v]->marked && verts[v]->locMaxGC && verts[v]->gc > maxGC && verts[v]->roi) //last condition ensures to process only the ROI vertices during this sampling
				maxGC = verts[v]->gc, winner = v;
		if (winner == -1 || samples.size() == N) //2nd condition to prevent exceeding the user-specified number of samples
			break;
		else
			samples.push_back(winner), verts[winner]->sample = verts[winner]->marked = true, dijkstraShortestPaths(winner), nLocMaxSamples++;
		if (spaced) //mark radius-vicinity of the newcomer winner as ineligible, i.e., another localMaxGC in that vicinity cannot be made sample by this function (it may be a sample later by FPS/FPSroi() depending on the density of sampling)
			for (size_t v = 0; v < verts.size(); v++)
				if (verts[winner]->geodesic[v] < radius * (denseSampling ? 1.5/*1.0 1.5*/ : 3)) //make a lot of locMaxGC ineligible (in !denseSampling mode, 1.5radius keeps high-curvature fingertips alive and fails to add relatively low-curvature nipples due to N size limit constraint so use 3 there); see the 1*radius comment in clearSamples() too
					verts[v]->marked = true; //ineligible 'cos in vicinity if(verts[v]->locMaxGC)cout<<winner<<" made locMaxGC "<<v<<" ineligible\n";
	}
	cout << nLocMaxSamples << " locMax added..";
}

void Mesh::pathSampling(size_t n)
{
	//samples each path running b/w samples uniformly (or nonuniformly based on gc) via n+1 vertices (there'll be n+1 (including path endpoints) samples on each path (small/short paths may remain under n+1))

	cout << "path sampling..";
	clock_t t = clock();
	size_t nPathsInUse = 0; //# valid paths (too short ones are invalid, i.e., the ones that are unable to provide n+1 uniform vertices)
	for (size_t v = 0; v < verts.size(); v++)
		verts[v]->area = verts[v]->patchGC = -1.0, verts[v]->patch.clear(), verts[v]->pathSamples.clear(); //will be made different than -1 below for pathSamples only
	for (size_t i = 0; i < samples.size(); i++)
	{
		//path info from samples[i] to all other vertices v is already stored in samples[i].pathVerts[v] so no dijkstra recomputation needed here
		for (size_t j = 0; j < samples.size(); j++)
		{
			vector< int > pathSamples; //path samples local to this loop (samples[i].pathVerts[v] is needed by nonuniform too)
			if (j < i) //symmetric part must be same as the previously computed i<j part, e.g., [3].pathSamples[0] will get its content from [0].pathSamples[3] where i=3, j=0
			{
				for (int k = (int) verts[ samples[j] ]->pathSamples[i].size() - 1; k >= 0; k--) //cannot use size_t due to -1 potential
					pathSamples.push_back(verts[ samples[j] ]->pathSamples[i][k]); //reverse(pathSamples.begin(), pathSamples.end()); could have been used but starting loop from size()-1 already makes the desired reverse effect
				verts[ samples[i] ]->pathSamples.push_back(pathSamples); //verts[ samples[i] ]->pathSamples[j] = verts[ samples[j] ]->pathSamples[i]; effect
				nPathsInUse += (pathSamples.size() == n+1);
				continue;
			}
			double factor = verts[ samples[i] ]->geodesic[ samples[j] ] / n;
			vector< int > pathVerts = verts[ samples[i] ]->pathVerts[j]; //all vertices on the path going from the ith sample to the jth sample are already stored in samples[i].pathVerts[j]
			int k = 0, a = pathVerts[k++]; //from a to a.prev, and then a.prev becomes the new a; recall path from sourceVert to t is stored such that [0] is t and [last] is sourceVert
			pathSamples.push_back(a); //includes end vertex (a=samples[j])
			while (a != samples[i])
			{
				a = pathVerts[k++]; //advance to the next path vertex
				//nth sample selection (whenever i exceed the currentSize*factor for the first time, i hit the next desired uniformly distributed nth sample on the path)
				if (verts[ samples[j] ]->geodesic[a] > pathSamples.size()*factor && a != samples[i]) //distance to samples[j]; 2nd condition: samples[i]'ll be added after the loop so not here
					pathSamples.push_back(a);
			}
			if (a != samples[j]) //a = samples[j] happens when i=j, i.e., path from sample to itself must only contain that sample itself (duplication occurs w/o this condition 'cos samples[j] already added above)
				pathSamples.push_back(a); //a = samples[i] for sure 
			verts[ samples[i] ]->pathSamples.push_back(pathSamples); //samples on the path going from the ith sample to the jth sample are given by samples[i].pathSamples[j]
			//for printing only
			if (pathSamples.size() > n+1)
				cout << "WARNING: excessive # of uniform path verts: " << pathSamples.size() << " > " << n << "\n\n";
			else if (pathSamples.size() == n+1) //long-enough path has provided desired amount of uniform verts
				nPathsInUse++; //else nTooShort++;		
			//define patches around each path samples
			for (size_t z = 0; z < pathSamples.size(); z++)
				patchAndDescriptors(pathSamples[z]);	
		}
	}
	//process path from samples[i] to samples[j] again to add new samples from samples[j] to coarseMap.samples[0, m-1] during denseMap().pathSampling() where coarseMap is available/non-empty (samplesToCoarse is empty if coarseMap is empty)
	size_t offset = (mesh2 ? 1 : 0);
	if (! coarseMap.empty() && extraSamples) //technically pathSamples added here are not on the path from samples[i] to samples[j] but very useful to represent the pairwise structure b/w samples[i] and samples[j]
		for (size_t i = 0; i < samples.size(); i++)
			for (size_t j = 0; j < samples.size(); j++) //path that goes from any other sample samples[j] to samples[i]
				for (size_t m = 0; m < verts[ samples[j] ]->samplesToCoarse.size(); m++) //path samples from this sample (currently samples[j]) to each coarseMap sample (currently the m'th sample since verts[ samples[j] ]->samplesToCoarse.size = coarseMap.size/2)					
					if (samples[i] != coarseMap[2*m + offset]) //verts[ samples[j] ]->samplesToCoarse[m=0] below is the set of samples from samples[j] to coarseMap[0+offset]; so skip this set if samples[i]=coarseMap[0+offset] 'cos samples[i]-samples[j] path is already sampled inside verts[ samples[i] ]->pathSamples[j]
						for (size_t k = 1; k < verts[ samples[j] ]->samplesToCoarse[m].size(); k++) //start from 1 to skip [0] = samples[j] which is already inside verts[ samples[i] ]->pathSamples[j] (added above)
							verts[ samples[i] ]->pathSamples[j].push_back(verts[ samples[j] ]->samplesToCoarse[m][k]), patchAndDescriptors(verts[ samples[j] ]->samplesToCoarse[m][k]); //define patch around this new path sample
	if (nonuniformPathSamplingHeuristic == 0)
	{
		cout << "uniform; # valid paths: " << nPathsInUse << " / " << samples.size() * samples.size() << " in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs due to " << (eucPatching ? "euclidean" : "geodesic") << " patching\n";
		return;
	}
//for (size_t i = 0; i < coarseMap.size(); i += 2) cout << coarseMap[i] << " -coarseMap- " << coarseMap[i+1] << endl;int v1 = samples[0], idx = 1;
//if(id==0&&samples.size()>20){v1 = samples[20], idx = 21; cout<<coarseMap.size()<<"|"<<v1<<"to"<<samples[idx]<<": pathSamples before: "; for (size_t k = 0; k < verts[v1]->pathSamples[idx].size(); k++) cout << verts[v1]->pathSamples[idx][k] << " "; cout << "\n\n";}
	//for each patch p (except the path endpoint patches since they represent actual/interesting samples), replace its centering sample w/ the max-curvature locMaxGC path vertex inside p which leads to nonuniform importance sampling (nonuniformPathSamplingHeuristic=1)
	//additional heuristic: replace the centering sample with the max-curvature locMaxGC patch (not path) vertex (nonuniformPathSamplingHeuristic=2) which makes us go off the path but not that off so sill ok (illogical to pick both closest and max-curvature so just use max-curvature)
	if (nonuniformPathSamplingHeuristic == 1)
		for (size_t i = 0; i < samples.size(); i++)
			for (size_t j = 0; j < samples.size(); j++)
			{
				for (size_t v = 0; v < verts.size(); v++)
					verts[v]->pathVert = false; //true for the vertices on the path from i to j (samples[i].pathVerts[j]); not using marked as patchAndDescriptors() updates marked values
				for (size_t v = 0; v < verts[ samples[i] ]->pathVerts[j].size(); v++)
					verts[ verts[ samples[i] ]->pathVerts[j][v] ]->pathVert = true;
				size_t k = 0; //endpoints (k=0 (always) and n-1 (if !extraSamples)) will not be processed by newSample loop below so make their patches unavailable for new sample selection via pathVerts=false
				//k = verts[ samples[i] ]->pathSamples[j].size() - 1; other endpoint (correct only if !extraSample so use the search loop below to be correct at all times)
				for (; k < verts[ samples[i] ]->pathSamples[j].size(); k++) //k=0 is an endpoint but after extraSamples the other endpoint is not in the end so use this loop to find it
					if (verts[ samples[i] ]->pathSamples[j][k] == samples[i]) //samples[i] end
						break;
				for (size_t m = 0; m < verts[ verts[ samples[i] ]->pathSamples[j][k] ]->patch.size(); m++)
					verts[ verts[ verts[ samples[i] ]->pathSamples[j][k] ]->patch[m] ]->pathVert = false; //now it cannot be selected by an upcoming maxGC loop iteration due to the first condition
				for (k = 0; k < verts[ samples[i] ]->pathSamples[j].size(); k++) //k=0 is an endpoint but after extraSamples the other endpoint is not in the end so use this loop to find it
					if (verts[ samples[i] ]->pathSamples[j][k] == samples[j]) //samples[j] end
						break;	
				for (size_t m = 0; m < verts[ verts[ samples[i] ]->pathSamples[j][k] ]->patch.size(); m++)
					verts[ verts[ verts[ samples[i] ]->pathSamples[j][k] ]->patch[m] ]->pathVert = false; //now it cannot be selected by an upcoming maxGC loop due to the first condition
//				for (k = 1; k < verts[ samples[i] ]->pathSamples[j].size() - 1; k++) //first (k=0) and last (size-1) skipped 'cos they are endpoints and hence actual samples (no replacement for them as they are already interesting fps samples)
				for (k = 0; k < verts[ samples[i] ]->pathSamples[j].size(); k++) //all path samples are replacable (including 2 path endpoints); irreplacable now thanks to the continue below
				{
					int newSample = -1, oldSample = verts[ samples[i] ]->pathSamples[j][k]; //may be replaced w/ the max-curvature locMaxGC path vertex inside patch (if no locMaxGC then just the max-curvature vertex)
					if (verts[oldSample]->sample)
						continue; //k=0 and k=last is definitely sample; an intermediate k might also be made sample during FPS/FPSroi() so don't update any existing samples (no replacement)
					size_t patchSize = verts[oldSample]->patch.size();
					if (patchSize == 0)
						cout << "WARNING: empty patch for " << oldSample << "\n\n";
					double maxGC = -INF;
					for (size_t m = 0; m < patchSize; m++)
//if (verts[ verts[oldSample]->patch[m] ]->pathVert && verts[ verts[oldSample]->patch[m] ]->locMaxGCrelaxed && abs(verts[ verts[oldSample]->patch[m] ]->gc) > maxGC)//abs(verts[ verts[oldSample]->patch[m] ]->gc) > maxGC //abs(gc) 'cos very negative curvatures also indicate interesting nonplanar regions, e.g., deep/radical saddle points
						if (verts[ verts[oldSample]->patch[m] ]->pathVert && verts[ verts[oldSample]->patch[m] ]->locMaxGCrelaxed && verts[ verts[oldSample]->patch[m] ]->gc > maxGC  //abs may block an informative high positive curvature vertex so no abs
							&& ! duplicateAhead(k, i, j, verts[oldSample]->patch[m])) //34662 34399 34706 35360 2309 --> 34399 could be replaced by 34706 ahead leading to duplicates (prevent that by this condition)
						{
//maxGC = abs(verts[ verts[oldSample]->patch[m] ]->gc), newSample = verts[oldSample]->patch[m];
								maxGC = verts[ verts[oldSample]->patch[m] ]->gc, newSample = verts[oldSample]->patch[m];
						}
					if (newSample == -1) //max-curvature locMaxGC vertex cannot be found so just take the max-curvature vertex as the winner
					{
						maxGC = -INF;
						for (size_t m = 0; m < patchSize; m++)
							if (verts[ verts[oldSample]->patch[m] ]->pathVert && verts[ verts[oldSample]->patch[m] ]->gc > maxGC //abs may block an informative high positive curvature vertex so no abs
								&& ! duplicateAhead(k, i, j, verts[oldSample]->patch[m])) //34662 34399 34706 35360 2309 --> 34399 could be replaced by 34706 ahead leading to duplicates (prevent that by this condition)
								maxGC = verts[ verts[oldSample]->patch[m] ]->gc, newSample = verts[oldSample]->patch[m];
					}
					newSample = (newSample != -1 ? newSample : oldSample); //necessary if pathVert=false loop below is enabled; otherwise redundant 'cos newSample cannot be -1 then
					verts[ samples[i] ]->pathSamples[j][k] = newSample; //if newSample = oldSample this assignment is redundant and patchAndDescriptors() returns immediately
					patchAndDescriptors(newSample); //if (i==0&&j==1)cout << oldSample << " --> " << newSample << endl;
					//to make nonuniform sampling more uniform, i.e., prevent too close newSamples, i mark the patch of newSample as unavailable for the next sample selection; to do so pathVert=false is sufficient
					for (size_t m = 0; m < verts[newSample]->patch.size(); m++) //too close samples may still validly arise, e.g., s+1 has newSample=-1 (all patch unavailable) and rolls back to oldSample which may be too close to previous sample s
						verts[ verts[newSample]->patch[m] ]->pathVert = false; //now it cannot be selected by an upcoming maxGC loop due to the first condition
				}
			}
	else //we can go off the path with this nonuniformPathSamplingHeuristic=2
		for (size_t i = 0; i < samples.size(); i++)
			for (size_t j = 0; j < samples.size(); j++)
			{
				for (size_t v = 0; v < verts.size(); v++)
					verts[v]->pathVert = false; //true for the vertices on the path from i to j (samples[i].pathVerts[j]) and also for the patches defined on this path
				for (size_t v = 0; v < verts[ samples[i] ]->pathVerts[j].size(); v++)
					for (size_t w = 0; w < verts[ verts[ samples[i] ]->pathVerts[j][v] ]->patch.size(); w++) //not just the pathVerts but also their patches are marked as available; pathVert=true means available (naming sucks 'cos pathc verts are not on the path, they're in the vicinity though)
						verts[ verts[ verts[ samples[i] ]->pathVerts[j][v] ]->patch[w] ]->pathVert = true;
				size_t k = 0; //endpoints (k=0 (always) and n-1 (if !extraSamples)) will not be processed by newSample loop below so make their patches unavailable for new sample selection via pathVerts=false
				for (; k < verts[ samples[i] ]->pathSamples[j].size(); k++) //k=0 is an endpoint but after extraSamples the other endpoint is not in the end so use this loop to find it
					if (verts[ samples[i] ]->pathSamples[j][k] == samples[i]) //samples[i] end
						break;
				for (size_t m = 0; m < verts[ verts[ samples[i] ]->pathSamples[j][k] ]->patch.size(); m++)
					verts[ verts[ verts[ samples[i] ]->pathSamples[j][k] ]->patch[m] ]->pathVert = false; //now it cannot be selected by an upcoming maxGC loop iteration due to the first condition
				for (k = 0; k < verts[ samples[i] ]->pathSamples[j].size(); k++) //k=0 is an endpoint but after extraSamples the other endpoint is not in the end so use this loop to find it
					if (verts[ samples[i] ]->pathSamples[j][k] == samples[j]) //samples[j] end
						break;
				for (size_t m = 0; m < verts[ verts[ samples[i] ]->pathSamples[j][k] ]->patch.size(); m++)
					verts[ verts[ verts[ samples[i] ]->pathSamples[j][k] ]->patch[m] ]->pathVert = false; //now it cannot be selected by an upcoming maxGC loop iteration due to the first condition
				for (k = 0; k < verts[ samples[i] ]->pathSamples[j].size(); k++) //all path samples are replacable (including 2 path endpoints); irreplacable now thanks to the continue below
				{
					int newSample = -1, oldSample = verts[ samples[i] ]->pathSamples[j][k]; //may be replaced below
//if(verts[oldSample]->sample && k != 0 && k != verts[ samples[i] ]->pathSamples[j].size() - 1)cout<<oldSample<<"intermediate is a sample\n";//prints and it's valid so no problem
					if (verts[oldSample]->sample)
						continue; //k=0 and k=last is definitely sample; an intermediate k might also be made sample during FPS/FPSroi() so don't update any existing samples (no replacement)
					size_t patchSize = verts[oldSample]->patch.size();
					if (patchSize == 0)
						cout << "WARNING: empty patch for " << oldSample << "\n\n";
					double maxGC = -INF;
					for (size_t m = 0; m < patchSize; m++) //1st condition below is necessary 'cos pathVert=false below helps with the evenly-spacing and the 1st condition makes me care about it
						if (verts[ verts[oldSample]->patch[m] ]->pathVert && verts[ verts[oldSample]->patch[m] ]->locMaxGCrelaxed && verts[ verts[oldSample]->patch[m] ]->gc > maxGC //max-curvature locMaxGC vertex is sought
							&& ! duplicateAhead(k, i, j, verts[oldSample]->patch[m])) //34662 34399 34706 35360 2309 --> 34399 could be replaced by 34706 ahead leading to duplicates (prevent that by this condition)
							maxGC = verts[ verts[oldSample]->patch[m] ]->gc, newSample = verts[oldSample]->patch[m];
					if (newSample == -1) //max-curvature locMaxGC vertex cannot be found so just take the max-curvature vertex as the winner (no such thing in replaceSample() 'cos sample replacement/selection is very selective)
					{
						maxGC = -INF;
						for (size_t m = 0; m < patchSize; m++)
							if (verts[ verts[oldSample]->patch[m] ]->pathVert && verts[ verts[oldSample]->patch[m] ]->gc > maxGC //max-curvature vertex is sought
								&& ! duplicateAhead(k, i, j, verts[oldSample]->patch[m])) //34662 34399 34706 35360 2309 --> 34399 could be replaced by 34706 ahead leading to duplicates (prevent that by this condition)
								maxGC = verts[ verts[oldSample]->patch[m] ]->gc, newSample = verts[oldSample]->patch[m];
					}
					newSample = (newSample != -1 ? newSample : oldSample); //necessary if pathVert=false loop below is enabled; otherwise redundant 'cos newSample cannot be -1 then
					verts[ samples[i] ]->pathSamples[j][k] = newSample; //if newSample = oldSample this assignment is redundant and patchAndDescriptors() returns immediately
					patchAndDescriptors(newSample); //if (i==0&&j==1)cout << oldSample << " --> " << newSample << endl;
					//to make nonuniform sampling more uniform, i.e., prevent too close newSamples, i mark the patch of newSample as unavailable for the next sample selection; to do so pathVert=false is sufficient
					for (size_t m = 0; m < verts[newSample]->patch.size(); m++) //too close samples may still validly arise, e.g., s+1 has newSample=-1 (all patch unavailable) and rolls back to oldSample which may be too close to previous sample s
						verts[ verts[newSample]->patch[m] ]->pathVert = false; //now it cannot be selected by an upcoming maxGC loop iteration due to the first condition
				}
			}
//if(id==0&&samples.size()>20){cout<<v1<<"to"<<samples[idx]<<": pathSamples after: "; for (size_t k = 0; k < verts[v1]->pathSamples[idx].size(); k++) cout << verts[v1]->pathSamples[idx][k] << " "; cout << "\n\n\n\n";}
	cout << "nonuniform; done in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs due to " << (eucPatching ? "euclidean" : "geodesic") << " patching\n";
}

bool Mesh::duplicateAhead(size_t k, size_t i, size_t j, int candidate)
{
	//returns true if candidate exists in the i.pathSamples[j] subarray i send here

	//34662 34399 34706 35360 2309 --> 34399 could be replaced by 34706 ahead leading to duplicates (this function helps me prevent that)
	for (size_t l = k + 1; l < verts[ samples[i] ]->pathSamples[j].size(); l++)
		if (verts[ samples[i] ]->pathSamples[j][l] == candidate)
			return true;
	return false;
}

void Mesh::pathSamplingCoarse(size_t n, bool print)
{
	//same as pathSamplingCoarseSlow() except this efficient version avoids dijkstraShortestPaths() which creates v.prev values for each path from scratch (uses the saved pathVerts vector instead)

	cout << "path sampling by coarseMap" << (print ? ".. " : " (for the upcoming pathSampling)\n");
	clock_t t = clock();
	for (size_t v = 0; v < verts.size(); v++)
		verts[v]->area = verts[v]->patchGC = -1.0, verts[v]->patch.clear(), verts[v]->pathSamples.clear(); //will be made different than -1 below for pathSamples only
	for (size_t i = 0; i < coarseMap.size(); i += 2) //even and odd idxs mesh1 and mesh2 side, respectively
	{
		int cSample = coarseMap[i + (mesh2 ? 1 : 0)]; //coarseMap sample is ready
		//to get the path from this coarseMap sample to each sample[j], use the precomputed cSample.pathVerts[ samples[j] ] and avoid the costly (memorywise costly when nDenseSamples >= 250) dijkstra recomputation that sets up the required v.prev values for the current path (includes start (cSample) and end (samples[j]) vertices)
		for (size_t j = 0; j < samples.size(); j++)
		{
			vector< int > pathSamples, pathVerts = verts[cSample]->pathVerts[j]; //all vertices on the path going from cSample to the jth sample are already stored in cSamples.pathVerts[j]
			double factor = verts[ samples[j] ]->geodesic[cSample] / n; //verts[cSample]->geodesic[ samples[j] ] is normalized and hence incompatible w/ verts[ samples[j] ]->geodesic[a] below
			int k = 0, a = pathVerts[k++]; //from a to a.prev, and then a.prev becomes the new a
			pathSamples.push_back(a); //includes end vertex (a=samples[j])
			while (a != cSample)
			{
				a = pathVerts[k++]; //advance to the next path vertex
				//nth sample selection (whenever i exceed the currentSize*factor for the first time, i hit the next desired uniformly distributed nth sample on the path)
				if (verts[ samples[j] ]->geodesic[a] > pathSamples.size()*factor && a != cSample) //distance to samples[j]; 2nd condition: cSample'll be added after the loop so not here
					pathSamples.push_back(a);
			}
			if (a != samples[j]) //a = samples[j] happens when i=j, i.e., path from sample to itself must only contain that sample itself (duplication occurs w/o this condition 'cos samples[j] already added above)
				pathSamples.push_back(a); //a = cSample for sure
			verts[cSample]->pathSamples.push_back(pathSamples); //samples on the path going from the ith coarseMap sample to the jth sample are given by cSample.pathSamples[j]
			//for printing only
			if (pathSamples.size() > n+1)
			{
				cout << "WARNING: excessive # of uniform path verts: " << pathSamples.size() << " > " << n+1 << " from " << cSample << " to " << samples[j] << "\n";
//				for (size_t k = 0; k < pathSamples.size(); k++)
//					cout << pathSamples[k] << "\t"; cout << "\n\n";
			}
//			else
//			{
//				cout << pathSamples.size() << " <= " << n+1 << " from " << cSample << " to " << samples[j] << "\t\t" << pathVerts.size() << endl;
//				for (size_t k = 0; k < pathSamples.size(); k++)
//					cout << pathSamples[k] << "ggggg\t"; cout << "\n\n";
//			}

			//else if (pathSamples.size() == n+1) nPathsInUse++; //else nTooShort++; //long-enough path has provided desired amount of uniform verts				
			//define patches around each path samples
			for (size_t z = 0; z < pathSamples.size(); z++)
				patchAndDescriptors(pathSamples[z]);
			//store samples from this sample (currently samples[j]) to each coarseMap sample (currently coarseMap[i]) separately to use in denseMap().pathSampling() later
//			if (i == 0) //using only 1 coarseMap sample is better for symmetric flip handling assuming that coarseMap does not have symmetric flip (disable this line to use all coarseMap samples which may fix the symmetric flip in the coarseMap)
				verts[ samples[j] ]->samplesToCoarse.push_back(pathSamples); //samplesToCoarse[0..m-1][0] is samples[j], [1] is next sample, .., [last] is cSample=coarseMap[i]
		}
	}
	if (print)
		cout << "uniform; done in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs due to " << (eucPatching ? "euclidean" : "geodesic") << " patching\n";
}

/*void Mesh::pathSamplingCoarseSlow(size_t n, bool print)
{
	//same as pathSampling() except paths run from coarseMap samples to other dense samples uniformly (no nonuniform option here to simplify/shorten coding; not that critical at this denseMap() level anyway)

	cout << "path sampling by coarseMap" << (print ? ".. " : " (for the upcoming pathSampling)\n");
	clock_t t = clock();
	for (size_t v = 0; v < verts.size(); v++)
		verts[v]->area = verts[v]->patchGC = -1.0, verts[v]->patch.clear(), verts[v]->pathSamples.clear(); //will be made different than -1 below for pathSamples only
	for (size_t i = 0; i < coarseMap.size(); i += 2) //even and odd idxs mesh1 and mesh2 side, respectively
	{
		int cSample = coarseMap[i + (mesh2 ? 1 : 0)]; //coarseMap sample is ready
		//to get the path, i need to re-set the v.prev values; so re-compute the dijkstra (to enter dijkstra, geodesic must be empty)
		verts[cSample]->geodesic.clear(), dijkstraShortestPaths(cSample);
		//now get the path that goes from any other sample to this coarseMap sample using the prev's set during dijkstra above (includes start (cSample) and end (samples[j]) vertices)
		for (size_t j = 0; j < samples.size(); j++)
		{
			vector< int > pathSamples; //path samples local to this loop
			double factor = verts[cSample]->geodesic[ samples[j] ] / n; //, minDiff = INF;
			int a = samples[j]; //from a to a.prev, and then a.prev becomes the new a
			pathSamples.push_back(a); //includes end vertex (a=samples[j])
			while (a != cSample)
			{
				if (verts[a]->prev == -1)
					cout << "WARNING: prev of " << a << " is not defined for " << cSample << " to " << samples[j] << "\n\n";
				a = verts[a]->prev;				
				//nth sample selection (whenever i exceed the k*factor for the first time, i hit the next desired uniformly distributed nth sample on the path)
				if (verts[ samples[j] ]->geodesic[a] > pathSamples.size()*factor && a != cSample) //distance to samples[j]; 2nd condition: cSample'll be added after the loop so not here
					pathSamples.push_back(a);
			}
			if (a != samples[j]) //a = samples[j] happens when i=j, i.e., path from sample to itself must only contain that sample itself (duplication occurs w/o this condition 'cos samples[j] already added above)
				pathSamples.push_back(a); //a = cSample for sure
			verts[cSample]->pathSamples.push_back(pathSamples); //samples on the path going from the ith coarseMap sample to the jth sample are given by cSample.pathSamples[j]
			//for printing only
			if (pathSamples.size() > n+1)
				cout << "WARNING: excessive # of uniform path verts: " << pathSamples.size() << " > " << n << "\n\n";
			//else if (pathSamples.size() == n+1) nPathsInUse++; //else nTooShort++; //long-enough path has provided desired amount of uniform verts				
			//define patches around each path samples
			for (size_t z = 0; z < pathSamples.size(); z++)
				patchAndDescriptors(pathSamples[z]);
			//store samples from this sample (currently samples[j]) to each coarseMap sample (currently coarseMap[i]) separately to use in denseMap().pathSampling() later
//			if (i == 0) //using only 1 coarseMap sample is better for symmetric flip handling assuming that coarseMap does not have symmetric flip (disable this line to use all coarseMap samples which may fix the symmetric flip in the coarseMap)
				verts[ samples[j] ]->samplesToCoarse.push_back(pathSamples); //samplesToCoarse[0..m-1][0] is samples[j], [1] is next sample, .., [last] is cSample=coarseMap[i]
		}
	}
	if (print)
		cout << "uniform; done in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs due to " << (eucPatching ? "euclidean" : "geodesic") << " patching\n";
}*/

void Mesh::patchAndDescriptors(int sourceVert)
{
	//sets the sourceVert.patch[] and its patch-based descriptors where sourceVert is a path sample

	if (! verts[sourceVert]->patch.empty())
		return; //don't recompute if patch was already computed/filled
	if (eucPatching) //for efficiency only
	{
		for (size_t v = 0; v < verts.size(); v++)
			if (distanceBetween(verts[sourceVert]->coords, verts[v]->coords) <= radius) //unnormalized euc distance vs. unnormalized radius
				verts[sourceVert]->patch.push_back(v); //sourceVert is included in its own neighborhood/patch
	}
	else //////////// modified Dijkstra (slower than modified BFS approximation: 4secs vs. 0.6secs) ////////////////////
	{
		for (size_t v = 0; v < verts.size(); v++)
			verts[v]->dSpanned = (v == sourceVert ? 0.0 : INF), verts[v]->marked = false;
		while (true)
		{
			//extract vert w/ min dSpanned value (no need to fill/maintain a heap 'cos this modified dijkstra will break early due to radius)
			int minDidx = -1, e, va;
			double minD = INF;
			for (size_t v = 0; v < verts.size(); v++)
				if (verts[v]->dSpanned < minD && ! verts[v]->marked)
					minD = verts[v]->dSpanned, minDidx = v;
			if (minDidx == -1)
				break; //no more verts w/ unknown shortest dist; so, we're done
			verts[sourceVert]->patch.push_back(minDidx); //sourceVert is included in its own neighborhood/patch (happens in the very first while loop iteration)
			verts[minDidx]->marked = true;
			if (minD > radius)  //far enough vert reached (unnormalized geo distance vs. unnormalized radius)
				break;
			//relax each edge incident to minDidx
			for (size_t ve = 0; ve < verts[minDidx]->edgeList.size(); ve++)
			{
				e = verts[minDidx]->edgeList[ve];
				va = (edges[e]->v1i == minDidx ? edges[e]->v2i : edges[e]->v1i);
				if (verts[minDidx]->dSpanned + edges[e]->length < verts[va]->dSpanned)
					verts[va]->dSpanned = verts[minDidx]->dSpanned + edges[e]->length; //relaxation (for length)
			}
			//relax each support edge, if any, incident to minDidx
			for (size_t ve = 0; ve < verts[minDidx]->sedgeList.size(); ve++)
			{
				e = verts[minDidx]->sedgeList[ve];
				va = (sedges[e]->v1i == minDidx ? sedges[e]->v2i : sedges[e]->v1i);
				if (verts[minDidx]->dSpanned + sedges[e]->length < verts[va]->dSpanned)
					verts[va]->dSpanned = verts[minDidx]->dSpanned + sedges[e]->length; //relaxation (for length)
			}//no support edge during patch areas 'cos it may rarely show disconnnections within patch due to support-based jumps*/
		} //end of while (true)		
	} //////////// modified Dijkstra (slower than modified BFS approximation: 4secs vs. 0.6secs) ends ////////////////////
	//////////// patch-based descriptors ////////////
	//patch area for the sourceVert is computed by adding 1/3 of each neighboring triangle area 'cos sourceVert.patchVert is responsible from 1/3 of the 3-vertex triangle
	verts[sourceVert]->area = verts[sourceVert]->patchGC = 0.0;
	size_t patchSize = verts[sourceVert]->patch.size();
	for (size_t v = 0; v < patchSize; v++)
		for (size_t tri = 0; tri < verts[ verts[sourceVert]->patch[v] ]->triNeighbors.size(); tri++)
			verts[sourceVert]->area += (tris[ verts[ verts[sourceVert]->patch[v] ]->triNeighbors[tri] ]->area / 3.0);
	//patch curvature/hksScalar for the sourceVert is computed by taking the avg of non-outlier gaussian curvatures/hksScalarOrgs of patch vertices; non-outlier is important 'cos at high-reso matches i have tiny triangles that produce gc>3000 while typical value is 0.9
	if (patchSize > 1) //size_t below cannot be negative
	{
		vector< double > curvs; //discard upper and lower 5% (or 10%) to eliminate outliers in the upcoming computation of patchGC
		curvs.reserve(patchSize); //preallocations to save tiny time (requesting more space than needed is also OK but i request the exact amount now)
		for (size_t v = 0; v < patchSize; v++)
			curvs.push_back(verts[ verts[sourceVert]->patch[v] ]->gc); //+ve -ve is important here as i'll discard from upper and lower ends
		size_t high = curvs.size() - 1, low = high / 5;//10;//20; //div by 10 (20) to discard 10% (5%); div by 5 to discard 20%
		high = high - 1 - low;
//low = 0, high = curvs.size(); //enable this to disable outlier elimination, i.e., all curvatures are used in patchGC computation
		sort(curvs.begin(), curvs.end()); //sorts in ascending order
		for (size_t i = low; i < high; i++)
		{
			verts[sourceVert]->patchGC += curvs[i];//abs(curvs[i]); //-ve and +ve gcs on the same patch yield ~0 patchGC which confuses w/ the planar patch's ~0 gc so use abs here to prevent that yield
			//verts[sourceVert]->patchGC += abs(curvs[i]); //-ve and +ve gcs on the same patch yield ~0 patchGC which confuses w/ the planar patch's ~0 gc so use abs here to prevent that yield
//if (sourceVert == 32865 || sourceVert == 3643) cout << curvs[i] << "\t";
		}
		verts[sourceVert]->patchGC /= (high - low);
//if (sourceVert == 32865 || sourceVert == 3643) cout << endl << high - low << " " << verts[sourceVert]->patchGC << "\n\n";
	}
	//////////// patch-based descriptors end ////////////
}

void Mesh::avgAreas(bool print)
{
	//sets avgPatchArea and avgRingArea for this mesh (for meaningful combination of patch areas and ring areas in qap forumlation, e.g., mapQualityPatch()); also sets avgPatchGC for info only

	//weightByGC=true updates v.area by multiplying it w/ patchGC (v.patchGC*v.area may be negative which does not affect the area>=0 condition below that makes me process only the pathSamples vertices (area gets -ve after that condition))
	double nAdds1 = 0.0, nAdds2 = 0.0;
	avgPatchArea = avgRingArea = avgPatchGC = 0.0;
	for (size_t v = 0; v < verts.size(); v++)
	{
		if (verts[v]->area >= 0.0) //area for nonsample verts is made -1 in pathSampling/pathSamplingCoarse() so ignore them here
			verts[v]->area *= (weightByGC ? verts[v]->patchGC : 1.0), avgPatchArea += verts[v]->area, nAdds1++;
		if (verts[v]->patchGC >= 0.0) //patchGC for nonsample verts is made -1 in pathSampling/pathSamplingCoarse() so ignore them here
			avgPatchGC += verts[v]->patchGC, nAdds2++;
	}
	avgPatchArea /= nAdds1, avgPatchGC /= nAdds2;
	for (size_t s = 0; s < samples.size(); s++)
		for (size_t a = 0; a < verts[ samples[s] ]->areasRing.size(); a++) //same size as areas[] which stores overlapping areas
			avgRingArea += verts[ samples[s] ]->areasRing[a]; //nonoverlapping rings
			//avgRingArea += verts[ samples[s] ]->areas[a]; //overlapping accumulated regions (not to be confused with the patches for path samples)
			//verts[ samples[s] ]->areasRing[a] *= (weightByGC ? verts[ samples[s] ]->patchGC : 1.0), avgRingArea += verts[ samples[s] ]->areasRing[a]; leads to bad initial coarse maps
	avgRingArea /= (samples.size() * verts[ samples[0] ]->areasRing.size()); //same for any sample so just use samples[0].areasRing.size() (areasRing.size() is same as areas.size() too)
	if (print)
		cout << "avg patch/ring area & gc: " << avgPatchArea << "/" << avgRingArea << " & " << avgPatchGC << endl;
}

void Mesh::fillQ(size_t pathSize, double areaWeight, bool descriptorAvailable)
{
	//fills the nxn Q2D matrix with pairwise geodesics between n samples (descriptor on diagonal) and nxnxpathSize Q3D where i,j entry is pathSize-dimensional, e.g.,
	//off-diagonal: set of pathSize many patch areas between i & j (instead of one length b/w i & j), diagonal: set of pathSize many areasRing as the pointwise descriptor of i
	
	//geodesic length based 2D Q matrix
	Q2D = new double*[ samples.size() ];
	for (size_t i = 0; i < samples.size(); i++)
	{
		Q2D[i] = new double[ samples.size() ];
		for (size_t j = 0; j < samples.size(); j++)
		{
			Q2D[i][j] = verts[ samples[i] ]->geodesic[ samples[j] ]; //sets to 0 if i==j
			if (i == j && descriptorAvailable) //0-value diagonal is replaced with a unary/pointwise descriptor value
				Q2D[i][j] = verts[ samples[i] ]->gc; //gaussian curvature descriptor (availability is not sufficient: descriptor must also be weighted appropriately like the areaWeight case below)
				//Q2D[i][j] = verts[ samples[i] ]->hksScalarOrg; //avg of the smaller half of time values for the HKS descriptor
//			cout << Q2D[i][j] << "\t";
		}
//		cout << endl;
	}
//	cout << "Q ready for mesh" << id << "\n" << (id > 0 ? "\n" : "");
	if (pathSize == 0)
		return; //no 3D Q matrix
	if (! samples.empty() && verts[ samples[0] ]->areasRing.size() > pathSize) //areasRing.size same for every sample so just use the first samples[0]
		cout << "WARNING: only the first " << pathSize << " entries in areasRing[] will be used (info loss)\n\n";
	//3D Q matrix that has pathSize-size arrays for the i-j'th entries (unlike 1-size geodesic scalars of Q2D matrix) and i-i'th entries (pointwise descriptor)
	Q = new double**[ samples.size() ];
	for (size_t i = 0; i < samples.size(); i++)
	{
		Q[i] = new double*[ samples.size() ];
		for (size_t j = 0; j < samples.size(); j++)
		{
			//size_t pathSizeLocal = verts[ samples[i] ]->pathSamples[j].size(); //pathSizeLocal < pathSize possible for short paths and pathSizeLocal < pathSize certainly for i=j case where pathSizeLocal=1
			Q[i][j] = new double[pathSize + 1]; //+1 for the last dimension storing pure geodesic length info from i to j
			for (size_t k = 0; k < pathSize; k++)
				if (i == j) //diagonal has the pointwise descriptor value
					//Q[i][j][k] = verts[ samples[i] ]->gc; //gaussian curvature descriptor
					//Q[i][j][k] = verts[ samples[i] ]->hksScalarOrg; //avg of the smaller half of time values for the HKS descriptor
					//Q[i][j][k] = (k < verts[ samples[i] ]->areas.size() ? areaWeight * verts[ samples[i] ]->areas[k] : 0.0); //overlapping/accumulated areas, i.e., previous areas included in the current area
					Q[i][j][k] = (k < verts[ samples[i] ]->areasRing.size() ? areaWeight * verts[ samples[i] ]->areasRing[k] : 0.0); //non-overlapping ring-based areas
				else if (k < verts[ samples[i] ]->pathSamples[j].size())
					Q[i][j][k] = verts[ verts[ samples[i] ]->pathSamples[j][k] ]->area;
				else
					Q[i][j][k] = 0.0; //pad w/ 0s 'till i have pathSize entries on this short path (that couldn't produce pathSize samples)
			//last dimension keeps geodesic length b/w i & j so that i can apply, e.g., maxGeoDiff filter during qapMap.computeQ2() (same way i do it in bruteForceMap.mapQualityVirtual())
			Q[i][j][pathSize] = verts[ samples[i] ]->geodesic[ samples[j] ]; //sets to 0 if i==j (Q[i][j][pathSize] = 0.0;//disable geoDiff > maxGeoDiff thresholding trick (setting maxGeoDiff=INF in main() has the same effect so this line is never used))
		}
	}
}

void Mesh::fillQdense(size_t pathSize, double areaWeight)
{
	//densified version of fillQ() where Qdense and Q2Ddense are filled by using ring areas (diagonals) and patch areas b/w samples (off-diagonals), and filled by lengthsToCoarseMap (diagonals) and lengthsBetweenSamples (off-diagonals)
	//no coarseMap at all in this function (except Q2Ddense's unary part: using 0 (hard to weight gc w/ geodesics) in fillQ() but here i've coarseMap info to be utilized)

	size_t offset = (mesh2 ? 1 : 0), cSize = coarseMap.size() / 2;
	//geodesic length based 2D Q matrix
	Q2Ddense = new double*[ samples.size() ];
	for (size_t i = 0; i < samples.size(); i++)
	{
		Q2Ddense[i] = new double[ samples.size() ];
		for (size_t j = 0; j < samples.size(); j++)
		{
			Q2Ddense[i][j] = verts[ samples[i] ]->geodesic[ samples[j] ]; //sets to 0 if i==j
			if (i == j) //0-value diagonal is replaced with a unary/pointwise descriptor value (geodesic to the frist coarseMap match; might have used average to all coarseMap matches as well)
				Q2Ddense[i][j] = (cSize > 0 ? verts[ samples[i] ]->geodesic[ coarseMap[0 + offset] ] : 0); //geodesic length to coarseMap (no such thing in mapQualityGeo() so expect different prints for qucoop objective and mapQualityGeo() distortion)
				//Q2Ddense[i][j] = 0.0; //same prints for qucoop objective and mapQualityGeo() distortion
		}
	}
	if (pathSize == 0)
		return; //no 3D Q matrix
	if (! samples.empty() && verts[ samples[0] ]->areasRing.size() > pathSize) //areasRing.size same for every sample so just use the first samples[0]
		cout << "WARNING: only the first " << pathSize << " entries in areasRing[] will be used (info loss)\n\n";	
	//3D Q matrix that has pathSize-size arrays for the i-j'th entries (unlike 1-size geodesic scalars of Q2D matrix) and i-i'th entries (pointwise descriptor)
	Qdense = new double**[ samples.size() ];
	for (size_t i = 0; i < samples.size(); i++)
	{
		Qdense[i] = new double*[ samples.size() ];
		for (size_t j = 0; j < samples.size(); j++)
		{
			Qdense[i][j] = new double[pathSize + 1 + cSize]; //+1 for the extra dimension storing pure geodesic length info from i to j and +cSize for the extra dimensions storing lengths from j to each coarseMap entries
			for (size_t k = 0; k < pathSize; k++)
				if (i == j) //diagonal has the pointwise descriptor value
					//Qdense[i][j][k] = (k < verts[ samples[i] ]->areas.size() ? areaWeight * verts[ samples[i] ]->areas[k] : 0.0); //overlapping/accumulated areas, i.e., previous areas included in the current area
					Qdense[i][j][k] = (k < verts[ samples[i] ]->areasRing.size() ? areaWeight * verts[ samples[i] ]->areasRing[k] : 0.0); //non-overlapping ring-based areas
				else if (k < verts[ samples[i] ]->pathSamples[j].size())
					Qdense[i][j][k] = verts[ verts[ samples[i] ]->pathSamples[j][k] ]->area;
				else
					Qdense[i][j][k] = 0.0; //pad w/ 0s 'till i have pathSize entries on this short path (that couldn't produce pathSize samples)
			//extra dimension keeps geodesic length b/w i & j so that i can apply, e.g., maxGeoDiff filter during qapMap.computeQ2() (same way i do it in bruteForceMap.mapQualityVirtual())
			Qdense[i][j][pathSize] = verts[ samples[i] ]->geodesic[ samples[j] ]; //sets to 0 if i==j
			//extra dimensions to keep geodesic lengths b/w j & coarseMap entries
			for (size_t e = 0; e < cSize; e++)
				Qdense[i][j][pathSize + 1 + e] = verts[ samples[j] ]->geodesic[ coarseMap[2*e + offset] ];
		}
	}
}

void Mesh::fillQdenseC2F(size_t pathSize, double areaWeight)
{
	//same as fillQdense() except this one uses the current set of samplesC2F (subset of global samples[]), e.g., Q2Ddense[0][1] = samples[ samplesC2F[0] ]->geodesic[ samples[ samplesC2F[1] ] ] instead of samples[0]->geodesic[ samples[1] ]

	size_t offset = (mesh2 ? 1 : 0), cSize = coarseMap.size() / 2;
	//geodesic length based 2D Q matrix
	Q2Ddense = new double*[ samplesC2F.size() ]; //note that size is based on subset samplesC2F[], not the global samples[]; global samples[] will be used to get the predefined correct rhs values for the unary/pairwise terms below
	for (size_t i = 0; i < samplesC2F.size(); i++)
	{
		Q2Ddense[i] = new double[ samplesC2F.size() ];
		for (size_t j = 0; j < samplesC2F.size(); j++)
		{
			Q2Ddense[i][j] = verts[ samplesC2F[i] ]->geodesic[ samplesC2F[j] ]; //sets to 0 if i==j
			if (i == j) //0-value diagonal is replaced with a unary/pointwise descriptor value (geodesic to the frist coarseMap match; might have used average to all coarseMap matches as well)
				Q2Ddense[i][j] = (cSize > 0 ? verts[ samplesC2F[i] ]->geodesic[ coarseMap[0 + offset] ] : 0); //geodesic length to coarseMap (no such thing in mapQualityGeo() so expect different prints for qucoop objective and mapQualityGeo() distortion)
				//Q2Ddense[i][j] = 0.0; //same prints for qucoop objective and mapQualityGeo() distortion
		}
	}
	if (pathSize == 0)
		return; //no 3D Q matrix
	if (! samplesC2F.empty() && verts[ samplesC2F[0] ]->areasRing.size() > pathSize) //areasRing.size same for every sample so just use the first samples[0]
		cout << "WARNING: only the first " << pathSize << " entries in areasRing[] will be used (info loss)\n\n";	
	//3D Q matrix that has pathSize-size arrays for the i-j'th entries (unlike 1-size geodesic scalars of Q2D matrix) and i-i'th entries (pointwise descriptor)
	Qdense = new double**[ samplesC2F.size() ];
	for (size_t i = 0; i < samplesC2F.size(); i++)
	{
		Qdense[i] = new double*[ samplesC2F.size() ];
		for (size_t j = 0; j < samplesC2F.size(); j++)
		{
			Qdense[i][j] = new double[pathSize + 1 + cSize]; //+1 for the extra dimension storing pure geodesic length info from i to j and +cSize for the extra dimensions storing lengths from j to each coarseMap entries
			for (size_t k = 0; k < pathSize; k++)
				if (i == j) //diagonal has the pointwise descriptor value
					//Qdense[i][j][k] = (k < verts[ samplesC2F[i] ]->areas.size() ? areaWeight * verts[ samplesC2F[i] ]->areas[k] : 0.0); //overlapping/accumulated areas, i.e., previous areas included in the current area
					Qdense[i][j][k] = (k < verts[ samplesC2F[i] ]->areasRing.size() ? areaWeight * verts[ samplesC2F[i] ]->areasRing[k] : 0.0); //non-overlapping ring-based areas
				else if (k < verts[ samplesC2F[i] ]->pathSamples[j].size())
					Qdense[i][j][k] = verts[ verts[ samplesC2F[i] ]->pathSamples[j][k] ]->area;
				else
					Qdense[i][j][k] = 0.0; //pad w/ 0s 'till i have pathSize entries on this short path (that couldn't produce pathSize samples)
			//extra dimension keeps geodesic length b/w i & j so that i can apply, e.g., maxGeoDiff filter during qapMap.computeQ2() (same way i do it in bruteForceMap.mapQualityVirtual())
			Qdense[i][j][pathSize] = verts[ samplesC2F[i] ]->geodesic[ samplesC2F[j] ]; //sets to 0 if i==j
			//extra dimensions to keep geodesic lengths b/w j & coarseMap entries
			for (size_t e = 0; e < cSize; e++)
				Qdense[i][j][pathSize + 1 + e] = verts[ samplesC2F[j] ]->geodesic[ coarseMap[2*e + offset] ];
		}
	}
}

int Mesh::extremeVert(bool curvatureBased)
{
	//returns idx of an extreme vertex, e.g., on finger or toe

	//vertex that is farthest away from arbitraryVert is definitely an extreme, where arbitraryVert can be anywhere on mesh (center, tips, head, ..)
	int arbitraryVert = 0, result = -1;
	while (verts[arbitraryVert]->edgeList.empty()) //skip isolated ones
		arbitraryVert++; //start w/ a nonisolated arbitrary vertex
	dijkstraShortestPaths(arbitraryVert); //arbitraryVert.geodesic[] is now ready, geodesic01 is bad due to unhealthy min/maxGeoDist values
	double maxDist = -INF;
	for (size_t s = 0; s < verts.size(); s++) //fps selects from amongst all vertices
		if (! verts[s]->sample && verts[arbitraryVert]->geodesic[s] > maxDist && ! verts[s]->edgeList.empty()) //3rd condition ensures isolated verts are not selected as extremes
		{
			if (fpsFromGT && verts[s]->gtIdx == -1) //for fair evaluation select amongs gt-equipped vertices
				continue;
			maxDist = verts[arbitraryVert]->geodesic[s];
			result = s;
		}
	verts[arbitraryVert]->geodesic.clear(); //redundant 'cos clear does not release or deallocate memory it just sets the size to 0 (useful for empty() tests which is not the case for arbitraryVert)
	//below is adapted from replaceSample() in order to replace result with an appropriate max-curvature locMaxGC vertex
	dijkstraShortestPaths(result);
	if (curvatureBased)
	{
		double maxGC = -INF;
		int replacement = -1;
		for (size_t v = 0; v < verts.size(); v++)
		{
			if (verts[v]->edgeList.empty() || //isolated verts should never be samples (even if they'd have dSapnning=INF values); they occur on microholes/holes/viewing classes of shrec11 dataset
				! verts[v]->locMaxGC) //only local maxima of gaussian curvatures can be the replacements (locMaxGCrelaxed not needed here 'cos this selects a single (very first) sample so be very selective/unrelaxed)
				continue;
			if (fpsFromGT && verts[v]->gtIdx == -1) //for fair evaluation select amongs gt-equipped vertices
				continue;
			if (verts[result]->geodesic[v] < radius * 1.5 && verts[v]->gc > maxGC)
				maxGC = verts[v]->gc, replacement = v;
		}
		if (replacement != -1 && replacement != result) //-1 means no replacement happened
			dijkstraShortestPaths(replacement), verts[result]->geodesic.clear(), result = replacement, cout << "extreme replaced..";
	}
	verts[result]->sample = true, samples.push_back(result);
	for (size_t v = 0; v < verts.size(); v++)
		if (verts[ samples[0] ]->geodesic[v] < radius) //extremeVert() creaets the very frist sample for sure hence use samples[0] here
			verts[v]->marked = true; //marked=true means v is in the vicinity of an existing FPS sample (for replaceSample() only; FPS() is already evenly-spaced (curvatureBased slightly distracts it but anyway) so this marked is for replaceSample() that has no evenly-spacing respect)
	return result;
}

void Mesh::curvaturesAll()
{
	//computes different types of curvatures for each vertex

	gaussianCurvatures();
//	meanCurvaturesApprox();
////	  meanCurvatures();
//	 principalCurvatures();
}

void Mesh::gaussianCurvatures()
{
	//computes Gaussian curvature at each vertex via 2pi - sumOfAngle div by area (saddle verts have sumOfAngles>2pi, i.e. negative gaussian curvature)

	int method = 1;//0;//2; //0: angle-deficit only, 1: divide deficit by approx voronoi area (On Surface Normal and Gaussian Curvature eq. 3.4), 2: divide deficit by exact voronoi area (Discrete Diff Geom Operators eq. 9 or Approximation of Gaussian Curvature by the Angular Defect), 3/4: variation1/2 of 2
	//(Discrete Schemes for Gaussian Curvature and Their Convergence eqs2/3); v.gc = 2pi - incident angle sum (method=0) is the angle-deficit formula; divide this by a suitable area to get methods1 & 2:
	//essentially, method=0 works for non-obtuse (no angle > 90) triangulations; divide deficit by areaSumMixed (method=2) or by simpler/approx areaSum/3 (method=1) to account for obtuse triangulations;
	//simple method=1 is recommended 'cos Discrete Schemes for Gaussian Curvature and Their Convergence paper states that all schemes converge only at degree-6 verts so no need to go crazy for the best scheme; voronoi area: region that is closer to x than to its neighbors
	//A surface has +ve curvature at x if the surface lives entirely on 1 side of the tangent plane at x (sphere pole), -ve curvature if it's saddle-shaped, 0 curvature if planar (see mphitchman.com/geometry/section7-1.html); x is called spherical, euclidean, saddle by SVG geodesic paper ying 2013
//method=2 not recommended 'cos using fixed areaVoro = 1.0 below

	double twoPI = 2*PI, PIover2 = PI/2.0, nl1[3], nl2[3], nl3[3], n = 0, gcMin = INF, gcMax = -INF;
	bool printOnce = true;
	avgGC = stdDevGC = 0.0;
	vector< double > curvs; //to sort curvatures in order to learn minPositiveCurvatureAllowed
	for (size_t v = 0; v < verts.size(); v++)
	{
		if (verts[v]->vertList1.empty()) //isolated vertex
		{
			verts[v]->gc = INF;
			if (printOnce)
				cout << "WARNING: gaussian curvature set to INF for a weird isolated vertex (there may be others; printing stopped)\n\n";
			printOnce = false;
			continue;
		}
		double angleSum = 0.0, angle, areaSum = 0.0, areaSumMixed2 = 0.0, areaSumMixed3 = 0.0, m4Extra = 0.0;//, areaVoro = (method == 2 ? voronoiArea(v) : 0); not now 'cos i need voronoi area of v in tn
		for (size_t tn = 0; tn < verts[v]->triNeighbors.size(); tn++)
		{
			int t = verts[v]->triNeighbors[tn];
			areaSum += tris[t]->area;
			//find 2 verts, v1 & v2, of t other than v, for angleSum
			int v1 = (tris[t]->v1i == v ? tris[t]->v2i : tris[t]->v1i), v2 = tris[t]->v1i;
			if (v2 == v1 || v2 == v)
				v2 = ((tris[t]->v2i != v1 && tris[t]->v2i != v) ? tris[t]->v2i : tris[t]->v3i);
			for (size_t c = 0; c < 3; c++)
			{
				nl1[c] = verts[v1]->coords[c] - verts[v]->coords[c]; //from v to neighbor vert v1
				nl2[c] = verts[v2]->coords[c] - verts[v]->coords[c]; //from v to neighbor vert v2
				nl3[c] = verts[v2]->coords[c] - verts[v1]->coords[c]; //from v1 to v2 (also needed to get the other angles of this t for its obtusity test)
			}
			angle = angleFromCotan(nl1, nl2); //first angle of t, i.e. the one at vertex v
			angleSum += angle;
			if (method == 2)
			{
double areaVoro = 1.0;//voronoiArea(v, t);
				double a2 = angleFromCotan(nl2, nl3); //2nd angle of t, i.e. the one at vertex v2
				double a3 = angleFromCotan(nl1, nl3); //3rd angle of t, i.e. the one at vertex v1
				if (angle >= PIover2 && angle <= PI) //angle of t at v is obtuse, i.e. between 90 and 180 degree
					areaSumMixed2 += tris[t]->area/2.0;
				else if ((a2 >= PIover2 && a2 <= PI) || (a3 >= PIover2 && a3 <= PI)) //t is obtuse, i.e. has an angle b/w 90 and 180 degree, but a1 is not that obtuse angle
					areaSumMixed2 += tris[t]->area/4.0;
				else //none of the angles between 90 & 180 degree, making t non-obtuse (= acute) so add Voronoi region (Fig 2b and eq 7 in desbrun) of v in t to areaSumMixed
areaSumMixed2 += areaVoro; //uses all 1-ring neighbors of v in the end of this tn loop
			}
			if (angle > 0.0 && method >= 3) //prevent div by zero by sin(angle) or cot(angle)=cos(angle)/sin(angle) below
			{
				double d1 = distanceBetween(verts[v]->coords, verts[v1]->coords), d2 = distanceBetween(verts[v]->coords, verts[v2]->coords), d3 = distanceBetween(verts[v1]->coords, verts[v2]->coords);
				if (method == 3)	
					areaSumMixed3 += ((d1*d2 - (cos(angle)/2.0) * (d1*d1 + d2*d2)) / (4.0*sin(angle)));
				else //if (method == 4)			
					m4Extra += (0.125*(cos(angle)/sin(angle))*d3*d3);
			}
		} //end of tn
		double totalAngle = (verts[v]->border ? PI : twoPI); //typicall 2PI but for border vertices total angle must be the half of it
		if (method == 0)
			verts[v]->gc = totalAngle - angleSum; //simply use angle-deficit as gaussian curvature
		else if (method == 1)
			verts[v]->gc = (totalAngle - angleSum) / (areaSum / 3.0); //divide angle-deficit by voronoi area approximation, i.e., take 1/3 of each triangle as v occupies 1/3 area (and other 2 verts occupy the remaining 2/3 area)
		else if (method == 2)
			verts[v]->gc = (totalAngle - angleSum) / areaSumMixed2; //divide angle-deficit by voronoi area (theoretically sound)
		else if (method == 3)
			verts[v]->gc = (totalAngle - angleSum) / areaSumMixed3; //divide angle-deficit by voronoi area approximation
		else //if (method == 4)
			verts[v]->gc = (totalAngle - angleSum) / ((areaSum / 2.0) - m4Extra); //divide angle-deficit by voronoi area approximation
		gcMax = (verts[v]->gc > gcMax ? verts[v]->gc : gcMax);
		gcMin = (verts[v]->gc < gcMin ? verts[v]->gc : gcMin);
		avgGC += verts[v]->gc; n++;
		curvs.push_back(verts[v]->gc);
	}
	avgGC /= n;
	sort(curvs.begin(), curvs.end()); //sorts in ascending order
	int nLocMaxGCs = localMaxGC(curvs);//, nLocMinNGCs = localMinNegGC(curvs);
	//curvature values' sorted quadrantHistogram: 13079 13 1 3 or 19885 0 0 2 so stdDev stays high due to the outliers in quadrants2-4 (13+1+3 or 2) which fails visualization (makes range too big); discard upper and lower 5% of curvature values to remedy this problem
	for (size_t i = 0; i < verts.size(); i++)
		if (verts[i]->gc != INF)
			stdDevGC += ((verts[i]->gc - avgGC) * (verts[i]->gc - avgGC));
	stdDevGC = sqrt(stdDevGC / n); //+= above counts n times for sure
//	if (printHeavy) //clipping values by discarding upper and lower 5% is a good option (not in use)
		cout << "min/avg/max Gaussian curvatures = " << gcMin << "/" << avgGC << "/" << gcMax << " w/ stdDev = " << stdDevGC << " (" << n << " verts w/ " << nLocMaxGCs << " local maxima)\n";//& " << nLocMinNGCs << " local maxima & minima)\n";
}

int Mesh::localMaxGC(const vector< double >& curvsSorted)
{
	//decides whether v is a local maxima of gaussian curvature values or not; locMaxGC pinpoints areas where a surface bends most intensely into a dome or pocket
	
	double nLocMaxGCs = 0, percentile = 0.9, minPositiveCurvatureAllowed = curvsSorted[ (int) (percentile * (curvsSorted.size() - 1)) ]; //top 10% (1.0-0.9 * 100 where precentile = 0.9) curvature vertices are eligible to be local maxima
	vector< int > candids; //local maxima gc vertices stored here are the candidates to survive the physical spacing enforcement below
	for (size_t v = 0; v < verts.size(); v++)
	{
		verts[v]->marked = false; //true means processed during physical spacing
		//threshold out flat/saddle/concave regions or weak peaks by specifying a minimum curvature threshold based on the mesh's curvature distribution, e.g., only top 10% curvature vertices are eligible to be local maxima
		if (verts[v]->gc < minPositiveCurvatureAllowed)
		{
			verts[v]->locMaxGC = false;
			continue; //2nd ring and minPositiveCurvatureAllowed enabled: 796/52565 and 102/5499 (recommended); only 1-ring and minPositiveCurvatureAllowed enabled: 1840/52565 and 186/5499
					  //796 and 102 reduce further to 31 and 37 (spacing=radius) or 14 and 17 (spacing=2radius) after the nice physical spacing enforcement action below
		}//*/

		verts[v]->locMaxGC = true;
		for (size_t va = 0; va < verts[v]->vertList1.size(); va++) //check 1-ring
			if (verts[ verts[v]->vertList1[va] ]->gc > verts[v]->gc || verts[v]->gc == INF) //2nd condition to discard weird isolated vertices
			{
				verts[v]->locMaxGC = false;
				break;
			}
		//5850/52565 local maxima or 643/5499 local maxima w/ 1-ring only which is too much; so check the 2nd ring below to decrease it to 2028/52565 and 241/5499 (minPositiveCurvatureAllowed disabled)
		if (verts[v]->locMaxGC) //check the 2nd ring if still local maxima
			for (size_t va = 0; va < verts[v]->vertList1.size(); va++)
				for (size_t va2 = 0; va2 < verts[ verts[v]->vertList1[va] ]->vertList1.size(); va2++)
					if (verts[ verts[ verts[v]->vertList1[va] ]->vertList1[va2] ]->gc > verts[v]->gc)
					{
						verts[v]->locMaxGC = false;
						break;
					}
/*		if (verts[v]->locMaxGC) //check the 3rd ring if still local maxima to reduce 2028/52565 and 241/5499 even further to 1098/52565 and 136/5499 (minPositiveCurvatureAllowed disabled)
			for (size_t va = 0; va < verts[v]->vertList1.size(); va++)
				for (size_t va2 = 0; va2 < verts[ verts[v]->vertList1[va] ]->vertList1.size(); va2++)
					for (size_t va3 = 0; va3 < verts[ verts[ verts[v]->vertList1[va] ]->vertList1[va2] ]->vertList1.size(); va3++)
					if (verts[ verts[ verts[ verts[v]->vertList1[va] ]->vertList1[va2] ]->vertList1[va3] ]->gc > verts[v]->gc)
					{
						verts[v]->locMaxGC = false;
						break;
					}*/
		if (verts[v]->locMaxGC)
			nLocMaxGCs++, candids.push_back(v); //to be sorted below
	}
	//greedily pick strongest (max gc) points and suppress/cancel nearby points, i.e., enforce minimum physical spacing between final locMaxGC vertices
	while (true)
	{
		double maxGC = -INF, spacing = radius * 1;//2; //recall radius = max(mesh1->avgEdgeLen, mesh2->avgEdgeLen) * 4.0 or 6.0
		int nextWinner = -1;
		for (size_t i = 0; i < candids.size(); i++)
			if (! verts[ candids[i] ]->marked && verts[ candids[i] ]->gc > maxGC)
				maxGC = verts[ candids[i] ]->gc, nextWinner = candids[i];
		if (nextWinner == -1)
			break; //all candidates are marked/processed by either surviving as nextWinner or being canceled so quit the loop
		for (size_t i = 0; i < candids.size(); i++)
			if (candids[i] != nextWinner && ! verts[ candids[i] ]->marked && distanceBetween(verts[nextWinner]->coords, verts[ candids[i] ]->coords) < spacing) //suppress/cancel this close locMaxGC vertex
				verts[ candids[i] ]->locMaxGC = false, nLocMaxGCs--, verts[ candids[i] ]->marked = true;
		verts[nextWinner]->marked = true; //indices processed in this iteration, by either being nextWinner or being canceled above, are marked
	}

	//locMaxGCrelaxed makes a relaxed/lighter trimming to accept all vertices that are at the local maxima of their 2-ring neighborhood, i.e., no minPositiveCurvatureAllowed and no physical spacing
	for (size_t v = 0; v < verts.size(); v++)
	{
		verts[v]->locMaxGCrelaxed = true; //used by replaceSample() and nonuniformPathSamplingHeuristic=2 only
		for (size_t va = 0; va < verts[v]->vertList1.size(); va++) //check 1-ring
			if (verts[ verts[v]->vertList1[va] ]->gc > verts[v]->gc || verts[v]->gc == INF) //2nd condition to discard weird isolated vertices
			{
				verts[v]->locMaxGCrelaxed = false;
				break;
			}
		if (verts[v]->locMaxGCrelaxed) //check the 2nd ring if still local maxima
			for (size_t va = 0; va < verts[v]->vertList1.size(); va++)
				for (size_t va2 = 0; va2 < verts[ verts[v]->vertList1[va] ]->vertList1.size(); va2++)
					if (verts[ verts[ verts[v]->vertList1[va] ]->vertList1[va2] ]->gc > verts[v]->gc)
					{
						verts[v]->locMaxGCrelaxed = false;
						break;
					}
	}

	return (int) nLocMaxGCs;
}

/*int Mesh::localMinNegGC(const vector< double >& curvsSorted)
{
	//decides whether v is a local maxima of negative gaussian curvature values or not; locMinNGC (N for negative) marks the deepest/most radical saddle areas on surface; note that locMin of gc can be 0, i.e., flat, which is an uninteresting point in my shape matching application
	//so go with locMinNGC; in other words, locMinGC (no N) could happen anywhere, e.g., if a surface has a deep, sharp bowl-like indentation, gc there is positive and if it transitions to a flat plane, that flat plane might be the local minimum
//almost the same result as localMaxGC() (tiny displacements on the ~same localMaxGC() samples) so disabling this function

	double nLocMinNGCs = 0, percentile = 0.9, maxNegativeCurvatureAllowed = curvsSorted[ (int) ((1.0 - percentile) * (curvsSorted.size() - 1)) ]; //bottom 10% curvature vertices are eligible to be local minima
	vector< int > candids; //local minima negative gc vertices stored here are the candidates to survive the physical spacing enforcement below
	for (size_t v = 0; v < verts.size(); v++)
	{
		verts[v]->marked = false; //true means processed during physical spacing
		//threshold out flat/polar regions or weak peaks by specifying a minimum curvature threshold based on the mesh's curvature distribution, e.g., only top 10% curvature vertices are eligible to be local minima
		if (verts[v]->gc > maxNegativeCurvatureAllowed)
		{
			verts[v]->locMinNGC = false;
			continue;
		}

		verts[v]->locMinNGC = true;
		for (size_t va = 0; va < verts[v]->vertList1.size(); va++) //check 1-ring
			if (verts[ verts[v]->vertList1[va] ]->gc < verts[v]->gc || verts[v]->gc == INF) //2nd condition to discard weird isolated vertices
			{
				verts[v]->locMinNGC = false;
				break;
			}
		if (verts[v]->locMinNGC) //check the 2nd ring if still local minima
			for (size_t va = 0; va < verts[v]->vertList1.size(); va++)
				for (size_t va2 = 0; va2 < verts[ verts[v]->vertList1[va] ]->vertList1.size(); va2++)
					if (verts[ verts[ verts[v]->vertList1[va] ]->vertList1[va2] ]->gc < verts[v]->gc)
					{
						verts[v]->locMinNGC = false;
						break;
					}
		if (verts[v]->locMinNGC)
			nLocMinNGCs++, candids.push_back(v); //to be sorted below
	}
	//greedily pick strongest (min gc) points and suppress/cancel nearby points, i.e., enforce minimum physical spacing between final locMinNGC vertices
	while (true)
	{
		double minGC = INF, spacing = radius * 1;//2; //recall radius = max(mesh1->avgEdgeLen, mesh2->avgEdgeLen) * 4.0 or 6.0
		int nextWinner = -1;
		for (size_t i = 0; i < candids.size(); i++)
			if (! verts[ candids[i] ]->marked && verts[ candids[i] ]->gc < minGC)
				minGC = verts[ candids[i] ]->gc, nextWinner = candids[i];
		if (nextWinner == -1)
			break; //all candidates are marked/processed by either surviving as nextWinner or being canceled so quit the loop
		for (size_t i = 0; i < candids.size(); i++)
			if (candids[i] != nextWinner && ! verts[ candids[i] ]->marked && distanceBetween(verts[nextWinner]->coords, verts[ candids[i] ]->coords) < spacing) //suppress/cancel this close locMinNGC vertex
				verts[ candids[i] ]->locMinNGC = false, nLocMinNGCs--, verts[ candids[i] ]->marked = true;
		verts[nextWinner]->marked = true; //indices processed in this iteration, by either being nextWinner or being canceled above, are marked
	}
	return (int) nLocMinNGCs;
}

double Mesh::voronoiArea(int v, int t)
{
	//returns voronoi area (region that si closer to v than its neighbors) of vertex v that lands in triangle t

	double areaVoro = 0.0, portion = 0.5; //for partial area case i use 2 triangles (t and its edge nghb tn) so return the portion for t, i.e., exclude tn
	bool fullArea = false; //false: partial voronoi area of v that lands in t (true: full voronoi area of v without any restriction)
	//look at each adjacent vertex of v
	for (size_t i = 0; i < verts[v]->vertList1.size(); i++)
	{
		int va = verts[v]->vertList1[i], vo1, vo2, t1, t2; //call-by-ref for the 2 opposite vertices of the edge vi-vj
		if (fullArea || tris[t]->v1i == va || tris[t]->v2i == va || tris[t]->v3i == va) //fullArea uses all neighboring triangles to v so if it's false i use only the 2 triangles that touch v-va edge
		{
			//va is a part of t so current edge v-va gives voronoi area that lands in t and t's edge neighbor tn; since 2 tris are involved divide the result based on the areas of t and tn, e.g., if same area then return half
			otherVertices(edgeSharedBy(v, va), vo1, vo2, t1, t2);	
			if (vo1 != -1)
			{
				double cotAlpha, cotBeta;
				cotAngles(v, va, vo1, vo2, cotAlpha, cotBeta); //dummy, dummy, dummy, dummy) fails so send 4 different dummies d1/2/3/4
				areaVoro += ((cotAlpha + cotBeta) * distanceBetween2(verts[v]->coords, verts[va]->coords));
				portion = (t == t1 ? tris[t1]->area / (tris[t1]->area + tris[t2]->area) : tris[t2]->area / (tris[t1]->area + tris[t2]->area));
			}
		}
	}
	if (! fullArea) //since 2 tris are involved divide the result based on the areas of t and tn, e.g., if same area then return half since portion=0.5 then
		areaVoro *= portion;
	//each v-va edge above considers 2 triangle so the same triangle is considered twice hence the areaVoro is doubled regardless of fullArea value; so div by 2 (not important 'cos doubled for all verts but anyway)
	return areaVoro / 2.0; //div by 2 to deal with the doubled area issue mentioned above
}

void Mesh::meanCurvaturesApprox()
{
	//approximates mean curvature via the edge-based solution in computergraphics.stackexchange.com/questions/1718

	double p12[3], n12[3], len, avgMC = 0, stdDev = 0, n = 0, mcMin = INF, mcMax = -INF; //difference in normals, projected along the edge, as a fraction of the length of the edge
	bool printOnce = true;
	for (size_t e = 0; e < edges.size(); e++)
	{
		for (size_t c = 0; c < 3; c++)
			p12[c] = verts[ edges[e]->v2i ]->coords[c] - verts[ edges[e]->v1i ]->coords[c], n12[c] = verts[ edges[e]->v2i ]->unormal[c] - verts[ edges[e]->v1i ]->unormal[c];
		len = p12[0]*p12[0] + p12[1]*p12[1] + p12[2]*p12[2]; //squared length hence no sqrt
		edges[e]->curvature = (p12[0]*n12[0] + p12[1]*n12[1] + p12[2]*n12[2]) / len;
	}
	//average of all absolute edge curvatures touching to v gives the curvature at v
	for (size_t v = 0; v < verts.size(); v++)
	{
		verts[v]->mca = 0.0; //mean curvature at v
		if (verts[v]->edgeList.empty())
		{
			verts[v]->mca = INF;
			if (printOnce)
				cout << "WARNING: mean curvature set to INF for a weird vertex (there may be others; printing stopped)\n\n"; //this also affects averaging etc so mesh should better not have an isolated vertex
			printOnce = false;
			continue;
		}
		for (size_t ve = 0; ve < verts[v]->edgeList.size(); ve++)
			verts[v]->mca += abs(edges[ verts[v]->edgeList[ve] ]->curvature);
		verts[v]->mca /= verts[v]->edgeList.size();
		mcMax = (verts[v]->mca > mcMax ? verts[v]->mca : mcMax);
		mcMin = (verts[v]->mca < mcMin ? verts[v]->mca : mcMin);
		avgMC += verts[v]->mca; n++;
	}
	avgMC /= n;
	//curvature values' sorted quadrantHistogram: 13027 44 7 18 or 19660 179 47 1 so stdDev stays high due to the outliers in quadrants2-4 which fails visualization (makes range too big); discard upper and lower 5% of curvature values to remedy this problem
	for (size_t i = 0; i < verts.size(); i++)
		if (verts[i]->mca != INF)
			stdDev += ((verts[i]->mca - avgMC) * (verts[i]->mca - avgMC));
	stdDev = sqrt(stdDev / n); //+= above counts n times for sure
	if (printHeavy) //clipping values by discarding upper and lower 5% is a good option (not in use)
		cout << "min/avg/max mean approx curvatures = " << mcMin << "/" << avgMC << "/" << mcMax << " w/ stdDev = " << stdDev << " (" << n << " verts)\n";
}

void Mesh::principalCurvatures()
{
	//computes min (kappa1) and max (kappa2) principal curvatures using their connections to gc and mc: gc = kappa1 * kappa2 and mc = (kappa1 + kappa2)/2

	double avgPr1 = 0.0, avgPr2 = 0.0, n = 0.0, stdDev1 = 0.0, stdDev2 = 0.0, pr1Min = INF, pr1Max = -INF, pr2Min = INF, pr2Max = -INF;
	bool printOnce = true;
	for (size_t v = 0; v < verts.size(); v++)
	{
		if (verts[v]->vertList1.empty()) //vert120 in minicooper is isolated so areaSum etc remains 0 for such verts
		{
			verts[v]->principal1 = verts[v]->principal2 = INF;
			if (printOnce)
				cout << "WARNING: principal curvatures set to INF for a weird vertex (there may be others; printing stopped)\n\n"; //this also affects averaging etc so mesh should better not have an isolated vertex
			printOnce = false;
			continue;
		}
		double delta = 
			
			
			
			//verts[v]->mc * verts[v]->mc - verts[v]->gc; //for eqs 10-11 in Discrete Differential Ops for Triangulated 2-Manifolds paper
			verts[v]->mca * verts[v]->mca - verts[v]->gc; //for eqs 10-11 in Discrete Differential Ops for Triangulated 2-Manifolds paper

		delta = (delta < 0.0 ? 0.0 : sqrt(delta)); //avoid numerical problems also perform square rooting
		verts[v]->principal2 = 
			
			
			//verts[v]->mc + delta;
			verts[v]->mca + delta;
		verts[v]->principal1 = 
			
			
			//verts[v]->mc - delta;
			verts[v]->mca - delta;
		pr1Max = (verts[v]->principal1 > pr1Max ? verts[v]->principal1 : pr1Max);
		pr2Max = (verts[v]->principal2 > pr2Max ? verts[v]->principal2 : pr2Max);
		pr1Min = (verts[v]->principal1 < pr1Min ? verts[v]->principal1 : pr1Min);
		pr2Min = (verts[v]->principal2 < pr2Min ? verts[v]->principal2 : pr2Min);
		avgPr1 += verts[v]->principal1; avgPr2 += verts[v]->principal2; n++;
	}
	avgPr1 /= n; avgPr2 /= n;
	for (size_t i = 0; i < verts.size(); i++)
	{
		if (verts[i]->principal1 != INF)
			stdDev1 += ((verts[i]->principal1 - avgPr1) * (verts[i]->principal1 - avgPr1));
		if (verts[i]->principal2 != INF)
			stdDev2 += ((verts[i]->principal2 - avgPr2) * (verts[i]->principal2 - avgPr2));
	}
	stdDev1 = sqrt(stdDev1 / n); stdDev2 = sqrt(stdDev2 / n); //+= above counts n times for sure
	if (printHeavy) //clipping values by discarding upper and lower 5% is a good option (not in use)
		cout << "min/avg/max principal1 curvatures = " << pr1Min << "/" << avgPr1 << "/" << pr1Max << " w/ stdDev = " << stdDev1 << " (" << n << " verts)\n",
		cout << "min/avg/max principal2 curvatures = " << pr2Min << "/" << avgPr2 << "/" << pr2Max << " w/ stdDev = " << stdDev2 << " (" << n << " verts)\n";
}*/

void Mesh::avgGeoDists()
{
	//computes average geodesic distance (agd) at each vertex (add distances to all other vertices then take the average)

	bool useGeoToSamples = false; //true: use geodesic distances to the samples (efficient), false: use Euclidean distance to all other vertices (no sampling trick)
	if (! samples.empty() && ! verts[ samples[0] ]->geodesic.empty())
		useGeoToSamples = true;
	double agdMin = INF, agdMax = -INF, avgAGD = 0.0;
	bool areaWeighted = false; //Mobius Transformations For Global Intrinsic Symmetry Analysis paper suggests using 1-ring-area/3 weighting for each distance being added (inconsistent w/ nonuniform triangulations, e.g., tosca-partial)
	for (size_t i = 0; i < verts.size(); i++) //if i is sample then i can compute dists to all verts as dijkstra already done for it; to keep everything consistent always use dists to all samples (# samples is 500 which is high so it's still very accurate)
	{
		verts[i]->agd = 0.0;
		if (useGeoToSamples)
		{
			for (size_t j = 0; j < samples.size(); j++) //sum of dists to all samples (for samples I already know geodesics)
				verts[i]->agd += (areaWeighted ? verts[ samples[j] ]->area / 3 : 1) * verts[ samples[j] ]->geodesic[i]; //geodesic defined for all samples for sure
			verts[i]->agd /= samples.size();
		}
		else
		{
			for (size_t j = 0; j < verts.size(); j++) //sum of dists to all vertices
				verts[i]->agd += distanceBetween(verts[j]->coords, verts[i]->coords); //euclidean distance defined everywhere
			verts[i]->agd /= verts.size();
		}
		agdMax = (verts[i]->agd > agdMax ? verts[i]->agd : agdMax);
		agdMin = (verts[i]->agd < agdMin ? verts[i]->agd : agdMin);
		avgAGD += verts[i]->agd;
	}
	avgAGD /= verts.size();
	if (printHeavy)
		cout << "min/avg/max agd = " << agdMin << "/" << avgAGD << "/" << agdMax << "\n";
}

void Mesh::areas()
{
	//computes areas of the set of triangles that are at a particular distance from samples[s] (areasRing[] stores non-overlapping ring shaped areas whereas areas[] stores overlapping/accumulated areas)

	cout << "area-related vectors per sample..";
	double minDist = 0.0, maxDist = radius; //particular geodesic distance is >=minDist and <maxDist (radius works as a stepSize here)
	for (size_t v = 0; v < verts.size(); v++)
		verts[v]->areas.clear(), verts[v]->areasRing.clear(); //for the second call (from denseMap())
	for (size_t s = 0; s < samples.size(); s++)
	{
		for (size_t a = 0; a < nRingAreas; a++)
		{
			verts[ samples[s] ]->areasRing.push_back(0.0), verts[ samples[s] ]->areas.push_back(0.0);
			for (size_t v = 0; v < verts.size(); v++)
			{
				if (verts[ samples[s] ]->geodesic[v] >= minDist && verts[ samples[s] ]->geodesic[v] < maxDist) //unnormalized geo distance vs. unnormalized radius-based maxDist
					for (size_t t = 0; t < verts[v]->triNeighbors.size(); t++)
						verts[ samples[s] ]->areasRing[a] += (tris[ verts[v]->triNeighbors[t] ]->area / 3.0); //patch area is computed by adding 1/3 of this related triangle area 'cos v is responsible from 1/3 of the 3-vertex triangle						
				if (verts[ samples[s] ]->geodesic[v] < maxDist) //unnormalized geo distance vs. unnormalized radius-based maxDist
					for (size_t t = 0; t < verts[v]->triNeighbors.size(); t++)
						verts[ samples[s] ]->areas[a] += (tris[ verts[v]->triNeighbors[t] ]->area / 3.0); //patch area is computed by adding 1/3 of this related triangle area 'cos v is responsible from 1/3 of the 3-vertex triangle						
			}
			minDist += radius, maxDist += radius; //advance distance thresholds to handle the next ring in the next iteration (minDist not used/advanced for areas[] hence it stores overlapping/accumulated areas)
		}
		minDist = 0.0, maxDist = radius; //reset for the next sample
	}
	cout << "done\n";
}

void Mesh::areasForDisplay(int i)
{
	//sets display-related binID values by computing areas of the set of triangles that are at a particular distant from samples[i]

//i = id == 0 ? 2 : 0;//must be consistent w/ prints in areas() above
	for (size_t v = 0; v < verts.size(); v++)
		verts[v]->binID = -1; //for visualization only
	double minDist = 0.0, maxDist = radius; //same parameters and structure as areas() but now i know the particular samples[i] as the endpoint of the first match (mesh1.m1 - mesh2.m1) (radius works as a stepSize here)
	for (size_t a = 0; a < nRingAreas; a++)
	{
		for (size_t v = 0; v < verts.size(); v++)
			if (verts[ samples[i] ]->geodesic[v] >= minDist && verts[ samples[i] ]->geodesic[v] < maxDist)
				verts[v]->binID = a; //binID set w.r.t. distances from samples[i] for visualization only
		minDist += radius, maxDist += radius; //advance distance thresholds to handle the next ring in the next iteration (minDist not used/advanced for areas[] hence it stores overlapping/accumulated areas)
	}
}

void Mesh::doubleBorderAreas()
{
	//doubles the areas/areasRing/area around the border samples (updates avgRingArea and avgPatchArea too)

	for (size_t s = 0; s < samples.size(); s++)
		for (size_t a = 0; a < verts[ samples[s] ]->areas.size(); a++)
			verts[ samples[s] ]->areasRing[a] *= (verts[ samples[s] ]->border ? 2 : 1), //double it for border vertices as the surface around them are missing/partial; doubled 'cos assuming symmetric misses (better alternative: exclude border samples during optimal map search)
			verts[ samples[s] ]->areas[a] *= (verts[ samples[s] ]->border ? 2 : 1);
	for (size_t v = 0; v < verts.size(); v++) //p.pathSamples[q][i] will access this area data so just update it over all vertices
		verts[v]->area *= (verts[v]->border ? 2 : 1);
	avgAreas(); //avgRingArea and avgPatchArea are updated as they're used by the areaWeight of the qap formulation
}

/*void Mesh::ROIold()
{
	//computes region of interest (ROI) which includes verts that are close to matched samples and distant to unmatched samples

	//mark v as ROI vertex if geodesic from matched sample s to v is < maxGeo for all s
	double maxGeo = -INF; //max geodesic distance b/w matched samples
	for (size_t i = 0; i < samples.size(); i++)
		if (verts[ samples[i] ]->matchIdx != -1) //matched sample
			for (size_t j = 0; j < samples.size(); j++)
				if (verts[ samples[j] ]->matchIdx != -1 && verts[ samples[i] ]->geodesic[ samples[j] ] > maxGeo)
					maxGeo = verts[ samples[i] ]->geodesic[ samples[j] ];
	for (size_t v = 0; v < verts.size(); v++)
	{
		verts[v]->roi = true; //a region of interest vertex
		for (size_t i = 0; i < samples.size(); i++)
			if (verts[ samples[i] ]->matchIdx != -1 && verts[ samples[i] ]->geodesic[v] >= maxGeo)
				verts[v]->roi = false, i = samples.size(); //latter assignment to break the inner loop
	}

	//unmark v (not roi anymore) that are at most s.geo[s']/2 apart from s' where s is the closest matched sample to the unmatched s'
	for (size_t i = 0; i < samples.size(); i++)
		if (verts[ samples[i] ]->matchIdx == -1) //unmatched sample s'
		{
			double minGeo = INF;
			for (size_t j = 0; j < samples.size(); j++)
				if (verts[ samples[j] ]->matchIdx != -1 && verts[ samples[i] ]->geodesic[ samples[j] ] < minGeo)
					minGeo = verts[ samples[i] ]->geodesic[ samples[j] ];//, s = samples[j]; //closest matched sample to the unmatched samples[i]
			maxGeo = minGeo / 2.0;
			//unmark based on maxGeo and the current unmatched samples[i]
			for (size_t v = 0; v < verts.size(); v++)
				if (verts[ samples[i] ]->geodesic[v] <= maxGeo)
					verts[v]->roi = false;
		}
	int nROIverts = 0;
	for (size_t v = 0; v < verts.size(); v++)
		nROIverts += verts[v]->roi;
	cout << nROIverts << " roi vertices\n";
}*/
double Mesh::ROI(const vector< int >& coarseMapFull, double upperBound)
{
	//computes region of interest (ROI) by taking the union of regions b/w each pair of matched samples
	
	clock_t t = clock();
	size_t offset = (mesh2 ? 1 : 0);
	bool showAndApplyProgress = true; //printing and also used to return early if the roiArea covers 95% of the totalArea; this large coverage cannot happen on complete models matching w/ partial models but upperBound can still help us return early on such cases
	double roiArea = 0.0, nFills = 0, coverPercent, coverPercent2;
	cout << (offset == 1 ? "\n" : "") << "roi..";
	for (size_t v = 0; v < verts.size(); v++)
		verts[v]->roi = verts[v]->roiBorder = verts[v]->marked = false; //roi=true if v is a region of interest vertex, roiBorder=true means at the border of the final roi, marked=true means used in roiArea computation (prevent reusage)
	for (size_t i = 0; i < coarseMapFull.size(); i += 2) //each 2 consecutive items make 1 match (even indices belong to mesh1, odd indices mesh2)
		for (size_t j = i + 2; j < coarseMapFull.size(); j += 2) //look for a different matched sample (pathSampling() requires j = 0 during sampling 'cos fillQ() uses that symmetric info but no such requirement here)
		{
//if (verts[ coarseMapFull[i + offset] ]->roi && verts[ coarseMapFull[j + offset] ]->roi) //reduces 20 secs to 4 secs in 52K-vertex with 10x10 coarseMapFull[] size but not perfectly robust (see north-south comment below) so disable this; robust coverPercent/2 brings this efficiency now
//continue; //for efficiency skip the coarseMapFull samples that are already inside the ROI (made roi by a previous ij iteration here); not that robust 'cos [i + offset] may be the south endpoint of a path covering north and [j + offset] be the north endpnt of a path covering south; north-south region uncovered now w/ this continue;

//if (! (i == 0 && j == 4))continue; else cout << "roi from coarseMapFull" << i/2 << " to " << j/2 << " shown\n\n\n\n\n"; //use and display a particular ROI
			ROIij(coarseMapFull[i + offset], coarseMapFull[j + offset]); //i.roiList updated
			for (size_t k = 0; k < verts[ coarseMapFull[i + offset] ]->roiList.size(); k++)
			{
				size_t v = verts[ coarseMapFull[i + offset] ]->roiList[k];				
				if (showAndApplyProgress && ! verts[v]->marked)
					for (size_t tri = 0; tri < verts[v]->triNeighbors.size(); tri++)
							roiArea += (tris[ verts[v]->triNeighbors[tri] ]->area / 3.0);
				verts[v]->roi = verts[v]->marked = true;
			}
			coverPercent = (roiArea / totalArea) * 100, coverPercent2 = (roiArea / upperBound) * 100; //the latter is the coverage based on the roiArea of the partial mesh2 model (no problem if mesh2 is not partial, i.e., complete matching scenario)
			verts[ coarseMapFull[i + offset] ]->roiList.clear(); //same coarseMapFull[i] will be partnered w/ the next coarseMapFull[j] so clear the existing roiList for efficiency only (redundant 'cos does not affect accuracy)
			if (showAndApplyProgress)
			{
				cout << roiArea << " / " << totalArea << " (" << coverPercent << "%) covered\n"; //58% 75% 88% 98% .. 99% on complete matching execution
				if (coverPercent > 95 || //regionGrowing below will cover the remaining 5% very efficiently so break early (not that critical if remained uncovered)
					coverPercent2 > 95) //early break on complete model as it reaches the coverage threshold based on the partial mesh2's roiArea
					i = j = coarseMapFull.size(); //break from this nested ij loop immediately
			}
		}
	nROIverts = 0;
	for (size_t v = 0; v < verts.size(); v++)
		if (verts[v]->roi)
		{
			nROIverts++;
			//patch area for the sourceVert is computed by adding 1/3 of each neighboring triangle area 'cos sourceVert.patchVert is responsible from 1/3 of the 3-vertex triangle
			if (! showAndApplyProgress) //already shown above so skip here
				for (size_t tri = 0; tri < verts[v]->triNeighbors.size(); tri++)
					roiArea += (tris[ verts[v]->triNeighbors[tri] ]->area / 3.0);
			for (size_t va = 0; va < verts[v]->vertList1.size(); va++)
				if (! verts[ verts[v]->vertList1[va] ]->roi)
					verts[v]->roiBorder = true; //used by FPSroi() in order to favor the selection of roiBorder vertices as dense samples and by regionGrowing() to fill small non-roi regions inside the main roi
		}
	//region growing to fill the small gaps on the current ROI set above
	for (size_t v = 0; v < verts.size(); v++)
		if (verts[v]->roiBorder)// && ! pfaust) //2nd condition to disable region growing on the pfaust dataset
		{
//break;//disable region growing
			int r = -1; //non-roi vertex adjacent to a roiBorder vertex
			for (size_t va = 0; va < verts[v]->vertList1.size(); va++)
				if (! verts[ verts[v]->vertList1[va] ]->roi)
					r = verts[v]->vertList1[va], va = verts[v]->vertList1.size(); //early break
			if (r != -1) //a previous v's r has already filled the small non-roi region so this current v.r cannot find a non-roi vertex adjacent to v, causing r=-1 validly
			{
				vector< int > region = regionGrowing(r, roiArea); //fill small non-roi regions inside the main roi (empty region returns here if it's leaked, i.e., covers too much area w.r.t. main roi)
				for (size_t i = 0; i < region.size(); i++)
				{
					nROIverts++;
					verts[ region[i] ]->roi = true;
					for (size_t tri = 0; tri < verts[ region[i] ]->triNeighbors.size(); tri++)
						roiArea += (tris[ verts[ region[i] ]->triNeighbors[tri] ]->area / 3.0);
				}
				nFills += region.size();
			}
		}//*/
	cout << nROIverts << " roi vertices (filled: " << nFills << ") w/ area " << roiArea << " in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs\n";
//for (size_t v = 0; v < verts.size(); v++) verts[v]->roi = false; //display a particular ROI given below
//size_t i = (id==0?0:15), j = (id==0?17:19); if(i<samples.size()&&j<samples.size()){ROIij(samples[i], samples[j]); for (size_t k = 0; k < verts[ samples[i] ]->roiList.size(); k++) verts[ verts[ samples[i] ]->roiList[k] ]->roi = true;}
	return roiArea;
}

vector< int > Mesh::regionGrowing(int sourceVert, double mainArea)
{
	//grows region, starting from sourceVert, via modified BFS that stops when a roiBorder vertex is hit

	double segArea = 0.0, tooLarge = mainArea / 10;//12.75
	queue< int > q; vector< int > region;
	q.push(sourceVert);
	for (size_t v = 0; v < verts.size(); v++)
		verts[v]->marked = (v != sourceVert ? false : true); //true means visited by modified BFS
	while (! q.empty())
	{
		int v = q.front();
		q.pop(), region.push_back(v);
		for (size_t tri = 0; tri < verts[v]->triNeighbors.size(); tri++)
			segArea += (tris[ verts[v]->triNeighbors[tri] ]->area / 3.0); //1 vertex makes 1/3 contribution to a given triangle area hence div by 3 (if v is adjacent to roi vert(s) then only 1/3 of this region tri'll be used which is wrong but not critical 'cos /10 threshold is not certain anyway)
		if (segArea > tooLarge) //segArea is too large so region growing is undone and an empty result is returned
		{
			region.clear();
			return region;
		}
		//look at each adjacent vertex of v
		for (size_t i = 0; i < verts[v]->vertList1.size(); i++)
		{
			int va = verts[v]->vertList1[i];
			if (! verts[va]->marked && ! verts[va]->roiBorder) //modified by the addition of the 2nd condition
			{
				verts[va]->marked = true;
				q.push(va);
			}
		}
	}
	return region;
}

void Mesh::ROIij(int vi, int vj)
{
	//computes region of interest (ROI) b/w vertices vi and vj (where vi and vj are samples w/ valid matches)

	//initial ROI as the region within roiWidth distance from the vi-vj geodesic path
	double roiWidth = 0.1;//0.3;//each v.geoToPath is divided by vi-vj path length so roiWidth becomes relative to the path length (and 0.3 makes sense no matter what the scale is)
	bool fuzzyScheme = true, useAllPathVerts = false; //use all path vertices as path samples (perfectly accurate but slow); fuzzyScheme faster than samplingScheme and gives similar ROIs (each ROI in 0.3secs vs. 1.5secs/6secs when n=20/100)
	if (fuzzyScheme)
	{
		//if condition uses fuzzy geodesic: how long is the shortest path from vi to vj through the diversion v, in comparison to the direct shortest path from vi to vj
		for (size_t v = 0; v < verts.size(); v++)
			if (abs(verts[vi]->geodesic[v] + verts[vj]->geodesic[v] - verts[vi]->geodesic[vj]) / verts[vi]->geodesic[vj] < roiWidth) //normalization so that roiWidth makes sense
				verts[vi]->roiList.push_back(v); //verts[vj]->roiList.push_back(v) is redundant so skipped (roiList assumed to be empty before the call)
	}
	else //samplingScheme
	{
		////////////// adapted from pathSampling() //////////////
		//to get the path, i need to re-set the v.prev values; so re-compute the dijkstra (to enter dijkstra, geodesic must be empty)
		verts[vi]->geodesic.clear(), dijkstraShortestPaths(vi);
		size_t n = (useAllPathVerts ? verts.size() : 20);//100; //using different # path samples than pathSampling() for the accuracy of ROI
		vector< int > pathSamples; //samples on the path going from sample vi to sample vj
		double factor = verts[vi]->geodesic[vj] / n;
		int a = vj; //from a to a.prev, and then a.prev becomes the new a
		pathSamples.push_back(a); //includes end vertex (a=vj)
		while (a != vi)
		{
			a = verts[a]->prev;				
			if (verts[vj]->geodesic[a] > pathSamples.size()*factor && a != vi)
				pathSamples.push_back(a);
		}
		if (a != vj) //happens when vi=vj, i.e., path from sample to itself must only contain that sample itself (duplication occurs w/o this condition 'cos vj already added above)
			pathSamples.push_back(a); //a = vi for sure
		////////////// adapted from pathSampling() ends //////////////					
		for (size_t p = 0; p < pathSamples.size(); p++)
			dijkstraShortestPaths(pathSamples[p]); //geoToPath below needs geodesic distances to path samples
		for (size_t v = 0; v < verts.size(); v++)
		{
			verts[v]->geoToPath = INF;
			for (size_t p = 0; p < pathSamples.size(); p++)
				if (verts[ pathSamples[p] ]->geodesic[v] < verts[v]->geoToPath)
					verts[v]->geoToPath = verts[ pathSamples[p] ]->geodesic[v];
			verts[v]->geoToPath /= verts[vi]->geodesic[vj]; //normalization so that roiWidth makes sense
			if (verts[v]->geoToPath < roiWidth)
				verts[vi]->roiList.push_back(v); //verts[vj]->roiList.push_back(v) is redundant so skipped (roiList assumed to be empty before the call)
		}
//cout << "roi by path: " <<vi<< "-" <<vj<< " used " << pathSamples.size() << " path samples\n";
	}
//return;//enable this to prevent excluding ROI part beyond vi and vj
	//filtering ROI as intersection of patches of i and j, i.e., all verts in at most i-j path length distance to i are in the patch of i
	vector< int > filteringRegion, finalRegion;
	for (size_t v = 0; v < verts.size(); v++)
		if (verts[vi]->geodesic[v] <= verts[vi]->geodesic[vj] + TINY && verts[vj]->geodesic[v] <= verts[vi]->geodesic[vj] + TINY) //51.7237 <= 51.7237 didn't return true for precision issues so use 51.7237 vs 51.7237+TINY
			filteringRegion.push_back(v);
	//intersect filtering ROI w/ the initial ROI to to exclude the portions of ROI that lie beyond/behind vi and vj
	for (size_t i = 0; i < filteringRegion.size(); i++)
		for (size_t j = 0; j < verts[vi]->roiList.size(); j++)
			if (filteringRegion[i] == verts[vi]->roiList[j])
				finalRegion.push_back(filteringRegion[i]);
	verts[vi]->roiList = finalRegion;
}

void Mesh::geoAll(vector< int > anchors, int nCloseAnchors)
{
	//sets the geoAll unit vector of all verts which stores all geodesic distances to anchor vertices

//nCloseAnchors = 10;
	bool makeUnit = true;//(pfaust ? false : true); //true: normalize geoAll[] so it has unit length (true/false good for shrec16/pfaust respectively)
	for (size_t v = 0; v < verts.size(); v++)
	{
		double vecLen = 0.0;	
		//dist from v to anchor is considered
		for (size_t y = 0; y < anchors.size(); y++)
		{
			double dist = verts[ anchors[y] ]->geodesic[v];
			verts[v]->geoAll.push_back(dist);
			vecLen += (dist * dist);
		}
		if (makeUnit)
		{
			vecLen = sqrt(vecLen);
			if (vecLen > 0)
				//geoAll appends many independent vectors (w/ commensurable quantities) so normalizing it to unit length makes no sense (no it's good 'cos this is a vector in feature space so get direction)
				for (size_t d = 0; d < verts[v]->geoAll.size(); d++) //normalize this vector so it has a unit length
					verts[v]->geoAll[d] /= vecLen;
				else
					cout << "WARNING: vertex produced a 0-vector geoAll\n\n";
		}
	}
	//keep only nCloseAnchors distances in geoAll that show distances from v to nCloseAnchors many anchor matches; this way i'm preventing tooFar and hence erroneous distances in v.geoAll
	if (nCloseAnchors > 0 && nCloseAnchors < (int) verts[0]->geoAll.size() / 2) //nonpositive is a special flag that tells me not to do this prevention
		for (size_t v = 0; v < verts.size(); v++)
		{	
			vector< double > tmpWinners;
			vector< int > tmpWinnersIDs;			
			for (int i = 0; i < nCloseAnchors; i++)
			{
				double minDist = INF;
				int minDistIdx = -1;
				for (size_t d = 0; d < verts[v]->geoAll.size(); d++)
				{
					//check whether the d'th entry is already inserted to tmpWinnersIDs
					bool alreadyInserted = false;
					for (size_t j = 0; j < tmpWinnersIDs.size(); j++)
						if (tmpWinnersIDs[j] == d)
							alreadyInserted = true;
					if (! alreadyInserted && verts[v]->geoAll[d] < minDist)
					{
						minDist = verts[v]->geoAll[d];
						minDistIdx = d;
					}					
				}
				tmpWinnersIDs.push_back(minDistIdx); //nCloseAnchors inserts are done to tmpWinnersIDs for this v
				tmpWinners.push_back(minDist);
			}
			verts[v]->geoAll.clear(); //refill this w/ nCloseAnchors many entries in tmpWinners
			double vecLen = 0.0;
			for (int i = 0; i < nCloseAnchors; i++)
			{
				verts[v]->geoAll.push_back(tmpWinners[i]);
				vecLen += (tmpWinners[i] * tmpWinners[i]);	
			}
			if (makeUnit)
			{
				vecLen = sqrt(vecLen);
				if (vecLen > 0)
					for (size_t d = 0; d < verts[v]->geoAll.size(); d++) //normalize this vector so it has a unit length
						verts[v]->geoAll[d] /= vecLen;
					else
						cout << "WARNING: vertex produced a 0-vector geoAll\n\n";
			}
		}
}

void Mesh::vertColors(int patchColoring, bool roiReady)
{
	//determines color values of verts; colormap is from paulbourke.net/miscellaneous/colourspace

	if (patchColoring == 1) //no descriptor coloring; just color patches of the samples on the path from samples[i] to samples[j] in red
	{
		bool onePathOnly = true; //draw i-j path only in red (no green i-k and white j-k paths on screen)
		for (size_t i = 0; i < verts.size(); i++)
			verts[i]->color = color; //default mesh color for the non-patch vertices
		float* patchColor = new float[3]; patchColor[0] = 1.0f, patchColor[1] = 0.0f, patchColor[2] = 0.0f; //float patchColor[3] = {1.0f, 0.0f, 0.0f}; does not work properly
		//sphered paths run from samples m1 to m2, m1 to m3, and m2 to m3 where m1 to m2 path also colors its patches via patchColor; note that mesh1.m1-mesh2.m1, mesh1.m2-mesh2.m2, mesh1.m3-mesh2.m3 are the 1st 3 matches
//if(id==0)m1=0,m2=2;else m1=0,m2=1;//see a particular path on mesh1/2 via manual update
//m1 = (mesh2 ? 3 : 29), m2 = (mesh2 ? 8 : 9); if (m1 >= samples.size() || m2 >= samples.size() || m3 >= samples.size()) cout << "wrong manual index\n", exit(0);
		size_t i = m1, j = m2, k = m3;
		for (size_t z = 0; ! verts[ samples[i] ]->pathSamples.empty() && z < verts[ samples[i] ]->pathSamples[j].size(); z++)
		{
//if (z == 5)break;
			for (size_t p = 0; p < verts[ verts[ samples[i] ]->pathSamples[j][z] ]->patch.size(); p++)
				verts[ verts[ verts[ samples[i] ]->pathSamples[j][z] ]->patch[p] ]->color = patchColor; //patch coloring
			sphereIdxs.push_back(verts[ samples[i] ]->pathSamples[j][z]); //patch center sphering
//cout << z << "\t" << verts[ verts[ samples[i] ]->pathSamples[j][z] ]->area << "area\t" << verts[ verts[ samples[i] ]->pathSamples[j][z] ]->patchGC << "gc\t" 
//	 << verts[ verts[ samples[i] ]->pathSamples[j][z] ]->patchHKSscalar << "hks\t" << verts[ verts[ samples[i] ]->pathSamples[j][z] ]->patchHKSscalar2 << "hks2\n";
		}
/*		cout << "patch areas shown on path from " << samples[i] << " to " << samples[j] << " (" << verts[ samples[i] ]->pathSamples[j].size() << " spheres) for mesh" << id << endl;
for (size_t h = 0; h < verts[ samples[i] ]->pathSamples[j].size(); h++) cout << verts[ verts[ samples[i] ]->pathSamples[j][h] ]->area << "a\t";cout<<"\n";
for (size_t h = 0; h < verts[ samples[i] ]->pathSamples[j].size(); h++) cout << verts[ verts[ samples[i] ]->pathSamples[j][h] ]->area * verts[ verts[ samples[i] ]->pathSamples[j][h] ]->patchGC << "a*gc\t";cout<<"\n\n\n\n\n";//*/
for (size_t v = 0; v < verts.size(); v++)verts[v]->marked=false; verts[ samples[i] ]->marked = verts[ samples[j] ]->marked = true;
//cout << i << "-" << j << "-" << k << " i-j-k in use\n";
//for (size_t i = 0; i < verts.size(); i++)verts[i]->color = color; //default mesh color for all vertices
		if (onePathOnly)
			return;
		sphereIdxs.push_back(-1); //finish red spheres; more sphere painting (green and white) by i to k and j to k paths below
		for (size_t z = 0; ! verts[ samples[i] ]->pathSamples.empty() && z < verts[ samples[i] ]->pathSamples[k].size(); z++) sphereIdxs.push_back(verts[ samples[i] ]->pathSamples[k][z]);sphereIdxs.push_back(-1); //i to k (i to j path already drawn above)
		for (size_t z = 0; ! verts[ samples[j] ]->pathSamples.empty() && z < verts[ samples[j] ]->pathSamples[k].size(); z++) sphereIdxs.push_back(verts[ samples[j] ]->pathSamples[k][z]);sphereIdxs.push_back(-1); //j to k (white)
		return;
	}
	else if (patchColoring == 2) //no descriptor coloring; just color area rings of based on v.binID values
	{		
		float ringColors[10][3] = {{1.0f, 0.68f, 0.79f}, {0.0f, 0.0f, 0.7f}, {0.0f, 0.7f, 0.0f}, {7.0f, 0.0f, 0.0f}, {0.7f, 0.0f, 0.7f}, {0.0f, 0.7f, 0.7f}, {0.7f, 0.7f, 0.0f}, {0.1f, 0.3f, 0.5f}, {0.5f, 0.3f, 0.9f}, {0.9f, 0.5f, 0.1f}}; //pink, blue, green, ..
		for (size_t i = 0; i < verts.size(); i++)
			if (verts[i]->binID == -1 || verts[i]->binID > 9)
				verts[i]->color = color; //default mesh color for the non-ring vertices or out-of-bound ring vertices (drawing such extra bins in random non-default colors might be better)
			else
				verts[i]->color = new float[3], verts[i]->color[0] = ringColors[ verts[i]->binID ][0], verts[i]->color[1] = ringColors[ verts[i]->binID ][1], verts[i]->color[2] = ringColors[ verts[i]->binID ][2];
		return;
	}
	else if (patchColoring == 3) //no descriptor coloring; just color rois based on v.roi values in yellow
	{		
		for (size_t i = 0; i < verts.size(); i++)
			if (! verts[i]->roi)
				verts[i]->color = color; //default mesh color for the non-roi vertices
			else
				verts[i]->color = new float[3], verts[i]->color[0] = verts[i]->color[1] = 0.7f, verts[i]->color[2] = 0.0f;
//for (size_t i = 0; i < verts.size(); i++)verts[i]->color = color; //default mesh color for all vertices
		return;
	}
	//descriptor coloring based on the descriptor selected here (tmp keeps the selected descriptor for the rest of this function)
	for (size_t i = 0; i < verts.size(); i++)
		verts[i]->tmp = verts[i]->geoToRef;//agd;//area,geoToRef,gpsScalar,eigfunc;
	double valMin = INF, valMax = -INF, dv, tiny2 = TINY*TINY, total = 0.0; //user could specify [valMin,valMax] as well in which case clipping below makes sense (now it's redundant)
	bool naiveColoring = false; //descriptor values must be in [0, 1] in naive mode (not recommended)
	for (size_t i = 0; i < verts.size(); i++)
	{
		if (verts[i]->tmp < valMin) valMin = verts[i]->tmp;
		if (verts[i]->tmp > valMax) valMax = verts[i]->tmp;
		total += verts[i]->tmp;
	}
	cout << "range for coloring: [" << valMin << "," << valMax << "], avg: " << total / verts.size() << endl;
	dv = (valMax - valMin > 0 ? valMax - valMin : tiny2);
	if (naiveColoring)
		for (size_t i = 0; i < verts.size(); i++)
			verts[i]->tmp = (verts[i]->tmp - valMin) / dv;
	//smooth coloring based on descriptor values
	for (size_t i = 0; i < verts.size(); i++)
	{
		verts[i]->color = new float[3];
		if (roiReady && ! verts[i]->roi) //if roi not available, i.e., NODENSEMAP mode, then no continue possibility as i.roi are not defined properly
		{
			for (int c = 0; c < 3; c++)
				verts[i]->color[c] = color[c];
			continue; //gray mesh color on this non-roi mesh2 vertex
		}
		//from descriptor=vmax red to descriptor=vmin blue: to make this happen do red (1) to green (0.5) to blue (0) where cyans are more visible
		double val = verts[i]->tmp;
		val = (val < valMin ? valMin : val); //clip val if it's outside [valMin,valMax]; currently redundant as valMin/Max are found automatically above
		val = (val > valMax ? valMax : val);		
		if (dv < tiny2) //very rare case where valMin=valMax (0 range for coloring; constant eigfunc from the 1st eigenfunction gives a range of 2.0435e-15 instead of 0 due to precision errors)
			verts[i]->color[0] = 0.0f, verts[i]->color[1] = verts[i]->color[2] = 1.0f;
		else if (val < (valMin + 0.25*dv))
			verts[i]->color[0] = 0.0f, verts[i]->color[1] = (float) (4 * (val - valMin) / dv), verts[i]->color[2] = 1.0f;
		else if (val < (valMin + 0.5*dv))
			verts[i]->color[0] = 0.0f, verts[i]->color[1] = 1.0f, verts[i]->color[2] = (float) (1 + 4 * (valMin + 0.25*dv - val) / dv);
		else if (val < (valMin + 0.75*dv))
			verts[i]->color[0] = (float) (4 * (val - valMin - 0.5*dv) / dv), verts[i]->color[1] = 1.0f, verts[i]->color[2] = 0.0f;
		else
			verts[i]->color[0] = 1.0f, verts[i]->color[1] = (float) (1 + 4 * (valMin + 0.75*dv - val) / dv), verts[i]->color[2] = 0.0f;

		//////////// alternative color map ////////////
		if (naiveColoring) //overwrites above by going from agd=1 red to agd=0 blue	(requires [0,1] interval for the descriptors; no such requirement above)
			for (int c = 0; c < 3; c++)
			{
				//overwrites above by going from agd=1 yellow to agd=0 blue
//				double blue[3] = {0.0, 0.0, 1.0}, yellow[3] = {1.0, 1.0, 0.0}; val = verts[i]->tmp; verts[i]->color[c] = (float) ((1.0 - val)*blue[c] + val*yellow[c]); continue;

				//from agd=1 red to agd=0 blue: to make this happen do red (1) to green (0.5) to blue (0) where cyans are less visible
				val = verts[i]->tmp*2.0;
				if (val < 1.0)
				{
					if (c == 0)
						verts[i]->color[c] = 0.0f; //no red
					else if (c == 1)
						verts[i]->color[c] = (float) val;
					else
						verts[i]->color[c] = 1.0f - verts[i]->color[1]; //1 - green
				}
				else
				{
					if (c == 0)
						verts[i]->color[c] = (float) val - 1.0f;
					else if (c == 1)
						verts[i]->color[c] = 1.0f - verts[i]->color[0]; //1 - red
					else
						verts[i]->color[c] = 0.0f; //no blue
				}
			}
		//////////// alternative color map ends ////////////
	}//*/
}
