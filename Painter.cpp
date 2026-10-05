#include "Painter.h"

SoSeparator* Painter::shapeSep(Mesh* mesh, double deltaX)
{
	//puts a shape (given as mesh.verts/tris) to the screen

	SoSeparator* sep = new SoSeparator();

	//transformation
	//not needed

	//material
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(mesh->color); //single coloring
	//multicoloring
	for (int i = 0; i < (int) mesh->verts.size(); i++)
		mat->diffuseColor.set1Value(i, (mesh->verts[i]->color ? mesh->verts[i]->color : mesh->color)); //multicoloring
		//mat->diffuseColor.set1Value(i, mesh->color); //single coloring
//	mat->transparency = 0.5f; //1: transparent, 0: opaque

	//a default creaseAngle of 0 (see SoShapeHints) does flat shading and of pi does Gouraud
	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = (float) PI;

	SoMaterialBinding* materialBinding = NULL;
//	if (multicoloring)
//	{
		materialBinding = new SoMaterialBinding; //for 2+ diffuse color usage on the same mesh
		materialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED;
//	}

	//shape
	SoCoordinate3* coords = new SoCoordinate3();
	for (int v = 0; v < (int) mesh->verts.size(); v++)
		coords->point.set1Value(v, (float) (mesh->verts[v]->coords[0] + deltaX), (float) mesh->verts[v]->coords[1], (float) (mesh->verts[v]->coords[2]));
	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	int nt = 0;
	for (int t = 0; t < (int) mesh->tris.size(); t++)
	{
		faceSet->coordIndex.set1Value(0 + 4*nt, mesh->tris[t]->v1i);
		faceSet->coordIndex.set1Value(1 + 4*nt, mesh->tris[t]->v2i);
		faceSet->coordIndex.set1Value(2 + 4*nt, mesh->tris[t]->v3i);
		faceSet->coordIndex.set1Value(3 + 4*nt, -1);
//		if (multicoloring)
//		{
			faceSet->materialIndex.set1Value(0 + 4*nt, mesh->tris[t]->v1i);
			faceSet->materialIndex.set1Value(1 + 4*nt, mesh->tris[t]->v2i);
			faceSet->materialIndex.set1Value(2 + 4*nt, mesh->tris[t]->v3i);
//		}
		nt++;
	}
	sep->addChild(mat);
	sep->addChild(hints); //Gouraud shading
//	if (multicoloring)
		sep->addChild(materialBinding);
	sep->addChild(coords);
	sep->addChild(faceSet);
	return sep;
}

int wbDrawing = 0; //worst or best match drawing
SoSeparator* Painter::sphereSep(Mesh* mesh, int v, double deltaX)
{
	//returns sphere on mesh.verts[v]

	SoSeparator* sphereSep = new SoSeparator;
	if (v < 0 || v > (int) mesh->verts.size())
		return sphereSep;
	double radius = 1.5 * mesh->avgEdgeLen * (mesh->verts.size() > 50000 ? 3 : 1); //1.95 * mesh->avgEdgeLen
radius *= 1.2; //slightly larger than spheresSep() spheres
//if(v==26554)radius*=1.6;
//if(v==662)radius*=2;
	//material
	SoMaterial* ma = new SoMaterial;
	ma->diffuseColor.setValue(SbColor(0.9f, 0.9f, (v == 34692 ? 0.0f : 0.9f)));
	if (wbDrawing == 1) //special flag for the worst match
		ma->diffuseColor.setValue(SbColor(0.7f, 0.0f, 0.0f)); //red
	else if (wbDrawing == 2) //special flag for the best match
		ma->diffuseColor.setValue(SbColor(0.0f, 0.7f, 0.0f)); //green
	else if (wbDrawing == 3) //special flag for the unmatched sample
		ma->diffuseColor.setValue(SbColor(0.0f, 0.0f, 0.7f)); //blue
	sphereSep->addChild(ma); //ok if i do this before SoTransform

	//transformation
	SoTransform* tra = new SoTransform();
	tra->translation.setValue((float) mesh->verts[v]->coords[0] + (float) deltaX, (float) mesh->verts[v]->coords[1], (float) mesh->verts[v]->coords[2]);
	sphereSep->addChild(tra);	
		
	//shape
	SoSphere* sph1 = new SoSphere();
	sph1->radius = (float) radius;
//sph1->radius = 0;//make mesh1 sphere invisible
	sphereSep->addChild(sph1); //whose position is decided by the translation applied above
	return sphereSep;
}

SoSeparator* Painter::spheresSep(Mesh* mesh, double deltaX)
{
	//put spheres to screen on (matched or all) mesh.samples locations

	SoSeparator* spheresSep = new SoSeparator;
//return spheresSep;
	double radius = (mesh->samples.size() >= 250 ? 0.5 : 1.5) * mesh->avgEdgeLen * (mesh->verts.size() > 50000 ? 3 : 1);
radius *= 1.3;
	for (size_t b = 0; b < mesh->samples.size(); b++)
	{
		if (mesh->verts[ mesh->samples[b] ]->matchIdx == -1)
			continue; //see only the matched samples

		//shape
		SoSeparator* sphere1Sep = new SoSeparator;

		//transformation
		SoTransform* tra = new SoTransform();
		tra->translation.setValue((float) mesh->verts[ mesh->samples[b] ]->coords[0] + (float) deltaX, (float) mesh->verts[ mesh->samples[b] ]->coords[1], (float) mesh->verts[ mesh->samples[b] ]->coords[2]);
		sphere1Sep->addChild(tra);

		//material
		SoMaterial* ma = new SoMaterial;
		//generate a random real number in [0,1] using (float) rand() / RAND_MAX
		SbColor color = SbColor((mesh->id == 6 ? 1.0f : 0.0f), 0.7f, 0.7f);//(float) rand() / RAND_MAX, (float) rand() / RAND_MAX, (float) rand() / RAND_MAX); //srand() already done in main()
		ma->diffuseColor.setValue(color);
		sphere1Sep->addChild(ma);

		//shape
		SoSphere* sph1 = new SoSphere();
		sph1->radius = (float) radius;
//sph1->radius = 0;//make mesh1 sphere invisible		
		sphere1Sep->addChild(sph1); //whose position is decided by the translation applied above
		spheresSep->addChild(sphere1Sep);
//cout<<mesh->samples[b]<<"spherized\n";
	}
	return spheresSep;
}

inline double distanceBetween(double* v1, double* v2)
{
	//returns distance b/w v1 and v2
	return sqrt((v2[0]-v1[0])*(v2[0]-v1[0]) +  (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]));
}

SoSeparator* Painter::spheresSep(Mesh* mesh1, Mesh* mesh2, double deltaX)
{
	//put spheres to screen on matched mesh1/2.samples locations w/ same colors for matching samples

	SoSeparator* spheresSep = new SoSeparator;
//	if (mesh1->samples.size() > 250)
//		return spheresSep; //too many spheres on screen may slow down rendering
	double radius1 = (mesh1->samples.size() > 300 ? 0.5 : 1.75) * mesh1->avgEdgeLen * (mesh1->verts.size() > 50000 ? 3 : 1),
		   radius2 = (mesh2->samples.size() > 300 ? 0.5 : 1.75) * mesh2->avgEdgeLen * (mesh2->verts.size() > 50000 ? 3 : 1);//tot = 0, n = 0;
//radius2 /= 10; //to see the manual spheres clearly
//radius2 *= 1.5;
	for (size_t b = 0; b < mesh1->samples.size(); b++)
	{
/*mesh1->verts[ mesh1->samples[b] ]->matchIdx = -1; //change them all
bool changed = (mesh1->verts[ mesh1->samples[b] ]->matchIdx == -1); int winner = -1; double minDist = INF;
mesh1->verts[ mesh1->samples[b] ]->matchIdx = (mesh1->verts[ mesh1->samples[b] ]->matchIdx != -1 ? mesh1->verts[ mesh1->samples[b] ]->matchIdx : mesh1->verts[ mesh1->samples[b] ]->gtIdx);
if (changed && mesh1->verts[ mesh1->samples[b] ]->matchIdx != -1) //alternative to assignment above where winner rhs is guaranteed to be a sample (above it can be a very close vertex)
{
	for (size_t j = 0; j < mesh2->samples.size(); j++)
		if (mesh2->verts[ mesh2->samples[j] ]->matchIdx == -1 && distanceBetween(mesh2->verts[ mesh2->samples[j] ]->coords, mesh2->verts[ mesh1->verts[ mesh1->samples[b] ]->matchIdx ]->coords) < minDist)
			minDist = distanceBetween(mesh2->verts[ mesh2->samples[j] ]->coords, mesh2->verts[ mesh1->verts[ mesh1->samples[b] ]->matchIdx ]->coords), winner = mesh2->samples[j];
	//if (winner != -1)
	if (winner != -1 && distanceBetween(mesh2->verts[ winner ]->coords, mesh2->verts[ mesh1->verts[ mesh1->samples[b] ]->gtIdx ]->coords) < 10)
		mesh1->verts[ mesh1->samples[b] ]->matchIdx = winner, mesh2->verts[winner]->matchIdx = mesh1->samples[b];
}//*/
//cout << mesh1->verts[ mesh1->samples[b] ]->gtDistortion << "\n"; tot += (mesh1->verts[ mesh1->samples[b] ]->gtDistortion > 0 ? mesh1->verts[ mesh1->samples[b] ]->gtDistortion : 0), n += (mesh1->verts[ mesh1->samples[b] ]->gtDistortion > 0);
//if (mesh1->verts[ mesh1->samples[b] ]->gtDistortion > 0.15 && mesh1->verts[ mesh1->samples[b] ]->gtIdx != -1) mesh1->verts[ mesh1->samples[b] ]->matchIdx = mesh1->verts[ mesh1->samples[b] ]->gtIdx, cout << "gtUpdate for " << mesh1->samples[b] << "\n\n";
//if ((float) rand() / RAND_MAX < 0.3) mesh1->verts[ mesh1->samples[b] ]->matchIdx = mesh1->verts[ mesh1->samples[b] ]->gtIdx;

//		if (mesh1->verts[ mesh1->samples[b] ]->matchIdx == -1)
//			continue; //see only the matched samples
		int v2 = mesh1->verts[ mesh1->samples[b] ]->matchIdx;

		//shape
		SoSeparator* sphere1Sep = new SoSeparator, * sphere2Sep = new SoSeparator;

		//transformation
		SoTransform* tra = new SoTransform();
		tra->translation.setValue((float) mesh1->verts[ mesh1->samples[b] ]->coords[0], (float) mesh1->verts[ mesh1->samples[b] ]->coords[1], (float) mesh1->verts[ mesh1->samples[b] ]->coords[2]);
		sphere1Sep->addChild(tra);
		SoTransform* tra2 = new SoTransform();
		if (v2 != -1) //prevent crash in case i disable continue above
			tra2->translation.setValue((float) mesh2->verts[v2]->coords[0] + (float) deltaX, (float) mesh2->verts[v2]->coords[1], (float) mesh2->verts[v2]->coords[2]);
		sphere2Sep->addChild(tra2);

		//material
		SoMaterial* ma = new SoMaterial;
		//generate a random real number in [0,1] using (float) rand() / RAND_MAX
		SbColor color = SbColor((float) rand() / RAND_MAX, (float) rand() / RAND_MAX, (float) rand() / RAND_MAX); //srand() already done in main()
		ma->diffuseColor.setValue(color);
		sphere1Sep->addChild(ma);
		SoMaterial* ma2 = new SoMaterial;
		ma2->diffuseColor.setValue(color);
		sphere2Sep->addChild(ma2);

//ma->diffuseColor.setValue( SbColor(1.0f, 0.7f, 0.7f) );ma2->diffuseColor.setValue( SbColor(0.0f, 0.7f, 0.7f) ); //pink vs. cyan coloring on mesh1 vs. mesh2

		//shape
		SoSphere* sph1 = new SoSphere(), * sph2 = new SoSphere();
		sph1->radius = (float) radius1, sph2->radius = (float) radius2;
sph1->radius = (float) radius2;
//sph1->radius1 = 0;//make mesh1 sphere invisible		
		sphere1Sep->addChild(sph1), sphere2Sep->addChild(sph2); //whose position is decided by the translation applied above
		spheresSep->addChild(sphere1Sep);
		if (v2 != -1) //no tra2 set above so prevent addChild() from producign weird spheres
			spheresSep->addChild(sphere2Sep);
//cout<<b<<" "<<mesh1->samples[b]<<"-"<<v2<<"spherized\n";
	}
//cout << n << "\tavg gt distortion: " << tot / n << "\n";
	return spheresSep;
}

SbColor color1 = SbColor(1, 0, 0), color2 = SbColor(0, 1, 0), color3 = SbColor(1, 1, 1);//(float) rand() / RAND_MAX, (float) rand() / RAND_MAX, (float) rand() / RAND_MAX);
SoSeparator* Painter::spheresSep2(Mesh* mesh, double deltaX)
{
	//put spheres to screen on the mesh.sphereIdxs locations (for path samples) separated by special -1 values

	SoSeparator* spheresSep = new SoSeparator;
//return spheresSep;
	double radius = 1.5 * mesh->avgEdgeLen * (mesh->verts.size() > 50000 ? 3 : 1);
	size_t b, c = 0;
	for (b = 0; b < mesh->sphereIdxs.size() && mesh->sphereIdxs[b] != -1; b++)
	{
		//cout << b << " " << mesh->sphereIdxs[b] << " sphere on red path\n";
		//shape
		SoSeparator* sphere1Sep = new SoSeparator;

		//transformation
		SoTransform* tra = new SoTransform();
		tra->translation.setValue((float) mesh->verts[ mesh->sphereIdxs[b] ]->coords[0] + (float) deltaX, (float) mesh->verts[ mesh->sphereIdxs[b] ]->coords[1], (float) mesh->verts[ mesh->sphereIdxs[b] ]->coords[2]);
		sphere1Sep->addChild(tra);

		//material
		SoMaterial* ma = new SoMaterial;
		ma->diffuseColor.setValue(color1);
		sphere1Sep->addChild(ma);

		//shape
		SoSphere* sph1 = new SoSphere();
		sph1->radius = (float) radius;
		sphere1Sep->addChild(sph1); //whose position is decided by the translation applied above
		spheresSep->addChild(sphere1Sep);
		c++;
	}
//cout << c << " red, ", c = 0;
	//second path's spheres
	for (b++; b < mesh->sphereIdxs.size() && mesh->sphereIdxs[b] != -1; b++)
	{
		//cout << b << " " << mesh->sphereIdxs[b] << " sphere on green path\n";
		//shape
		SoSeparator* sphere1Sep = new SoSeparator;

		//transformation
		SoTransform* tra = new SoTransform();
		tra->translation.setValue((float) mesh->verts[ mesh->sphereIdxs[b] ]->coords[0] + (float) deltaX, (float) mesh->verts[ mesh->sphereIdxs[b] ]->coords[1], (float) mesh->verts[ mesh->sphereIdxs[b] ]->coords[2]);
		sphere1Sep->addChild(tra);

		//material
		SoMaterial* ma = new SoMaterial;
		ma->diffuseColor.setValue(color2);
		sphere1Sep->addChild(ma);

		//shape
		SoSphere* sph1 = new SoSphere();
		sph1->radius = (float) radius;
		sphere1Sep->addChild(sph1); //whose position is decided by the translation applied above
		spheresSep->addChild(sphere1Sep);
		c++;
	}
//cout << c << " green, ", c = 0;
	//third path's spheres
	for (b++; b < mesh->sphereIdxs.size() && mesh->sphereIdxs[b] != -1; b++)
	{
		//cout << b << " " << mesh->sphereIdxs[b] << " sphere on white path\n";
		//shape
		SoSeparator* sphere1Sep = new SoSeparator;

		//transformation
		SoTransform* tra = new SoTransform();
		tra->translation.setValue((float) mesh->verts[ mesh->sphereIdxs[b] ]->coords[0] + (float) deltaX, (float) mesh->verts[ mesh->sphereIdxs[b] ]->coords[1], (float) mesh->verts[ mesh->sphereIdxs[b] ]->coords[2]);
		sphere1Sep->addChild(tra);

		//material
		SoMaterial* ma = new SoMaterial;
		ma->diffuseColor.setValue(color3);
		sphere1Sep->addChild(ma);

		//shape
		SoSphere* sph1 = new SoSphere();
		sph1->radius = (float) radius;
		sphere1Sep->addChild(sph1); //whose position is decided by the translation applied above
		spheresSep->addChild(sphere1Sep);
		c++;
	}
//cout << c << " white on screen\n";
	return spheresSep;
}

SoSeparator* Painter::spheresSep3(Mesh* mesh, bool pathVerts, double deltaX)
{
	//put spheres to screen on the v.pathVert/locMaxGC=true vertex locations

	SoSeparator* spheresSep = new SoSeparator;
//return spheresSep;
	double radius = (pathVerts ? 0.2 : 1.2) * mesh->avgEdgeLen * (mesh->verts.size() > 50000 ? 3 : 1);

//	double radius = (mesh->samples.size() >= 250 ? 0.5 : 1.5) * mesh->avgEdgeLen * (mesh->verts.size() > 50000 ? 3 : 1);
//radius *= 1.3;

	int counts[9] = {0};
	for (size_t v = 0; v < mesh->verts.size(); v++)
		if ((pathVerts && mesh->verts[v]->pathVert) || //i may validly see some gaps on path vertex drawing due to support
			(! pathVerts && mesh->verts[v]->locMaxGC)) //a special vertex in red, e.g., localMaximaOfGC
			//(! pathVerts && mesh->verts[v]->locMinNGC)) //a special vertex in red, e.g., localMinimaOfNegativeGC
			//mesh->verts[v]->sample) //only samples will be drawn with colors given by tmp below
		{
			//shape
			SoSeparator* sphere1Sep = new SoSeparator;

			//transformation
			SoTransform* tra = new SoTransform();
			tra->translation.setValue((float) mesh->verts[v]->coords[0] + (float) deltaX, (float) mesh->verts[v]->coords[1], (float) mesh->verts[v]->coords[2]);
			sphere1Sep->addChild(tra);

			//material
			SoMaterial* ma = new SoMaterial;
			//generate a random real number in [0,1] using (float) rand() / RAND_MAX
			SbColor color = SbColor((mesh->id == 0 ? 1.0f : 0.0f), 0.7f, 0.7f);//(float) rand() / RAND_MAX, (float) rand() / RAND_MAX, (float) rand() / RAND_MAX); //srand() already done in main()
			ma->diffuseColor.setValue(color);
			sphere1Sep->addChild(ma);
			if (mesh->verts[v]->tmp == 0)
				ma->diffuseColor.setValue( SbColor(1.0, 0.0, 0.0) ), counts[ (int) mesh->verts[v]->tmp ]++; //30
			else if (mesh->verts[v]->tmp == 1)
				ma->diffuseColor.setValue( SbColor(0.0, 1.0, 0.0) ), counts[ (int) mesh->verts[v]->tmp ]++;//60
			else if (mesh->verts[v]->tmp == 2)
				ma->diffuseColor.setValue( SbColor(0.0, 0.0, 1.0) ), counts[ (int) mesh->verts[v]->tmp ]++;//90
			else if (mesh->verts[v]->tmp == 3)
				ma->diffuseColor.setValue( SbColor(1.0, 1.0, 0.0) ), counts[ (int) mesh->verts[v]->tmp ]++;//120
			else if (mesh->verts[v]->tmp == 4)
				ma->diffuseColor.setValue( SbColor(1.0, 0.0, 1.0) ), counts[ (int) mesh->verts[v]->tmp ]++;//150
			else if (mesh->verts[v]->tmp == 5)
				ma->diffuseColor.setValue( SbColor(0.0, 1.0, 1.0) ), counts[ (int) mesh->verts[v]->tmp ]++;//180
			else if (mesh->verts[v]->tmp == 6)
				ma->diffuseColor.setValue( SbColor(1.0, 1.0, 1.0) ), counts[ (int) mesh->verts[v]->tmp ]++;//210
			else if (mesh->verts[v]->tmp == 7)
				ma->diffuseColor.setValue( SbColor(0.7f, 0.7f, 0.7f) ), counts[ (int) mesh->verts[v]->tmp ]++;//240
			else if (mesh->verts[v]->tmp == 8)
				ma->diffuseColor.setValue( SbColor(0.0, 0.0, 0.0) ), counts[ (int) mesh->verts[v]->tmp ]++;//250
//color = SbColor((mesh->id == 6 ? 1.0f : 0.0f), 0.7f, 0.7f), ma->diffuseColor.setValue(color);
			
			//shape
			SoSphere* sph1 = new SoSphere();
			sph1->radius = (float) radius;
//sph1->radius = 0;//make mesh1 sphere invisible		
			sphere1Sep->addChild(sph1); //whose position is decided by the translation applied above
			spheresSep->addChild(sphere1Sep);
//cout<<v<<"spherized\n";
		}
//for(int i=0;i<9;i++)cout << counts[i] << " for color " << i << "\n";
	return spheresSep;
}

SoSeparator* Painter::spheresSep4(Mesh* mesh1, Mesh* mesh2, int worst, int best, double deltaX)
{
	//put spheres to screen on the endpoints of the worst/best matches (and optionally the unmatched samples)

	SoSeparator* spheresSep = new SoSeparator;
	if (best < 0 || best > (int) mesh1->verts.size() || worst < 0 || worst > (int) mesh1->verts.size())
		return spheresSep;
	
	wbDrawing = 1, spheresSep->addChild( sphereSep(mesh1, worst) ), spheresSep->addChild( sphereSep(mesh2, mesh1->verts[worst]->matchIdx, deltaX) );//cout << "worst: " << worst << "-" << mesh1->verts[worst]->matchIdx << " " << mesh1->verts[worst]->distortion << endl;
	wbDrawing = 2, spheresSep->addChild( sphereSep(mesh1, best) ), spheresSep->addChild( sphereSep(mesh2, mesh1->verts[best]->matchIdx, deltaX) );//cout << "best: " << best << "-" << mesh1->verts[best]->matchIdx << " " << mesh1->verts[best]->distortion << endl;
	wbDrawing = 3; //display unmatcheds in blue (optional)
	for (size_t i = 0; i < mesh1->samples.size(); i++)
		if (mesh1->verts[ mesh1->samples[i] ]->matchIdx == -1)
			spheresSep->addChild( sphereSep(mesh1, mesh1->samples[i]) ), cout << mesh1->samples[i] << " "; cout << "unmatched1\n";
	for (size_t i = 0; i < mesh2->samples.size(); i++)
		if (mesh2->verts[ mesh2->samples[i] ]->matchIdx == -1)
			spheresSep->addChild( sphereSep(mesh2, mesh2->samples[i], deltaX) ), cout << mesh2->samples[i] << " "; cout << "unmatched2\n";//*/
	wbDrawing = 0; //for potential future calls to sphereSep()
	/*double radius1 = 1.95 * mesh1->avgEdgeLen * (mesh1->verts.size() > 50000 ? 3 : 1), radius2 = 1.95 * mesh2->avgEdgeLen * (mesh2->verts.size() > 50000 ? 3 : 1);
	for (size_t v = 0; v < mesh1->verts.size(); v++) marked-based coloring, e.g., show coarseMap which has v.marked=true
		if (mesh1->verts[v]->marked && mesh1->verts[v]->matchIdx != -1) //marked vertices have a bad match for sure (so matchIdx!=-1 redundant but anyway)
		{
			//shape
			SoSeparator* sphere1Sep = new SoSeparator;
			//transformation
			SoTransform* tra = new SoTransform();
			tra->translation.setValue((float) mesh1->verts[v]->coords[0], (float) mesh1->verts[v]->coords[1], (float) mesh1->verts[v]->coords[2]);
			sphere1Sep->addChild(tra);
			//material (to be shared by v.matchIdx below)
			SoMaterial* ma = new SoMaterial;
			//generate a random real number in [0,1] using (float) rand() / RAND_MAX
			SbColor color = SbColor((float) rand() / RAND_MAX, (float) rand() / RAND_MAX, (float) rand() / RAND_MAX); //srand() already done in main()
			ma->diffuseColor.setValue(color);
			sphere1Sep->addChild(ma);
			//shape
			SoSphere* sph1 = new SoSphere();
			sph1->radius = (float) radius1;
			sphere1Sep->addChild(sph1); //whose position is decided by the translation applied above
			spheresSep->addChild(sphere1Sep);

			//shape
			SoSeparator* sphere2Sep = new SoSeparator;
			SoTransform* tra2 = new SoTransform();
			tra2->translation.setValue((float) mesh2->verts[ mesh1->verts[v]->matchIdx ]->coords[0] + (float) deltaX, (float) mesh2->verts[ mesh1->verts[v]->matchIdx ]->coords[1], (float) mesh2->verts[ mesh1->verts[v]->matchIdx ]->coords[2]);
			sphere2Sep->addChild(tra2);
			sphere2Sep->addChild(ma);
			SoSphere* sph2 = new SoSphere();
			sph2->radius = (float) radius2;
			sphere2Sep->addChild(sph2); //whose position is decided by the translation applied above
			spheresSep->addChild(sphere2Sep);
		}//*/
	return spheresSep;
}

SoSeparator* Painter::matchingLinesSep(Mesh* mesh1, Mesh* mesh2, double deltaX)
{
	//returns the lines going from mesh1.v to mesh1.v.matchIdx in mesh2 (mesh1 is source, mesh2 is target)
	
	SoSeparator* linesSep = new SoSeparator();

	SoMaterial* ma = new SoMaterial;
	for (int i = 0; i < (int) mesh1->verts.size(); i++)
		ma->diffuseColor.set1Value(i, 0.0f, 0.0f, 0.0f); //black lines
	linesSep->addChild(ma);
//SoDrawStyle* style = new SoDrawStyle; style->lineWidth = 1.5f; linesSep->addChild(style); //thick lines
	
	SoIndexedLineSet* lines = new SoIndexedLineSet;
	SoCoordinate3* co = new SoCoordinate3;	
	int lc = 0, everyOtherLine = 1; //make it 1 to draw all lines and x>1 to draw every other xth line
	for (int b = 0; b < (int) mesh1->verts.size(); b += everyOtherLine)
	{
		int corresp = mesh1->verts[b]->matchIdx;
		if (corresp == -1 || ! mesh1->verts[b]->sample)
			continue; //no match for this mesh1 vertex in mesh2
		SbVec3f from((float) mesh1->verts[b]->coords[0], (float) mesh1->verts[b]->coords[1], (float) mesh1->verts[b]->coords[2]),
				to((float) (mesh2->verts[corresp]->coords[0] + deltaX), (float) mesh2->verts[corresp]->coords[1], (float) (mesh2->verts[corresp]->coords[2]));
		co->point.set1Value(2*lc, from);
		co->point.set1Value(2*lc + 1, to);
		lines->coordIndex.set1Value(3*lc, 2*lc);
		lines->coordIndex.set1Value(3*lc + 1, 2*lc + 1);
		lines->coordIndex.set1Value(3*lc + 2, -1); //end this line with -1
		lc++;
//cout<<b<<"-"<<corresp<<"lined\n";
	}
cout << lc << " lines on screen" << (everyOtherLine > 1 ? " (some matches skipped)" : "\n");
	linesSep->addChild(co);
	linesSep->addChild(lines);

	return linesSep;
}
