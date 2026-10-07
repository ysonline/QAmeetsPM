//Author: ANONYMIZED
//Paper: ANONYMIZED (submitted)

#define _CRT_SECURE_NO_WARNINGS //i won't use fopen_s, sprintf_s, etc; so, don't warn me

#include "Mesh.h"
#include "Correspondence.h"

#ifdef GRAPHICS_STUFF
#define _SECURE_SCL 0 //applications build on top of Open Inventor in release mode MUST define the _SECURE_SCL=0 preprocessor variable in order to avoid alignment problems. Using this 
					  //mode reduce greatly memory footprint of STL objects and increase performance by a factor 2 to 10 depending on operations
//#define SOWIN_DLL #define COIN_DLL //doing here not sufficient 'cos there're inventor include's in other files as well; so put these preprocessor definitions to Projects->C/C++->Preprocessor
#include "Painter.h"

//viewer
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#endif
double setAreaWeight(Mesh* mesh1, Mesh* mesh2, double& val);
int main(int argc, char ** argv)
{
	srand( (unsigned) time(NULL) ); //without srand the program would have generated the same number each time we ran it because the generator would have been seeded with the same default value each time; now it'll be seeded w/ sys time and generate different numbers
	int idSrc = 100, idTgt = 101, m = 5/*3 5 10*/, n = 10, //match n points on source to m points on target where m <= n
		nDenseSamples = 250/*30 100 250 350*/, qapEnergy = POINTWISEPAIRWISE, qapSolver = QUCOOPMATCH, nSamplesPerSubregion = 10/*5 10 15*/, nMaxC2FLevels = 10, //more parameters
		dissimType = POINTWISEPAIRWISE, //PAIRWISENOPATCH, PAIRWISE, POINTWISE, //terms of dissimilarity function that guides the QAP and mapQualityPatch/Virtual()
		denseMapType = (nDenseSamples == -1 ? NODENSEMAP : C2FQAP),//QAP,LAP,NODENSEMAP //lap or qap or coarse-to-fine qap, use NODENSEMAP for no dense matching, i.e., only coarse map is computed
		nPathSamples = 9/*9 49*/, //there'll be n+1 (including path endpoints) samples on each path (small paths may remain under n+1)
		borderTreat = 1,//0; //treat border samples by 1: doubling the areas around them, 2: excluding them during map search, 0: doing nothing special about borders
		nonuniformPathSampling = 2; //each path is sampled uniformly (0) or nonuniformly (1 or 2)
	bool curvFPS = true, eucPatching = false, //pathing always comes with patching which is based on fast eucldian distance or accurate geodesic distance (as patrches are supposed to be local euclidean is a good approximation)
		 //other params nRingAreas & radius set below as mesh globals (nRingAreas: # ring areas, radius: patch radius or stepSize during ring computation)
		 //more ablations: overlapping areas[] vs. non-overlapping areasRing[] for the unary term and also patchCostDistortion = true vs. false during the C2F process
		 loadP3D = false, //load the results created by the program of the paper: Partial 3-D Correspondence from Shape Extremities
		 loadDL = false, //load the results created by the program of 3 deep learning papers: HFM, ULRSS, EchoMatch
		 loadQM = false, //qmatch-style algorithm that removes k worst matches during its iterative updates and applies qmatch/qucoop/qucm solver along the way		 
		 fileFPSforQM = false, //use qmatch-style algorithm's samples from file for a fair comparison
		 fastC2F = false, //fast coarse-to-fine uses less nMaxC2FLevels and nSamplesPerSubregion (nPathSamples is always 5 in dense C2F mode as set by denseMap())
		 fmRefinement = (denseMapType != NODENSEMAP),//false //refine c2f matches via fullMap() that is based on c2f matches
		 displayFullMap = false; //show the interpolated dense, e.g., denseMapType=C2FQAP, or coarse, i.e., denseMapType=NODENSEMAP, map which is a full map between all vertices	
	double maxAreaRatio = 0.5/*0.45 0.5*/; //partial (or complete) mesh2 area / complete mesh1 area
	if (fastC2F)
		nMaxC2FLevels = nSamplesPerSubregion = 5;		
	if (dissimType == POINTWISE || dissimType == PAIRWISENOPATCH)
		nPathSamples = -1; //no pathing and patching so don't waste time on pathSampling()
	cout << "<<<<<<<<<< Quantum Annealing Meets Partial Matching >>>>>>>>>>\n";
	cout << "n = m recommended for complete-to-complete matching and m <= n for complete-to-partial matching where\nn: # of initial samples on the complete source model and\nm: # of initial samples on the partial (can also be complete) target model"
		 << "\nthere is also nDenseSamples: # of dense ROI samples to be matched in coarse-to-fine manner\n\n";
	///////////////// params from command prompt /////////////////
#ifndef GRAPHICS_STUFF //public mode does not support my openinventor grapics unit for max portability so do below if GRAPHICS_STUFF is not defined		
	if (argc != 12)
	{
		cout << "11 parameters expected (not provided): PartialCorrespQAP.exe <sourceID> <targetID> <nInitialSourceSamples> <mInitialTargetSamples> <nDenseSamples> <qapEnergy> <nPathSamples> <qapSolver> <nSamplesPerSubregion> <nMaxC2FLevels> <fullMap>\n\n";
		cout << "<*ID>: nonnegative integer (meshes must be in off format named as <*ID>.off, default sourceID: 10 and targetID: 20)\n";
		cout << "<nInitialSourceSamples>: # of initial samples on the complete source model: n (paper default: 10)\n";
		cout << "<mInitialSourceSamples>: # of initial samples on the partial (can also be complete) target model: m where m <= n (paper default: 5)\n";
		cout << "<nDenseSamples>: # of dense ROI samples to be matched in coarse-to-fine manner; set to -1 to disable this stage (paper default: 30)\n";
		cout << "<qapEnergy>: QAP energy: 1: paper's Eq. 7 with vectoral terms (3D Q matrix created), 2: lighter scalar terms (2D Q matrix created) (paper default: 1)\n";
		cout << "<nPathSamples>: QAP energy = 1 samples this many vertices (plus 1) along paths (paper default: 9)\n";
		cout << "<qapSolver>: drop-in QAP solver: 1: QuCM (ours), 2: QuCOOP, 3: QMATCH, 4: QGM\n";
		cout << "<nSamplesPerSubregion>: # of samples inside each subregion (of ROI) during the coarse-to-fine process (paper default: 5)\n";
		cout << "<nMaxC2FLevels>: # of coarse-to-fine levels (paper default: 10)\n";
		cout << "<fullMap>: 1: convert nDenseSamples-size map to full map, 0: no conversion (paper default: 1)\n";
		cout << "defaults are in use\n\n";
	}
	else
		cout << "11 parameters expected and provided: PartialCorrespQAP.exe <sourceID> <targetID> <nInitialSourceSamples> <mInitialTargetSamples> <nDenseSamples> <qapEnergy> <nPathSamples> <qapSolver> <nSamplesPerSubregion> <nMaxC2FLevels> <fullMap>\n\n",
		idSrc = atoi(argv[1]), idTgt = atoi(argv[2]), n = atoi(argv[3]), m = atoi(argv[4]), nDenseSamples = atoi(argv[5]),
		qapEnergy = atoi(argv[6]), nPathSamples = atoi(argv[7]), qapSolver = atoi(argv[8]), nSamplesPerSubregion = atoi(argv[9]), nMaxC2FLevels = atoi(argv[10]), fmRefinement = (atoi(argv[11]) == 1);
	dissimType = qapEnergy, nPathSamples = (dissimType == POINTWISE || dissimType == PAIRWISENOPATCH ? -1 : nPathSamples), denseMapType = (nDenseSamples == -1 ? NODENSEMAP : C2FQAP);
#endif
	///////////////// params from command prompt /////////////////

	Mesh* mesh1 = new Mesh(idSrc, borderTreat, false, nonuniformPathSampling), * mesh2 = new Mesh(idTgt, borderTreat, true, nonuniformPathSampling); //mesh1 complete, mesh2 complete or partial model (so always make mesh1 the complete one)
	//source
	char fName1[250], fName2[250];
	sprintf(fName1, "%d.off", idSrc), sprintf(fName2, "%d.off", idTgt); 	
	if (loadDL) //load competitors' models and results
		//sprintf(fName1, "hybridfmaps-results-v2\\outputs\\shrec16\\raw-outputs\\scape200_decimated_david\\visualization\\%d.off", idSrc), mesh1->loadOff(fName1), sprintf(fName2, "hybridfmaps-results-v2\\outputs\\shrec16\\raw-outputs\\scape200_decimated_david\\visualization\\%d.off", idTgt), mesh2->loadOff(fName2); //HFM partial david run
		//mesh1->id = 120, mesh2->id = idTgt = 119, sprintf(fName1, "hybridfmaps-results-v2\\outputs\\shrec16\\raw-outputs\\scape200_decimated_david\\visualization\\%d.off", idSrc), mesh1->loadOff(fName1), sprintf(fName2, "hybridfmaps-results-v2\\outputs\\shrec16\\raw-outputs\\scape200_decimated_david\\visualization\\%d.off", idTgt), mesh2->loadOff(fName2); //HFM complete david run
		//sprintf(fName1, "toscahorse_quadric_05.off"), mesh1->loadOff(fName1), sprintf(fName2, "toscahorse_quadric_02.off"), mesh2->loadOff(fName2); //HFM horse run
		//sprintf(fName1, "ulrssm-shrec16\\shrec16_null_to_partial\\visualization\\%d.off", idSrc), mesh1->loadOff(fName1), sprintf(fName2, "ulrssm-shrec16\\shrec16_null_to_partial\\visualization\\%d.off", idTgt), mesh2->loadOff(fName2); //ULRSS david run
		//sprintf(fName1, "ulrssm-shrec16\\shrec16_null_to_partial\\visualization\\horse%d.off", idSrc), mesh1->loadOff(fName1), idTgt = mesh2->id = 6, sprintf(fName2, "ulrssm-shrec16\\shrec16_null_to_partial\\visualization\\horse%d.off", idTgt), mesh2->loadOff(fName2); //ULRSS horse run
		sprintf(fName1, "EchoMatchResults\\shrec16\\echo_match_shrec16_benchmark_dino50_\\visualization\\david%d.off", idSrc), mesh1->loadOff(fName1), sprintf(fName2, "EchoMatchResults\\shrec16\\echo_match_shrec16_benchmark_dino50_\\visualization\\david%d_partial3.off", idTgt), mesh2->loadOff(fName2); //EchoMatch partial david run
	else
		mesh1->loadOff(fName1, eucPatching), mesh2->loadOff(fName2, eucPatching);
//	if (! loadDL)
//		mesh2->loadGT(fName1, mesh1); //evaluation inactive in this public release

	Correspondence* corresp = NULL;
if (loadP3D) mesh1->resultFromFile(mesh2); //bad indentation on purpose as resultFromFile() is not a critical function for this project
else if (loadDL) mesh1->resultFromFile3(mesh2, 3); //bad indentation on purpose as resultFromFile() is not a critical function for this project
else if (loadQM) mesh1->id = 0, mesh2->id = 2, mesh1->FPSfile(250/*50*/, true), mesh2->FPSfile(250/*50*/, true), mesh1->resultFromFile2(mesh2);
else if (fileFPSforQM) mesh1->id = 0, mesh2->id = 2, corresp = new Correspondence(mesh1, mesh2), corresp->denseMap(denseMapType, true), corresp->mapQualityGT(); //do nothing (no qapMap() for coarseMap matching below) as fair comparison mode doesn't need it
else //common block for this paper's execution
	{
	double r1 = max(mesh1->avgEdgeLen, mesh2->avgEdgeLen) * 4.0, r2 = max(mesh1->avgEdgeLen, mesh2->avgEdgeLen) * 6.0 * (idTgt % 10 != 0/*partialMatching*/ ? 1 : 2.5), areaRatio = mesh2->totalArea / mesh1->totalArea;
	if (areaRatio < maxAreaRatio) //davidcomplete-partial1/2/3/4/5/6 = .52/.55/.45/.22/.147/.153 (<.16 means m=3 in partial model is in a very small region; to make the corresponding region on complete model sampled, use n=30 which makes bruteForce impossible but qap possible)
		r1 /= (nPathSamples > 20 ? 3 : 1), /*m = 3,*/mesh1->radius = mesh2->radius = r1;//n = (areaRatio < .16 ? 30 : n); //mesh2 is too small and big radius may make big overlaps/small distinctiveness
	else
		r2 /= (nPathSamples > 20 ? 3 : 1), mesh1->radius = mesh2->radius = r2; //patch radius (also works as the stepSize during areas())
	mesh1->nRingAreas = mesh2->nRingAreas = 5; //nRingAreas many ring areas (or overlapping accumulated areas) around each sample as its descriptor (local coverage (nRingAreas=5) recommended for partial matching 'cos as i go global i may hit missing parts that ruins consistency)

//	if (idSrc == 10 && idTgt == 13) //david0-david10partial3
//		m = 3, nonuniformPathSampling = mesh1->nonuniformPathSamplingHeuristic = mesh2->nonuniformPathSamplingHeuristic = 1, cout<<"\nm=3 for pfaust or m=3&pathSampler=1 davidpartial3\n\n";
	if(idTgt % 10 == 0) //completeMatching, e.g., david0-david1 or david0-david3
	{
		n = m = 8;//10; //completely similar so n can be equal to m and they can also be large
//		if (idSrc == 10 && idTgt == 20)
//			cout<<"\nn=m=8(10 makes incompatible sampling)&nPathSamples=19 for complete10-20\n\n", n = m = 8, nPathSamples = 19;
	}
	if (m > n || n < 2 || m < 2)
		cout << "WARNING: use m <= n where n and m >= 2\n\n"; //m=1 cannot define pairwise terms and m<1 is just meaningless

	//spectral dense descriptors
	//mesh1->lap->hksAndGPS(), mesh2->lap->hksAndGPS(); cout << endl; //usage=1, i.e., orthonormalWrtAreaMtrx=false	
	//spatial dense descriptors
	mesh1->curvaturesAll(), mesh2->curvaturesAll();
//	if (max(mesh1->stdDevGC / mesh2->stdDevGC, mesh2->stdDevGC / mesh1->stdDevGC) > 100 || max(mesh1->avgGC / mesh2->avgGC, mesh2->avgGC / mesh1->avgGC) > 10)
//		cout << "WARNING: gaussian curvatures incompatible; considering smoothing the bigGC-mesh so that its gc values get smaller\n\n"; //o/w fingertip.gc=0.7 vs. fingertip.gc=5799 prevents good fingertip matches

	//sampling (geodesics too)
	mesh1->FPS(n, curvFPS), mesh2->FPS(m, curvFPS); //curvFPS=true makes our modified FPS where curvature info is injected to the process
	if (nPathSamples != -1) //patches are centered at path samples
		mesh1->pathSampling(nPathSamples), mesh2->pathSampling(nPathSamples); //paths b/w samples are sampled	
	//spatial sparse descriptors (defined for samples only)
	if (dissimType == POINTWISE || dissimType == POINTWISEPAIRWISE) //pointwise descriptor is set by areas() so dissimType must include pointwise (to save tiny time)
	{
		mesh1->areas(), mesh2->areas();
		if (borderTreat == 1)
			mesh1->doubleBorderAreas(), mesh2->doubleBorderAreas(); //avgAreas() called inside doubleBorderAreas()
		else
			mesh1->avgAreas(), mesh2->avgAreas(); //set avgPatch/RingArea values used by areaWeight
	}

	//avgRingArea tends to be bigger than avgPatchArea, e.g. 700 vs. 300, so make them similar by areaWeight normalization
	double avgPairwise = 0, maxGeoDiff = max(mesh1->avgEdgeLen, mesh2->avgEdgeLen)*(idTgt % 10 != 0/*partialMatching*/?20:60); //20/60 is not about partiality: in tosca partial models are very low-reso so mesh2.avgEdgeLen is too high w.r.t. mesh1.avgEdgeLen (so no 20/60-update in denseMap() where partialMatching problem reduced to completeMatching)
//maxGeoDiff = INF; //disable geoDiff > maxGeoDiff thresholding trick in bruteForce and qapMap optimizations (maxGeoDiff variable in use only if GEODIFFTHRESHOLDING is defined)
	double d = max(mesh1->avgRingArea, mesh2->avgRingArea), awObserve = (d > 0 ? max(mesh1->avgPatchArea, mesh2->avgPatchArea) / d : 1.0), //weighting by observation
		   areaWeight = (dissimType == POINTWISEPAIRWISE ? setAreaWeight(mesh1, mesh2, avgPairwise) : awObserve); //automatic weighting that makes average unary and pairwise terms in Q (and in mapQualityPatch/Virtual()) approximately the same
//areaWeight = 1.0; //no unary weight, i.e., pairwise and unary features are uniformly weighted with 1.0
	corresp = new Correspondence(mesh1, mesh2, dissimType, maxGeoDiff, areaWeight, avgPairwise, nDenseSamples, nSamplesPerSubregion, nMaxC2FLevels); //unitGeoScaling() called inside
#ifndef GRAPHICS_STUFF //command prompt mode
	corresp->qapSolver = qapSolver;
#else
	corresp->qapSolver = qapSolver;//(m == 3 && partialID != 4 && strcmp(subject, "david") == 0 && ! pfaust ? QUCOOP : QUCOOPMATCH);//davidpartial3 requires QUCOOP solver initially to prevent sym flip (denseMap() part later always switches to QUCOOPMATCH (currently no switch))
#endif
//corresp->qapSolver = QUCOOP;
//corresp->qapSolver = QMATCH;
//corresp->qapSolver = QGMATCH;
	sprintf(corresp->qapSolverStr, "%s", (corresp->qapSolver == QUCOOP ? "qucoop" : "qucoopmatch"));
	if (corresp->qapSolver == QMATCH)
		sprintf(corresp->qapSolverStr, "qmatch");

	//fill Q matrix for the quadratic assignment problem (QAP); must come after Correspondence constructor where i normalize geodesic[] to [0, 1]
	mesh1->fillQ(nPathSamples + 1, areaWeight, dissimType != PAIRWISENOPATCH), mesh2->fillQ(nPathSamples + 1, areaWeight, dissimType != PAIRWISENOPATCH);
#ifdef GEODIFFTHRESHOLDING
	cout << "PARAMS: radius = " << mesh1->radius << ", nRingAreas = " << mesh1->nRingAreas << ", nPathSamples = " << nPathSamples+1 << ", unaryWeight = " << areaWeight << ", maxGeoDiff = " << maxGeoDiff << (mesh1->weightByGC ? ", wByGC" : "") << ", areasRing, qapEnergy: " << dissimType << ", pathSampler: " << nonuniformPathSampling << "\n"; //non-overlapping areasRing[] alternative is overlapping is areas[] for the unary term
#else
	cout << "PARAMS: radius = " << mesh1->radius << ", nRingAreas = " << mesh1->nRingAreas << ", nPathSamples = " << nPathSamples+1 << ", unaryWeight = " << areaWeight << ", barrier = " << (corresp->barrier > 0 ? corresp->barrier : 0) << (mesh1->weightByGC ? ", wByGC" : "") << ", areasRing, qapEnergy: " << dissimType << ", pathSampler: " << nonuniformPathSampling << "\n"; //barrier n/a if negative
#endif
	do
	{
//		if (qap)
			corresp->qapMap(nPathSamples + 2, false); //nPathSamples=-1 sends 1 which is a special flag for Q2D usage; otherwise nPathSamples + 1 gives # path samples and another + 1 gives access to the additional dimension storing pure geodesic length info (for maxGeoDiff filtering or geoDiff-based barrier)
			//corresp->qapMap2();
/*		else if (virtualSamples) //C(10,4)*C(10,4)*4!=1.058.400 trials (2secs) of bruteForceMapQuad() are replaced by 10!=3.628.800 (11secs); just 86.400 trials in bruteForceMapTrip()
			corresp->bruteForceMap(m); //force n samples match w/ m samples and n-m virtual samples (virtual matches are then discarded since virtual samples do not even belong to mesh2 surface)
		else if (m == 3)
			corresp->bruteForceMapTrip(); //compare C(n, 3) mesh1 triplets w/ each mesh2 triplet (instead of forcing n source to match w/ n target, which is often meaningless in partial matching)
		else if (m == 4)
			corresp->bruteForceMapQuad(); //compare C(n, 4) mesh1 quadruples w/ each mesh2 quadruple (instead of forcing n source to match w/ n target, which is often meaningless in partial matching)
		else if (m == 5)
			corresp->bruteForceMapQuin(); //compare C(n, 5) mesh1 quintuples w/ each mesh2 quintuple (instead of forcing n source to match w/ n target, which is often meaningless in partial matching)
		else
			cout << "WARNING: select m=3,4,5\n\n";*/
		if (corresp->mapOK) //definitely true after qapMap() above
			break;
		else
			corresp->maxGeoDiff *= 2, cout << "maxGeoDiff doubled (brute)\n"; //try mapping w/ a relaxed maxGeoDiff in case corresp->mapOK stayed false in this iteration
	} while (true);
	bool roiReady = false;
	if (denseMapType != NODENSEMAP)
		corresp->denseMap(denseMapType), roiReady = true; //unaryWeight = areaWeight changed, i.e., PARAMS print on screen is deprecated now		
	else
		corresp->mapQualityGeo(true); //print matches' individual geodesic distortions and the overall map distortion //corresp->mapQualityPatch(true), corresp->mapQualityVirtual(true); //print using different criteria (..yVirtual() prints the same as ..yPatch())
	if (fmRefinement)
		cout << "before fullMap: ", corresp->mapQualityGT(), corresp->fullMap(roiReady, true); //additional refinement via anchors based on current few coarse matches (denseMapType=NODENSEMAP) or denser c2f matches (denseMapType=C2FQAP||QAP||LAP)
	if (displayFullMap)
		corresp->transferColors(roiReady);
	else
		corresp->compatibleColoring(false, roiReady); //display map related colors/patched paths (roiReady=false) or rois (roiReady=true)
	corresp->bestWorstMatches(), cout << (fmRefinement ? "after fullMap: " : ""), corresp->mapQualityGT(); //after fullMap() refinement
	} //common block hence the bad indentation
	if (corresp)
		corresp->resultToFile(fmRefinement);
#ifdef GRAPHICS_STUFF
	/////////////////// visuals ///////////////////////
	HWND window = SoWin::init(argv[0]);
	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);
	Painter* painter = new Painter();
	SoSeparator * root = new SoSeparator();
	root->ref();	
	//mesh1
	root->addChild( painter->shapeSep(mesh1) );
	//mesh2
	double deltaX = mesh1->avgEdgeLen * (idTgt != 101 ? 210 : 110);
	if (loadDL)
		deltaX = mesh1->avgEdgeLen * 50;//70;
	root->addChild( painter->shapeSep(mesh2, deltaX) );
	//sampling
//	root->addChild( painter->spheresSep(mesh1) ), root->addChild( painter->spheresSep(mesh2, deltaX) ); //samples
	if (! displayFullMap && ! loadDL)
		root->addChild( painter->spheresSep(mesh1, mesh2, deltaX) ); //samples w/ matching colors
//	root->addChild( painter->spheresSep2(mesh1) ), root->addChild( painter->spheresSep2(mesh2, deltaX) ); //path samples
//	root->addChild( painter->spheresSep3(mesh1, false) ); root->addChild( painter->spheresSep3(mesh2, false, deltaX) ); //path vertices (true) or special vertices (false), e.g., localMaximaofGC
//	root->addChild( painter->spheresSep4(mesh1, mesh2, corresp->worst, corresp->best, deltaX) ); //worst and best matches in red and green, respectively
	//miscellaneous
//	root->addChild( painter->sphereSep(mesh1, 24726) ), root->addChild( painter->sphereSep(mesh2, 13, deltaX) );
	//painter->count++;
//	root->addChild( painter->sphereSep(mesh1, 42897   ) ), root->addChild( painter->sphereSep(mesh1, 2490    ) ), root->addChild( painter->sphereSep(mesh1, 36459   ) ), root->addChild( painter->sphereSep(mesh1, 2402    ) ), root->addChild( painter->sphereSep(mesh1, 36612) );
//	root->addChild( painter->sphereSep(mesh2, 5367    , deltaX) ), root->addChild( painter->sphereSep(mesh2, 3902    , deltaX) ), root->addChild( painter->sphereSep(mesh2, 4519    , deltaX) ), root->addChild( painter->sphereSep(mesh2, 4280    , deltaX) ), root->addChild( painter->sphereSep(mesh2, 894, deltaX) );
//root->addChild( painter->sphereSep(mesh1, 34692) ), root->addChild( painter->sphereSep(mesh2, 5234, deltaX) );//6th sampleeeeeeeeeeeeeeeeeeee
	//painter->count++;
//	root->addChild( painter->sphereSep(mesh1, 5800    ) ), root->addChild( painter->sphereSep(mesh1, 36612   ) ), root->addChild( painter->sphereSep(mesh1, 26579   ) ), root->addChild( painter->sphereSep(mesh1, 2402    ) ), root->addChild( painter->sphereSep(mesh1, 25378) );
//	root->addChild( painter->sphereSep(mesh2, 1308    , deltaX) ), root->addChild( painter->sphereSep(mesh2, 3071    , deltaX) ), root->addChild( painter->sphereSep(mesh2, 4794    , deltaX) ), root->addChild( painter->sphereSep(mesh2, 14      , deltaX) ), root->addChild( painter->sphereSep(mesh2, 105, deltaX) );
	//painter->count++;
//	root->addChild( painter->sphereSep(mesh1, 47193) ), root->addChild( painter->sphereSep(mesh2, 47105, deltaX) ); //24877-2175 rightnipple
//	root->addChild( painter->sphereSep(mesh1, 31674) ), root->addChild( painter->sphereSep(mesh2, 47699, deltaX) ); //35186-3902 leftnipple
//	root->addChild( painter->sphereSep(mesh1, 51914 ) ), root->addChild( painter->sphereSep(mesh2, 51917, deltaX) ); //13719 mouth //34692-4337 left shoulder
//	root->addChild( painter->sphereSep(mesh1, 47824) ), root->addChild( painter->sphereSep(mesh2, 31639, deltaX) ); //13719 mouth //34692-4337 left shoulder
	//correspondence
	if (! displayFullMap && ! loadDL)
		root->addChild( painter->matchingLinesSep(mesh1, mesh2, deltaX) );
	//viewer stuff
	viewer->setBackgroundColor(SbColor(1.0f, 1.0f, 1.0f));
	viewer->setSize(SbVec2s(1024, 768));
	viewer->setSceneGraph(root);
	viewer->setTitle("QAmeetsPM");
	//viewer->setDrawStyle(SoWinViewer::STILL, SoWinViewer::VIEW_WIREFRAME_OVERLAY); //VIEW_POINT);//VIEW_HIDDEN_LINE);//VIEW_LINE);
	//viewer->setWireframeOverlayColor( SbColor(0, 0, 0) );
	viewer->setDrawStyle(SoWinViewer::INTERACTIVE, SoWinViewer::VIEW_LINE); //VIEW_POINT);//VIEW_HIDDEN_LINE);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
#endif

	return 0;
}

double setAreaWeight(Mesh* mesh1, Mesh* mesh2, double& avgPairwise)
{
	//returns areaWeight that leads to meaningful contributions of unary and pairwise terms, i.e., avgUnary and avgPairwise in Q are close after this weighting (areaWeight controls the unary term and is suggested to be [2, 10] by dubrovina-kimmel paper when terms are normalized by their max values, which is not the case here)

	double avgUnary = 0.0, nu = 0, np = 0, u = 1;//0.5;//1.0; //make unary term even more (or less) powerful via u > 1 (u < 1)
	size_t n = mesh1->samples.size(), m = mesh2->samples.size();
	for (size_t s = 0; s < n; s++) //compare samples[s] w/ each mesh2.samples[p] to consider each potential unary term value
		for (size_t p = 0; p < m; p++)
		{
			double unary = 0.0;
			for (size_t a = 0; a < mesh1->verts[ mesh1->samples[s] ]->areasRing.size(); a++)
				//unary += pow(mesh1->verts[ mesh1->samples[s] ]->areas[a] - mesh2->verts[ mesh2->samples[p] ]->areas[a], 2.0); //overlapping/accumulated areas
				unary += pow(mesh1->verts[ mesh1->samples[s] ]->areasRing[a] - mesh2->verts[ mesh2->samples[p] ]->areasRing[a], 2.0); //non-overlapping ring-shaped areas
			avgUnary += sqrt(unary), nu++;
//s=p=n;//samples[0] vs mesh2.samples[0] to make a faster but worse approximation (unary421.456 vs. pairwise65.4292 --> unary90.0984 vs. pairwise74.8811 w/ this approx)
		}
	if (nu == 0)
		cout << "WARNING: areasRing empty for areaWeight; div by 0 coming\n\n"; //pathSamples is also empty as well so segmentation fault is coming before div by 0
	for (size_t s = 0; s < n; s++)
		for (size_t t = s + 1; t < n; t++) //compare samples[s] - samples[t] w/ each mesh2.samples[p] - mesh.samples[q] pair to consider each potential pairwise term value
				for (size_t p = 0; p < m; p++)
					for (size_t q = p + 1; q < m; q++)
					{
						double pairwise = 0.0;
						size_t pathSize = max(mesh1->verts[ mesh1->samples[s] ]->pathSamples[t].size(), mesh2->verts[ mesh2->samples[p] ]->pathSamples[q].size()); //10 vs. 9 samples uses 10 here so that the 10th.area vs. 0 makes a nonzero l2 contribution to the pairwise term
						for (size_t a = 0; a < pathSize; a++) //l2-norm comparison of descriptors desc1 and desc2 along mesh1 & mesh2 paths
						{
							double desc1 = (a < mesh1->verts[ mesh1->samples[s] ]->pathSamples[t].size() ? mesh1->verts[ mesh1->verts[ mesh1->samples[s] ]->pathSamples[t][a] ]->area : 0.0),
								   desc2 = (a < mesh2->verts[ mesh2->samples[p] ]->pathSamples[q].size() ? mesh2->verts[ mesh2->verts[ mesh2->samples[p] ]->pathSamples[q][a] ]->area : 0.0);
							pairwise += pow(desc1 - desc2, 2.0);
						}
						avgPairwise += sqrt(pairwise), np++; //l2-norm
//s=t=p=q=n;//samples[0]-samples[1] vs mesh2.samples[0]-mesh2.samples[1] to make a faster but worse approximation (unary421.456 vs. pairwise65.4292 --> unary90.0984 vs. pairwise74.8811 w/ this approx)
					}
	avgUnary /= nu, avgPairwise /= (np > 0 ? np : 1); //call-by-reference
	cout << "avgUnary: " << avgUnary << " vs. avgPairwise: " << avgPairwise << " before coarse " << m << "x" << m << " map\n"; //" w/ " << nu << " & " << np << " counts
	return u * (m - 1) * avgPairwise / avgUnary; //w/o u & m-1: avgUnary=400, avgPairwise=80 --> areaWeight becomes 80/400=0.2 so that unary and pairwise units are almost equal when added in objective areaWeight*unary + pairwise
												 //via Q matrix (barrier will be added later); *m-1 'cos 1 unary + m-1 pairwise are added from each Q row so make unary more powerful via *m-1
}
