#ifndef _CORRESPONDENCE_H_
#define _CORRESPONDENCE_H_

#include "Mesh.h"

//1: pointwise only, 2: pairwise only (3: patchless pairwise only), 4: pointwise + pairwise combination for the dissimilarity function in mapQualityPatch/Virtual(); pairwise requires pathing and patching (patchless pairwise just compares geodesics so no pathing/patching there)
#define POINTWISE			4 //set of area rings (path/patch disabled)
#define PAIRWISE			3 //path/patch active
#define PAIRWISENOPATCH		2 //path/patch disabled
#define POINTWISEPAIRWISE	1 //path/patch active
//#define GEODIFFTHRESHOLDING //use geoDiff as a threshold to filter out very bad matches (defined) or as a weight for the pairwise term (undefined)
#define LAP			1 //dense map via linear assignment problem (lap)
#define QAP			2 //dense map via quadratic assignment problem (qap)
#define C2FQAP		3 //dense map via coarse-to-fine qap
#define NODENSEMAP	-1 //no dense map just return the coarse initial map (that defines the initial patches)
#define QUCOOP		2 //qap solver is qucoop
#define QUCOOPMATCH	1 //qap solver is qucoopmatch
#define QMATCH		3 //qap solver is qmatch
#define QGMATCH		4 //qap solver is qgm
//#define CANCELONMANYTO1 //define this to cancel 2+ matches targeting the same mesh2 sample (manyTo1 case); undefine to keep only the best/minDistortion of those 2+ matches

class Correspondence
{
public:
	Mesh* mesh1, * mesh2; //shape correspondence b/w these 2 meshes
	double distortion, maxGeoDiff, g, areaWeight, barrier;
	bool mapOK, noBarrier, printQAP, brute, multipleC2Fs;
	int dissimType, qapSolver, nDenseSamples, nSamplesPerSubregion, nMaxLevels, nQAPs, ** perms5, best, worst; //best and worst matches for visualization only
	char dType[250], qapSolverStr[250];
	vector< int > trimmedSamples; //list of samples whose matches are trimmed
	Correspondence(Mesh* m1, Mesh* m2, int disType = 0, double geo = 0, double area = 0, double ap = 0, int nDense = 0, int nSubrSamples = 0, int nMaxLev = 0) : mesh1(m1), mesh2(m2), mapOK(false), dissimType(disType), maxGeoDiff(geo), areaWeight(area), nDenseSamples(nDense), nSamplesPerSubregion(nSubrSamples), nMaxLevels(nMaxLev), brute(false), multipleC2Fs(false), printQAP(true)
	{
		dissimTypeString();
		best = worst = -1, g = 1, barrier = 0.0;
#ifndef GEODIFFTHRESHOLDING //normalize mesh1/2.geodesic[] in [0, 1] as per tanh/cubic function request (maxGeoDiff'll not be used in this GEODIFFTHRESHOLDING-undefined mode so no update)
		noBarrier = false; //true makes no barrier stuff at all (achieved by setting barrier=0; unitGeoScaling() redundant if noBarrier=true but anyway)
		if (geo > 0) //geo=0 means fileFPSforQM is active and hence maxGeoDist is undefined so don't bother doing unitGeoScaling()
			g = mesh1->maxGeoDist/*max(mesh1->maxGeoDist, mesh2->maxGeoDist)*/, mesh1->unitGeoScaling(g), mesh2->unitGeoScaling(g); //max is always mesh1.maxGeoDist of the complete mesh1 model (if 2 models are complete then max is necessary but not that critical; actually pfaust partial models may give a bigger maxGeoDist due to lack of cropped regions/paths so remove max())
		if (! noBarrier && (dissimType == PAIRWISE || dissimType == POINTWISEPAIRWISE)) //average of possible valid pairwise term values affects barrier term since barrier term is added to each pairwise term (ap = 0 means this value is n/a so use geodesic based barrier g*10 verified by observation)
			barrier = (ap > 0 ? ap * 1 : g * 10); //make barrier negative to use geoDiff as a multiplier of pairwise term, not as a barrier term; make barrier 0 to completely disable barrier stuff
//barrier = 0; //no barrier at all
#endif
		//fread permutations from file to perms5[][] mtrx in case bruteForceMapQuinC2F() is used which happens if brute=true
		if (brute)
		{
			int n = 5;//10; //use 10! = 3628800 (instead of 5! = 120) for stress test: took 209 secs to do a single 10x10 brute-force map (10! = 3628800) and i need many of them during c2f so this is simply intractable
			FILE* fPtr;
			char fName[250];
			sprintf(fName, "permutations\\%d!.txt", n); //obtained using matlab's perms function, e.g., perms(0:3) for 4!.txt	
			if (! (fPtr = fopen(fName, "r")))
			{
				cout << fName << " not found\n";
				exit(0);
			}
			//matlab perms(0:10) took 1 sec and fprinting the resulting 10! entries to the 10!.txt file took 2 mins
			perms5 = new int*[(n == 5 ? 120 : 3628800)]; //naming sucks if n=10 'cos perms5[] would store perms10[] then
			for (int i = 0; i < (n == 5 ? 120 : 3628800); i++) //r! reorderings/permutations
			{
				perms5[i] = new int[n];
				for (int j = 0; j < n; j++)
					fscanf(fPtr, "%d\t", &perms5[i][j]);
			}
			fclose(fPtr);
		}
	};
	void bruteForceMap(size_t m);
	void bruteForceMapTrip();
	void bruteForceMapQuad();
	void bruteForceMapQuin();
	bool bruteForceMapQuinC2F();
	void qapMap(int k, bool dense, bool unmatchedCall = false);
	void qapMap2(int nMaxIters = 50, int nSweeps = 100, double beta_min = 0.1, double beta_max = 5.0);	
	int mapQualityGeo(bool print = false);
	int mapQualityPatch(bool print = false);
	int mapQualityVirtual(bool print = false);
	void mapQualityGT();
	void compatibleColoring(bool agdOnly, bool roiReady);
	void denseMap(int mapType, bool fileFPS = false, bool trim = true);
	void denseMapLargeN(size_t k);
	void fullMap(bool roiReady, bool forceSample = false);
	bool printMap(clock_t t, char* mapType);
	void bestWorstMatches();
	void resultToFile(bool full);
	void transferColors(bool roiReady);
private:
	void dissimTypeString();
	double distortionByPermutation(const int* perm, const int& n, const int& m);
	int trimMap(double f, const vector< int >& samples, bool coarseMapMatchesExempt, bool patchCostDistortion, bool print = false);
	void matchUnmatcheds(int nPathSamples, int nPathSamplesExtra);
	void c2fHelper(int src, int tgt, double radius, bool curvFPS, int nPathSamples, int nPathSamplesExtra, size_t cSize, bool patchCostDistortion, int& qual, vector<int>& rd, bool second);
};

#endif
