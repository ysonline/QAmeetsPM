#include "Correspondence.h"
//#include "BloMatch.h" //needed due to blossomMap() below
#include "qucoop.h" //needed due to qapMap() below
#include "qucoopmatch.h" //alternative to qucoop.h
#include "qmatch.h" //needed due to qapMap() below
#include "qgm.h" //needed due to qapMap() below

inline int factorial(int k)
{
	if (k > 10)
		return -1; //don't bother with large numbers
	int res = 1;
	for (int i = k; i > 0; i--)
		res *= i;
	return res;
}

void Correspondence::bruteForceMapTrip()
{
	//same as bruteForceMap() except each C(n, 3) source triplet is compared w/ each C(n, 3) target triplet here (so output is a 3x3 bijective map)

	clock_t t = clock();
	cout << "Brute-force triplet map computation..\n";
	int nb1 = 0, nb2 = 0, n1 = (int) mesh1->samples.size(), n2 = (int) mesh2->samples.size(); //n1=n2 in my case but this equality is not important for this function
	if (n1 < 3 || n2 < 3)
		cout << "WARNING: insufficient # samples for triplets\n\n";
	for (int i1 = 0; i1 < n1; i1++)
		nb1 += (! mesh1->verts[ mesh1->samples[i1] ]->border);
	for (int i2 = 0; i2 < n2; i2++)
		nb2 += (! mesh2->verts[ mesh2->samples[i2] ]->border);
	if (nb1 < 3 || nb2 < 3 && mesh1->borderTreat == 2) //borderTreat=1 already called doubleBorderAreas() that doubled the areas so cannot switch to =0 in =1 mode
		cout << "WARNING: insufficient # non-border samples so switching to borderTreat = 0\n\n", mesh1->borderTreat = mesh2->borderTreat = 0; //cannot switch to borderTreat=1 'cos areas() already called
	//fread permutations from file to perms[][] mtrx
	FILE* fPtr;
	char fName[] = "permutations\\3!.txt"; //obtained using matlab's perms function, e.g., perms(0:2) for 3!.txt
	int factR = factorial(3); 
	if (! (fPtr = fopen(fName, "r")))
	{
		cout << fName << " not found\n";
		exit(0);
	}
	//matlab perms(0:10) took 1 sec and fprinting the resulting 10! entries to the 10!.txt file took 2 mins
	int** perms = new int*[factR];
	for (int i = 0; i < factR; i++) //r! reorderings/permutations
	{
		perms[i] = new int[3];
		for (int j = 0; j < 3; j++)
			fscanf(fPtr, "%d\t", &perms[i][j]);
	}
	fclose(fPtr);

	//use permutations to perform the brute-force map computation
	vector< int > optimalMap;
	double minDistortion = INF;
	//for each triplet in source mesh1
	for (int i1 = 0; i1 < n1; i1++)
		for (int j1 = i1 + 1; j1 < n1; j1++) //each pair
			for (int k1 = j1 + 1; k1 < n1; k1++) //each triplet
			{
				//mesh1.samples[i1], mesh1.samples[j1], mesh1.samples[k1] versus the mesh2 triplet below
				vector< int > mesh1Samples;
				mesh1Samples.push_back(mesh1->samples[i1]), mesh1Samples.push_back(mesh1->samples[j1]), mesh1Samples.push_back(mesh1->samples[k1]);
				for (int i2 = 0; i2 < n2; i2++)
					for (int j2 = i2 + 1; j2 < n2; j2++) //each pair
						for (int k2 = j2 + 1; k2 < n2; k2++) //each triplet
						{
							//compare the current mesh1 triplet w/ each 3!=6 permutation of this current mesh2 triplet: mesh2.samples[i2], mesh2.samples[j2], mesh2.samples[k2]
							vector< int > mesh2Samples;
							mesh2Samples.push_back(mesh2->samples[i2]), mesh2Samples.push_back(mesh2->samples[j2]), mesh2Samples.push_back(mesh2->samples[k2]);
							for (int i = 0; i < factR; i++) //factR = 6 here for sure
							{
								vector< int > tmpMap; //2 consecutive entries (mesh1-mesh2) gives 1 match of the entire map/correspondence
								int idx1 = 0, idx2 = 0;
								for (int y = 0; y < 6; y++) //tmpMap will always have 3 matches in this function so 2*3 = 6 iterations (3 for mesh1 side 3 for mesh2 side)
									if (y % 2 == 0) //even idxs from the fixed mesh1 side
										tmpMap.push_back(mesh1Samples[idx1++]);
									else			//odd idxs for permutations from mesh2 side
									{
										tmpMap.push_back(mesh2Samples[ perms[i][idx2] ]);
										mesh1->verts[ mesh1Samples[idx1-1] ]->matchIdx = mesh2Samples[ perms[i][idx2++] ]; //mapQuality() needs matchIdx values (-1 'cos idx1 was already incremented)
									}
								if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
									mapQualityPatch(); //global distortion is updated
								else
									mapQualityGeo(); //global distortion is updated
								if (distortion < minDistortion)
									minDistortion = distortion, optimalMap = tmpMap;//, cout << i1 << " " << j1 << " " << k1 << " vs. " << i2 << " " << j2 << " " << k2 << endl;
								mesh1->verts[ mesh1Samples[0] ]->matchIdx = mesh1->verts[ mesh1Samples[1] ]->matchIdx = mesh1->verts[ mesh1Samples[2] ]->matchIdx = -1; //reset to unmatched for the next iteration
							}
						}
			}
	//optimalMap to matchIdxs
	for (int i1 = 0; i1 < n1; i1++)
		mesh1->verts[ mesh1->samples[i1] ]->matchIdx = -1; //matchIdx have the very last configuration above so reset them to -1

/*optimalMap[0] = mesh1->samples[5], optimalMap[1] = mesh2->samples[6]; //lhand-lhand
optimalMap[2] = mesh1->samples[4], optimalMap[3] = mesh2->samples[4]; //rhand-rhand
optimalMap[4] = mesh1->samples[1], optimalMap[5] = mesh2->samples[3];//ground-truth for some particular pair*/
	for (int s = 0; s < 3; s++) //optimalMap will always have 3 matches in this function
		mesh1->verts[ optimalMap[2*s] ]->matchIdx = optimalMap[2*s + 1], mesh2->verts[ optimalMap[2*s + 1] ]->matchIdx = optimalMap[2*s]; //dual operation
	//use the latest/optimal matchIdx values to assign individual v.distortion scores for printing only
	if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
		mapQualityPatch(); //global distortion is updated
	else
		mapQualityGeo(); //global distortion is updated
	for (int s = 0; s < 3; s++) //optimalMap will always have 3 matches in this function
	{
		cout << optimalMap[2*s] << " -bruteTrip- " << optimalMap[2*s + 1] << "\t" << mesh1->verts[ optimalMap[2*s] ]->distortion << "\t";
		if (mesh1->verts[ optimalMap[2*s] ]->patchGC == 0.0) //not in use so no print
			cout << endl;
		else
			cout << abs(mesh1->verts[ optimalMap[2*s] ]->patchGC - mesh2->verts[ optimalMap[2*s + 1] ]->patchGC) << "\n"; //difference b/w patch gc also printed after the v.distortion score
	}
	//recapture memory
	for (int i = 0; i < factR; i++)
		delete [] perms[i];
	delete [] perms;
	cout << ">>>>>>>> " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs for matching samples (bruteForceTriplet distortion=" << dType << ": " << distortion << ")\n";
	mapOK = true;
}

int Correspondence::mapQualityGeo(bool print)
{
	//quantifies the quality of the current map/correspondence (implied by matchIdx values) in terms of geodesic distortion, by updating global distortion

	size_t n = mesh1->samples.size(), nMatches = 0;
	distortion = 0.0;
	for (size_t s = 0; s < n; s++)
		if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
		{
			double nAdds = 0.0, individual = 0.0; //individual distortion of this match
			for (size_t t = 0; t < n; t++) //traverse samples[s] - samples[s].matchIdx match above over all matches to get its individual distortion
				if (t != s && mesh1->verts[ mesh1->samples[t] ]->matchIdx != -1) //current matches being traversed: samples[t] - samples[t].matchIdx where t != s
				{
					double val = mesh1->verts[ mesh1->samples[s] ]->geodesic[ mesh1->samples[t] ] - mesh2->verts[ mesh1->verts[ mesh1->samples[s] ]->matchIdx ]->geodesic[ mesh1->verts[ mesh1->samples[t] ]->matchIdx ];
//					individual += (val * val); //use squared L2-norm at all times while comparing distance/geodesic[] vectors as done in Mobius Voting and Biharmonic Distance papers
					individual += abs(val); //disable += val^2 above for this l1-norm effect (use l1 norm for consistency w/ the Q2D-based Q used by qucoop; not an exact consistency)
					nAdds++;
				}
//			individual = sqrt(individual);
//			individual /= (nAdds != 0.0 ? nAdds : 1.0); //always same size maps are evaluated in this project so no division to be consistent with qucoop.objective()
			distortion += individual;
			nMatches++;
			mesh1->verts[ mesh1->samples[s] ]->distortion = individual;
//cout << individual << " added to distortion due to " << mesh1->samples[s] << " - " << mesh1->verts[ mesh1->samples[s] ]->matchIdx << " traversal over all other " << nAdds << " matches\n";
		}
//	if (nMatches > 0.0) //always same size maps are evaluated in this project so no division to be consistent with qucoop.objective()
//		distortion /= nMatches; //average of l2-norms
//	else
//		distortion = INF;
	if (print)
	{
		for (size_t s = 0; s < n; s++)
			if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
				cout << mesh1->samples[s] << " - " << mesh1->verts[ mesh1->samples[s] ]->matchIdx << "\t" << mesh1->verts[ mesh1->samples[s] ]->distortion << "puregeo\n";
		cout << nMatches << "-size map distortion=puregeo: " << distortion << "\n"; //corresponds to PAIRWISENOPATCH
	}
	return nMatches;
}

bool printOnce2 = true;
int Correspondence::mapQualityPatch(bool print)
{
	//same as mapQualityGeo() except more than lengths are compared here using corresponding pathSamples (at the same index by definition)

	size_t n = mesh1->samples.size(), m = mesh2->samples.size(), nMatches = 0;
	distortion = 0.0;
	for (size_t s = 0; s < n; s++)
		if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
		{
			double individual = 0.0; //individual distortion of this match
			int ss = mesh1->samples[s], ssm = mesh1->verts[ mesh1->samples[s] ]->matchIdx;
			if ((mesh1->verts[ss]->border || mesh2->verts[ssm]->border) && mesh1->borderTreat == 2) //treating border samples in way2 by excluding them during optimal map search
			{
				distortion = INF;
				return nMatches; //exclude border samples during optimal map search
			}
			if (dissimType == PAIRWISE || dissimType == POINTWISEPAIRWISE) //pairwise term or pointwise+pairwise terms for individual distortion
				for (size_t t = 0; t < n; t++) //traverse samples[s] - samples[s].matchIdx match above over all matches to get its individual distortion
					if (t != s && mesh1->verts[ mesh1->samples[t] ]->matchIdx != -1) //current matches being traversed: samples[t] - samples[t].matchIdx where t != s
					{
						int ttm = mesh1->verts[ mesh1->samples[t] ]->matchIdx, k = -1;
						double geoDiff = abs(mesh1->verts[ss]->geodesic[ mesh1->samples[t] ] - mesh2->verts[ssm]->geodesic[ttm]), overallDissim = 0.0; //overall dissimilarity b/w s-t and s.matchIdx-t.matchIdx paths; rename them as ss-tt and ssm-ttm here for readability
#ifdef GEODIFFTHRESHOLDING
						if (geoDiff > maxGeoDiff) //very different length paths imply high overall dissimilarity
							overallDissim = INF;
						else
#endif
							for (size_t kk = 0; kk < m; kk++) //learn the index of ttm in samples[] so i can use it as an index to pathSamples below
								if (mesh2->samples[kk] == ttm)
								{
									k = kk;
									break;
								}
						if (k != -1) //k=-1 means overallDissim = INF for sure so don't recompute overallDissim for this current bad match pair
						{
							size_t pathSize = max(mesh1->verts[ss]->pathSamples[t].size(), mesh2->verts[ssm]->pathSamples[k].size()); //size difference rarely happens when samples are too close; 10 vs. 9 samples then use 10 here so that the 10th.area vs. 0 makes a nonzero l2 contribution to overallDissim (consistent w/ with qucoop.objective())
							for (size_t p = 0; p < pathSize; p++) //l2-norm comparison of descriptors desc1 and desc2 along mesh1 & mesh2 paths
							{
								double desc1 = (p < mesh1->verts[ss]->pathSamples[t].size() ? mesh1->verts[ mesh1->verts[ss]->pathSamples[t][p] ]->area : 0.0), desc2 = (p < mesh2->verts[ssm]->pathSamples[k].size() ? mesh2->verts[ mesh2->verts[ssm]->pathSamples[k][p] ]->area : 0.0); //worse alternative: patchGC
								overallDissim += pow(desc1 - desc2, 2.0);
							}
							overallDissim = sqrt(overallDissim); //l2-norm
#ifdef GEODIFFTHRESHOLDING
							individual += overallDissim; //no geoDiff-based weighting 'cos geoDiff is already used by maxGeoDiff filtering above
#else						//geoDiff in [0, 1] for sure
							if (barrier < 0)
								//individual += overallDissim * geoDiff; //linear weighting y=x line where x is geoDiff
								individual += overallDissim * geoDiff * (geoDiff*geoDiff - 3*geoDiff + 3); //non-linear cubic weighting that stays above y=x line where x is geoDiff (endpoints 0 to 1)
								//individual += overallDissim * sqrt(geoDiff); //non-linear weighting that stays above y=x line for x < 1 (endpoints 0 to 1)
								//individual += overallDissim * tanh(geoDiff); //non-linear weighting that stays below y=x line for x < 1 (endpoints 0 to 0.76) (very close to 1-e^-x which ends at 0.63)
								//individual += overallDissim * (geoDiff*geoDiff); //non-linear weighting that stays below y=x line for x < 1 (endpoints 0 to 1)
								//individual += overallDissim; //no weighting
							else
							{
								individual += (overallDissim + barrier * geoDiff * (geoDiff*geoDiff - 3*geoDiff + 3)); //geoDiff (from source i to j vs. from target k to l) is added to the unweighted pairwise term (overallDissim)

if (print) cout << mesh1->coarseMap.size() << ": " << mesh1->verts[ mesh1->samples[t] ]->geodesic.size() << " orrrrr " << mesh2->verts[ttm]->geodesic.size() << endl;
								for (size_t e = 0; e < mesh1->coarseMap.size() / 2; e++)
									geoDiff = abs(mesh1->verts[ mesh1->samples[t] ]->geodesic[ mesh1->coarseMap[2*e] ] - mesh2->verts[ttm]->geodesic[ mesh1->coarseMap[2*e + 1] ]), //from source j to coarseMap source sides vs. from target l to coarseMap target sides
									individual += (barrier * geoDiff * (geoDiff*geoDiff - 3*geoDiff + 3)); //geoDiffs to coarseMap entries are added too
							}
#endif
						}
						else
							individual += overallDissim; //overallDissim=INF here
						//nAdds++;
					} //individual /= (nAdds != 0.0 ? nAdds : 1.0); //always same size maps are evaluated in this project so no division to be consistent with qucoop.objective(); many pairwise terms plus 1 unary term below makes the latter ineffective so made unary powerful via large areaWeight
			if (dissimType == POINTWISE || dissimType == POINTWISEPAIRWISE) //pointwise term or pointwise+pairwise terms for individual distortion
			{
				//areasRing descriptor difference b/w ss and ssm should affect the individual distortion of ss-ssm too; although area above and areasRing[a] below are same units i still need to weight them for a meaningful combination				
				double descDissim = 0.0; //descDissim is pointwise (not pairwise) descriptor dissimilarity
				for (size_t a = 0; a < mesh1->verts[ss]->areasRing.size(); a++)
					//descDissim += pow(mesh1->verts[ss]->areas[a] - mesh2->verts[ssm]->areas[a], 2.0); //overlapping/accumulated areas
					descDissim += pow(mesh1->verts[ss]->areasRing[a] - mesh2->verts[ssm]->areasRing[a], 2.0); //non-overlapping ring-shaped areas
				descDissim = sqrt(descDissim);
				if (dissimType == POINTWISE)
				{
					individual = descDissim; //no combination just the pointwise descDissim term
					if (printOnce2)
						cout << "WARNING: Q is filled w/ gc but mapQuality() checks areasRing so it's OK to get different objective & distortion prints\n\n", printOnce2 = false;
				}
				else
					individual += (descDissim * areaWeight); //combining pairwise term individual w/ unary/pointwise term descDissim using areaWeight
			}
			distortion += individual;
			nMatches++;
			mesh1->verts[ mesh1->samples[s] ]->distortion = individual;
		}
//	if (nMatches > 0.0) //always same size maps are evaluated in this project so no division to be consistent with qucoop.objective()
//		distortion /= nMatches; //average of l2-norms
//	else
//		distortion = INF;
	if (print)
	{
		for (size_t s = 0; s < n; s++)
			if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
				cout << mesh1->samples[s] << " - " << mesh1->verts[ mesh1->samples[s] ]->matchIdx << "\t" << mesh1->verts[ mesh1->samples[s] ]->distortion << dType << endl;
		cout << nMatches << "-size map dissimType=" << dType << " distortion: " << distortion << ", avg: " << distortion / nMatches << "\n";
	}
	return nMatches;
}

int Correspondence::mapQualityVirtual(bool print)
{
	//same as mapQualityPatch() except virtual samples (id=-2) are handled here

	size_t n = mesh1->samples.size(), m = mesh2->samples.size(), nMatches = 0, ok = 1;
	distortion = 0.0;
	for (size_t s = 0; s < n; s++)
		if (mesh1->verts[ mesh1->samples[s] ]->matchIdx > -1) //-1 (no match) and -2 (virtual sample match) are discarded
		{
			double individual = 0.0; //individual distortion of this match
			int ss = mesh1->samples[s], ssm = mesh1->verts[ mesh1->samples[s] ]->matchIdx;
			if ((mesh1->verts[ss]->border || mesh2->verts[ssm]->border) && mesh1->borderTreat == 2 && dissimType != PAIRWISENOPATCH) //PAIRWISENOPATCH can robustly work w/ border samples 'cos pure geodesic lengths are not affected by the missing/cut data around border samples
			{
				distortion = INF;
				return nMatches; //exclude border samples during optimal map search
			}
			if (dissimType == PAIRWISE || dissimType == POINTWISEPAIRWISE || dissimType == PAIRWISENOPATCH) //pairwise term or pointwise+pairwise terms or patchless pairwise term for individual distortion
				for (size_t t = 0; t < n; t++) //traverse samples[s] - samples[s].matchIdx match above over all matches to get its individual distortion
					if (t != s && mesh1->verts[ mesh1->samples[t] ]->matchIdx > -1) //current matches being traversed: samples[t] - samples[t].matchIdx where t != s (-1 (no match) and -2 (virtual sample match) are discarded)
					{
						int ttm = mesh1->verts[ mesh1->samples[t] ]->matchIdx, k = -1;
						double geoDiff = abs(mesh1->verts[ss]->geodesic[ mesh1->samples[t] ] - mesh2->verts[ssm]->geodesic[ttm]), overallDissim = 0.0; //overall dissimilarity b/w s-t and s.matchIdx-t.matchIdx paths; rename them as ss-tt and ssm-ttm here for readability
#ifdef GEODIFFTHRESHOLDING
						if (geoDiff > maxGeoDiff) //very different length paths imply high overall dissimilarity
							overallDissim = INF;
						else
#endif
							for (size_t kk = 0; kk < m; kk++) //learn the index of ttm in samples[] so i can use it as an index to pathSamples below
								if (mesh2->samples[kk] == ttm)
								{
									k = kk;
									break;
								}
						if (k == -1) //different size pathSamples imply different paths so high dissimilarity (geoDiff handles this condition and forcing it here may be wrong for nonuniform triangulations); k==-1 condition may be true due to geoDiff
							individual += overallDissim, ok = -1; //geoDiff test made an INF so sampling not OK (if sampling is ok/compatible then there'll eventually be a permutation w/o INF so ok=1 will survive and notify samplingOK/mapOK global) (overallDissim=INF here)
						else if (dissimType == PAIRWISE || dissimType == POINTWISEPAIRWISE)
						{
							size_t pathSize = max(mesh1->verts[ss]->pathSamples[t].size(), mesh2->verts[ssm]->pathSamples[k].size()); //size difference rarely happens when samples are too close; 10 vs. 9 samples then use 10 here so that the 10th.area vs. 0 makes a nonzero l2 contribution to overallDissim (consistent w/ with qucoop.objective())
							for (size_t p = 0; p < pathSize; p++) //l2-norm comparison of descriptors desc1 and desc2 along mesh1 & mesh2 paths
							{
								double desc1 = (p < mesh1->verts[ss]->pathSamples[t].size() ? mesh1->verts[ mesh1->verts[ss]->pathSamples[t][p] ]->area : 0.0), desc2 = (p < mesh2->verts[ssm]->pathSamples[k].size() ? mesh2->verts[ mesh2->verts[ssm]->pathSamples[k][p] ]->area : 0.0); //worse alternative: patchGC
								overallDissim += pow(desc1 - desc2, 2.0);
							}
							overallDissim = sqrt(overallDissim); //l2-norm
#ifdef GEODIFFTHRESHOLDING
							individual += overallDissim; //no geoDiff-based weighting 'cos geoDiff is already used by maxGeoDiff filtering above
#else						//geoDiff in [0, 1] for sure
							if (barrier < 0)
								//individual += overallDissim * geoDiff; //linear weighting y=x line where x is geoDiff
								individual += overallDissim * geoDiff * (geoDiff*geoDiff - 3*geoDiff + 3); //non-linear cubic weighting that stays above y=x line where x is geoDiff (endpoints 0 to 1)
								//individual += overallDissim * sqrt(geoDiff); //non-linear weighting that stays above y=x line for x < 1 (endpoints 0 to 1)
								//individual += overallDissim * tanh(geoDiff); //non-linear weighting that stays below y=x line for x < 1 (endpoints 0 to 0.76) (very close to 1-e^-x which ends at 0.63)
								//individual += overallDissim * (geoDiff*geoDiff); //non-linear weighting that stays below y=x line for x < 1 (endpoints 0 to 1)
								//individual += overallDissim; //no weighting
							else
							{
								individual += (overallDissim + barrier * geoDiff * (geoDiff*geoDiff - 3*geoDiff + 3)); //geoDiff (from source i to j vs. from target k to l) is added to the unweighted pairwise term (overallDissim)
								for (size_t e = 0; e < mesh1->coarseMap.size() / 2; e++)
									geoDiff = abs(mesh1->verts[ mesh1->samples[t] ]->geodesic[ mesh1->coarseMap[2*e] ] - mesh2->verts[ttm]->geodesic[ mesh1->coarseMap[2*e + 1] ]), //from source j to coarseMap source sides vs. from target l to coarseMap target sides
									individual += (barrier * geoDiff * (geoDiff*geoDiff - 3*geoDiff + 3)); //geoDiffs to coarseMap entries are added too
							}
#endif
						}
						else //PAIRWISENOPATCH
							individual += geoDiff; //overallDissim cannot be INF here since k!=-1 for sure (overallDissim=geoDiff here but overallDissim not updated so use += geoDiff instead)
						//nAdds++;
					} //individual /= (nAdds != 0.0 ? nAdds : 1.0); //always same size maps are evaluated in this project so no division to be consistent with qucoop.objective(); many pairwise terms plus 1 unary term below makes the latter ineffective so made unary powerful via large areaWeight
			if (dissimType == POINTWISE || dissimType == POINTWISEPAIRWISE) //pointwise term or pointwise+pairwise terms for individual distortion
			{
				//areasRing descriptor difference b/w ss and ssm should affect the individual distortion of ss-ssm too; although area above and areasRing[a] below are same units i still need to weight them for a meaningful combination				
				double descDissim = 0.0; //descDissim is pointwise (not pairwise) descriptor dissimilarity
				for (size_t a = 0; a < mesh1->verts[ss]->areasRing.size(); a++)
					//descDissim += pow(mesh1->verts[ss]->areas[a] - mesh2->verts[ssm]->areas[a], 2.0); //overlapping/accumulated areas
					descDissim += pow(mesh1->verts[ss]->areasRing[a] - mesh2->verts[ssm]->areasRing[a], 2.0); //non-overlapping ring-shaped areas
				descDissim = sqrt(descDissim);
				if (dissimType == POINTWISE)
				{
					individual = descDissim; //no combination just the pointwise descDissim term
					if (printOnce2)
						cout << "WARNING: Q is filled w/ gc but mapQuality() checks areasRing so it's OK to get different objective & distortion prints\n\n", printOnce2;
				}
				else
					descDissim *= areaWeight, individual += descDissim; //combining pairwise term individual w/ unary/pointwise term descDissim using areaWeight
			}
			distortion += individual;
			nMatches++;
			mesh1->verts[ mesh1->samples[s] ]->distortion = individual;
		}
//	if (nMatches > 0.0) //always same size maps are evaluated in this project so no division to be consistent with qucoop.objective()
//		distortion /= nMatches; //average of l2-norms
//	else
//		distortion = INF;
	if (ok == 1)
		mapOK = true;
	if (print)
	{
		for (size_t s = 0; s < n; s++)
			if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
				cout << mesh1->samples[s] << " - " << mesh1->verts[ mesh1->samples[s] ]->matchIdx << "\t" << mesh1->verts[ mesh1->samples[s] ]->distortion << dType << endl;
		cout << nMatches << "-size map dissimType=" << dType << " distortion: " << distortion << ", avg: " << distortion / nMatches << "\n";
	}
	return nMatches;
}

void Correspondence::mapQualityGT()
{
	//computes ground-truth distortion of the current map for evaluation only

cout<<"no gt eval in public release\n";return;//no gt evaluation in this public release
	double nAdds = 0, maxGT = -INF, minGT = INF, gtDistortion = 0.0;
	worst = best = -1;
//for (int i = 0; i < (int) mesh1->samples.size(); i++) cout << mesh1->samples[i] << ".gtIdx = " << mesh1->verts[ mesh1->samples[i] ]->gtIdx << endl; //not all mesh1 samples have gtIdx info so nAdds < nMatches is possible
	for (int i = 0; i < (int) mesh1->samples.size(); i++)
		if (mesh1->verts[ mesh1->samples[i] ]->gtIdx != -1 && mesh1->verts[ mesh1->samples[i] ]->matchIdx != -1) //vert whose gt is known has a match
		{
			int v2 = mesh1->verts[ mesh1->samples[i] ]->matchIdx;
			if (mesh2->verts[v2]->geodesic.empty()) //v2 is definitely a sample whose geodesics to all other verts are known but if called after fullMap() v2 may not be a sample so add this condition
				continue;
			mesh1->verts[ mesh1->samples[i] ]->gtDistortion = (mesh2->verts[v2]->geodesic[ mesh1->verts[ mesh1->samples[i] ]->gtIdx ] / 1.0);//mesh2->avgEdgeLen);			
			gtDistortion += mesh1->verts[ mesh1->samples[i] ]->gtDistortion;
//cout << mesh1->samples[i] << " - " << v2 << "\t" << mesh1->verts[ mesh1->samples[i] ]->gtDistortion << "gtDistortion\ttrue gt: " << mesh1->verts[ mesh1->samples[i] ]->gtIdx << endl;
			nAdds++;
			if (mesh1->verts[ mesh1->samples[i] ]->gtDistortion > maxGT)
				maxGT = mesh1->verts[ mesh1->samples[i] ]->gtDistortion, worst = mesh1->samples[i];
			if (mesh1->verts[ mesh1->samples[i] ]->gtDistortion < minGT)
				minGT = mesh1->verts[ mesh1->samples[i] ]->gtDistortion, best = mesh1->samples[i];
		}
	if (nAdds > 0)
		cout << "gt-distortion of this map: " << gtDistortion / nAdds << " (avg over " << nAdds << " adds); worst: " << worst << "-" << mesh1->verts[worst]->matchIdx << " " << maxGT << "; best: " << best << "-" << mesh1->verts[best]->matchIdx << " " << minGT << "\n";

	//51458 -qap- 4166 is a match but 51458.gtIdx=-1; 4166.gtIdx=51421 so i use 51458.geodesic[51421] cost over mesh1
	//47444 -qap- 2220 is a match w/ 47444.gtIdx = 2221; 2220.gtIdx=47445 so i could use 2220.geodesic[2221] over mesh2 (loop above) or 47444.geodesic[47445] over mesh2 (loop below) (similar but different costs)
	nAdds = 0, maxGT = -INF, minGT = INF, gtDistortion = 0.0;
	worst = best = -1;	
//for (size_t i = 0; i < mesh2->samples.size(); i++) cout << mesh2->verts[ mesh2->samples[i] ]->gtIdx << ".gtIdx = " << mesh2->samples[i] << endl; //all mesh2 samples have gtIdx info for sure (see loadGT())
	for (size_t i = 0; i < mesh2->samples.size(); i++) //each mesh2 sample has a valid gtIdx so use them to trace distortion over mesh1 (unlike above which traces through mesh2)
		if (mesh2->verts[ mesh2->samples[i] ]->matchIdx != -1 && //vert whose gt is known (always the case for mesh2 vert) has a match
			mesh2->verts[ mesh2->samples[i] ]->gtIdx != -1) //just be safe in case no gt info is available for this dataset
		{
			int v1 = mesh2->verts[ mesh2->samples[i] ]->matchIdx;
			if (mesh1->verts[v1]->geodesic.empty()) //v1 is definitely a sample whose geodesics to all other verts are known but if called after fullMap() v1 may not be a sample so add this condition
				continue;
			mesh1->verts[v1]->gtDistortion = mesh2->verts[ mesh2->samples[i] ]->gtDistortion = (mesh1->verts[v1]->geodesic[ mesh2->verts[ mesh2->samples[i] ]->gtIdx ] / 1.0);//mesh2->avgEdgeLen);			
			gtDistortion += mesh2->verts[ mesh2->samples[i] ]->gtDistortion;
			nAdds++;
			if (mesh1->verts[v1]->gtDistortion > maxGT)
				maxGT = mesh1->verts[v1]->gtDistortion, worst = v1;
			if (mesh1->verts[v1]->gtDistortion < minGT)
				minGT = mesh1->verts[v1]->gtDistortion, best = v1;
		}
	if (nAdds > 0)
		cout << "gt-distortion of this map: " << gtDistortion / nAdds << " (avg over " << nAdds << " adds); worst: " << worst << "-" << mesh1->verts[worst]->matchIdx << " " << maxGT << "; best: " << best << "-" << mesh1->verts[best]->matchIdx << " " << minGT << " (over mesh1)\n";
}

void Correspondence::resultToFile(bool full)
{
	//fprints the current map to file as 2 ints per row where int1 is source sample id and int2 is matching target sample id

//return;

	char fName[250], fName2[250];
	sprintf(fName, "output\\%d-%d%s m=%d n=%d%s.txt", mesh1->id, mesh2->id, (full ? " full" : ""), nSamplesPerSubregion, nDenseSamples, (brute ? " brute" : ""));
	FILE* fPtr = fopen(fName, "w");
	int nMatches = 0;
//	for (int i = 0; i < (int) mesh1->samples.size(); i++)
//		if (mesh1->verts[ mesh1->samples[i] ]->matchIdx != -1)
//			fprintf(fPtr, "%d\t%d\n", mesh1->samples[i], mesh1->verts[ mesh1->samples[i] ]->matchIdx), n++;
	for (int i = 0; i < (int) mesh1->verts.size(); i++)
		if (mesh1->verts[i]->matchIdx != -1)
			fprintf(fPtr, "%d\t%d\n", i, mesh1->verts[i]->matchIdx), nMatches++;
	cout << "# matches fprinted to " << fName << " is ";
	//csv format for the texture transfer script, e.g., use the same index 0 even if matchIdx=-1 (0 is arbitary; just make it constant in the range [0,V] across the file)	
	sprintf(fName, "output\\%d-%d%s m=%d n=%d%s.csv", mesh1->id, mesh2->id, (full ? " full" : ""), nSamplesPerSubregion, nDenseSamples, (brute ? " brute" : "")), sprintf(fName2, "output\\%d-%d%s m=%d n=%d%s.csv", mesh2->id, mesh1->id, (full ? " full" : ""), nSamplesPerSubregion, nDenseSamples, (brute ? " brute" : ""));
	FILE* fPtrCsv = fopen(fName, "w"), * fPtrCsv2 = fopen(fName2, "w");
	for (int i = 0; i < (int) mesh1->verts.size(); i++) //mesh1 to mesh2
		if (i == mesh1->verts.size() - 1) //no comma after the last match
			fprintf(fPtrCsv, "%d\n", (mesh1->verts[i]->matchIdx != -1 ? mesh1->verts[i]->matchIdx : 0));				
		else //comma after each match
			fprintf(fPtrCsv, "%d,", (mesh1->verts[i]->matchIdx != -1 ? mesh1->verts[i]->matchIdx : 0));	
	for (int i = 0; i < (int) mesh2->verts.size(); i++) //mesh2 to mesh1
		if (i == mesh2->verts.size() - 1) //no comma after the last match
			fprintf(fPtrCsv2, "%d\n", (mesh2->verts[i]->matchIdx != -1 ? mesh2->verts[i]->matchIdx : 0));
		else //comma after each match
			fprintf(fPtrCsv2, "%d,", (mesh2->verts[i]->matchIdx != -1 ? mesh2->verts[i]->matchIdx : 0));
	fclose(fPtr), fclose(fPtrCsv), fclose(fPtrCsv2), cout << nMatches << "\n";
}

void Correspondence::compatibleColoring(bool agdOnly, bool roiReady)
{
	//determines compatible color values on mesh1 & mesh2 based on the agreed nColors value (applicable to matching projects)

	if (agdOnly || mesh1->samples.empty())
	{
		mesh1->vertColors(4), mesh2->vertColors(4); //1: patch coloring, 2: area ring coloring, 3: roi coloring, others: descriptor coloring, e.g., agd descriptor (bourke colormap)
		return;
	}
	if (roiReady) //region of interest coloring and no further action
	{		
		mesh1->vertColors(3), mesh2->vertColors(3); //1: patch coloring, 2: area ring coloring, 3: roi coloring, others: descriptor coloring (bourke colormap)
		return;
	}
	bool firstMatches = false; //true: first 3 matches, false: random 3 matches provide the patched path visualization
	if (firstMatches) //find the first matched sample on mesh1 (mesh1.m1) and find the mesh2 end of that match (mesh2.m1); do the same for the second and third matched samples too
	{
		for (mesh1->m1 = 0; mesh1->m1 < mesh1->samples.size(); mesh1->m1++)
			if (mesh1->verts[ mesh1->samples[ mesh1->m1 ] ]->matchIdx != -1) //found the first matched sample on mesh1
			{
				for (mesh2->m1 = 0; mesh2->m1 < mesh2->samples.size(); mesh2->m1++)
					if (mesh2->samples[ mesh2->m1 ] == mesh1->verts[ mesh1->samples[ mesh1->m1 ] ]->matchIdx) //found the mesh2 side of the first matched sample on mesh1
						break;
				mesh2->m1 = (mesh2->m1 == mesh2->samples.size() ? 0 : mesh2->m1); //fullMap() may make mesh1.m1.matchIdx a non-sample in which case mesh2.m1 would be invalid so reset it to 0
				break;
			}
		//do the same for second matched and third samples
		for (mesh1->m2 = mesh1->m1 + 1; mesh1->m2 < mesh1->samples.size(); mesh1->m2++)
			if (mesh1->verts[ mesh1->samples[ mesh1->m2 ] ]->matchIdx != -1) //found the second matched sample on mesh1
			{
				for (mesh2->m2 = 0; mesh2->m2 < mesh2->samples.size(); mesh2->m2++)
					if (mesh2->samples[ mesh2->m2 ] == mesh1->verts[ mesh1->samples[ mesh1->m2 ] ]->matchIdx) //found the mesh2 side of the second matched sample on mesh1
						break;
				mesh2->m2 = (mesh2->m2 == mesh2->samples.size() ? 0 : mesh2->m2); //fullMap() may make mesh1.m1.matchIdx a non-sample in which case mesh2.m1 would be invalid so reset it to 0
				break;
			}
		for (mesh1->m3 = mesh1->m2 + 1; mesh1->m3 < mesh1->samples.size(); mesh1->m3++)
			if (mesh1->verts[ mesh1->samples[ mesh1->m3 ] ]->matchIdx != -1) //found the third matched sample on mesh1
			{
				for (mesh2->m3 = 0; mesh2->m3 < mesh2->samples.size(); mesh2->m3++)
					if (mesh2->samples[ mesh2->m3 ] == mesh1->verts[ mesh1->samples[ mesh1->m3 ] ]->matchIdx) //found the mesh2 side of the third matched sample on mesh1
						break;
				mesh2->m3 = (mesh2->m3 == mesh2->samples.size() ? 0 : mesh2->m3); //fullMap() may make mesh1.m1.matchIdx a non-sample in which case mesh2.m1 would be invalid so reset it to 0
				break;
			}
		if (mesh1->m1 >= mesh1->samples.size() || mesh1->m2 >= mesh1->samples.size() || mesh1->m3 >= mesh1->samples.size() || mesh2->m1 >= mesh2->samples.size() || mesh2->m2 >= mesh2->samples.size() || mesh2->m3 >= mesh2->samples.size())
			cout << "WARNING: <3 matches not OK; use at least 3 matches for a nice visualization\n\n";
	}
	else //find a random matched sample on mesh1 (mesh1.m1) and find the mesh2 end of that match (mesh2.m1); do the same for the second and third random matched samples too
	{
		int nTry = 0, nMaxTry = 10000;
		mesh1->m1 = (int) (rand() % mesh1->samples.size()); //random matched sample in [0, size-1) interval
		while (mesh1->verts[ mesh1->samples[ mesh1->m1 ] ]->matchIdx == -1 && nTry++ < nMaxTry)
			mesh1->m1 = (int) (rand() % mesh1->samples.size());
		for (mesh2->m1 = 0; mesh2->m1 < mesh2->samples.size(); mesh2->m1++)
			if (mesh2->samples[ mesh2->m1 ] == mesh1->verts[ mesh1->samples[ mesh1->m1 ] ]->matchIdx) //found the mesh2 side of the first random matched sample on mesh1
				break;
		mesh2->m1 = (mesh2->m1 == mesh2->samples.size() ? 0 : mesh2->m1); //fullMap() may make mesh1.m1.matchIdx a non-sample in which case mesh2.m1 would be invalid so reset it to 0
		//do the same for second matched and third samples
		mesh1->m2 = (int) (rand() % mesh1->samples.size()); //random matched sample in [0, size-1) interval
		while ((mesh1->m1 == mesh1->m2 || mesh1->verts[ mesh1->samples[ mesh1->m2 ] ]->matchIdx == -1) && nTry++ < nMaxTry) //ensure m1 != m2
			mesh1->m2 = (int) (rand() % mesh1->samples.size());
		for (mesh2->m2 = 0; mesh2->m2 < mesh2->samples.size(); mesh2->m2++)
			if (mesh2->samples[ mesh2->m2 ] == mesh1->verts[ mesh1->samples[ mesh1->m2 ] ]->matchIdx) //found the mesh2 side of the second random matched sample on mesh1
				break;
		mesh2->m2 = (mesh2->m2 == mesh2->samples.size() ? 0 : mesh2->m2); //fullMap() may make mesh1.m1.matchIdx a non-sample in which case mesh2.m1 would be invalid so reset it to 0
		mesh1->m3 = (int) (rand() % mesh1->samples.size()); //random matched sample in [0, size-1) interval
		while ((mesh1->m1 == mesh1->m3 || mesh1->m2 == mesh1->m3 || mesh1->verts[ mesh1->samples[ mesh1->m3 ] ]->matchIdx == -1) && nTry++ < nMaxTry) //ensure m1 != m3 and m2 != m3
			mesh1->m3 = (int) (rand() % mesh1->samples.size());
		for (mesh2->m3 = 0; mesh2->m3 < mesh2->samples.size(); mesh2->m3++)
			if (mesh2->samples[ mesh2->m3 ] == mesh1->verts[ mesh1->samples[ mesh1->m3 ] ]->matchIdx) //found the mesh2 side of the third random matched sample on mesh1
				break;
		mesh2->m3 = (mesh2->m3 == mesh2->samples.size() ? 0 : mesh2->m3); //fullMap() may make mesh1.m1.matchIdx a non-sample in which case mesh2.m1 would be invalid so reset it to 0
		if (nTry >= nMaxTry)
			cout << "WARNING: <3 matches not OK; use at least 3 matches for a nice visualization\n\n";
	}
	for (size_t v = 0; v < mesh1->verts.size(); v++) //patched path visualization (end points of the patched path(s) as spheres)
		mesh1->verts[v]->marked = (v == mesh1->samples[ mesh1->m1 ] || v == mesh1->samples[ mesh1->m2 ]);//|| v == s[mesh1->m3]); //true means display it as sphere on screen (for visualization only)
	size_t nColors = 0; //this many smoothly changing colors in colors[]
	double maxGeoToRef1 = -INF, maxGeoToRef2 = -INF;
	if (dissimType == PAIRWISENOPATCH) //pure geodesic distances are used so display geoToReft where reference ref is the endpoints of the first match (mesh1.m1 - mesh2.m1)
	{
		//geodesic to reference samples mesh1.m1 and mesh2.m1
		for (size_t v = 0; v < mesh1->verts.size(); v++)
			mesh1->verts[v]->geoToRef = mesh1->verts[ mesh1->samples[ mesh1->m1 ] ]->geodesic[v] * g; //undo the scaling in the constructor
		for (size_t v = 0; v < mesh2->verts.size(); v++)
			mesh2->verts[v]->geoToRef = mesh2->verts[ mesh2->samples[ mesh2->m1 ] ]->geodesic[v] * g;
		//maxGeoToRef1/2 decision for normalization
		for (size_t i = 0; i < mesh1->verts.size(); i++)
			maxGeoToRef1 = (mesh1->verts[i]->geoToRef > maxGeoToRef1 ? mesh1->verts[i]->geoToRef : maxGeoToRef1);
		for (size_t i = 0; i < mesh2->verts.size(); i++)
			maxGeoToRef2 = (mesh2->verts[i]->geoToRef > maxGeoToRef2 ? mesh2->verts[i]->geoToRef : maxGeoToRef2);
		nColors = (size_t) max(maxGeoToRef1, maxGeoToRef2) + 1; //+1 to handle, e.g., [0,172] range where i don't do -1 during picking below so [172], and not [171], is possible; +1 here mallocs [172] as well (for geoColors=true only)
	}
	else if (dissimType == POINTWISE) //areasRing in use as pointwise descriptors, i.e., no patching, so use vertColors(2) here which uses v.binID values where bins are rings (mesh1.m1 - mesh2.m1 set above used here)
	{
		mesh1->areasForDisplay(mesh1->m1), mesh2->areasForDisplay(mesh2->m1); //binID ready
		mesh1->vertColors(2), mesh2->vertColors(2); //1: patch coloring, 2: area ring coloring, others: descriptor coloring (bourke colormap)
		return;
	}
	else //patch coloring 'cos patch-based PAIRWISE or POINTWISEPAIRWISE in use (m1/2/3 set above will be used by vertColors(1))
	{
		//typical parameter is 1 which colors patches of the samples on the path from samples[m1] to samples[m2] in red (m3 not in use)
		mesh1->vertColors(1), mesh2->vertColors(1); //1: patch coloring, 2: area ring coloring, others: descriptor coloring (bourke colormap)
		return;
	}
	float** colors = new float*[nColors], valMin = 0.0, valMax = (float) nColors - 1, dv = valMax - valMin; //only PAIRWISENOPATCH reaches here
	for (size_t i = 0; i < nColors; i++)
	{
		colors[i] = new float[3];
		double val = i;
		if (val < (valMin + 0.25*dv))
			colors[i][0] = 0.0f, colors[i][1] = (float) (4 * (val - valMin) / dv), colors[i][2] = 1.0f;
		else if (val < (valMin + 0.5*dv))
			colors[i][0] = 0.0f, colors[i][1] = 1.0f, colors[i][2] = (float) (1 + 4 * (valMin + 0.25*dv - val) / dv);
		else if (val < (valMin + 0.75*dv))
			colors[i][0] = (float) (4 * (val - valMin - 0.5*dv) / dv), colors[i][1] = 1.0f, colors[i][2] = 0.0f;
		else
			colors[i][0] = 1.0f, colors[i][1] = (float) (1 + 4 * (valMin + 0.75*dv - val) / dv), colors[i][2] = 0.0f;
/*		//colors above are smooth but uses the range dv which is based on partiality, i.e., hks values under partial matching, so use random colors instead
		for (size_t c = 0; c < 3; c++)
			colors[i][c] = (float) rand() / RAND_MAX;//*/
	}
	//pick mesh1 and mesh2 colors from colors[] according to the index based on hksScalar/geoToRef descriptor value
	for (size_t i = 0; i < mesh1->verts.size(); i++)
		mesh1->verts[i]->color = new float[3], mesh1->verts[i]->color[0] = colors[ (int) mesh1->verts[i]->geoToRef ][0], mesh1->verts[i]->color[1] = colors[ (int) mesh1->verts[i]->geoToRef ][1], mesh1->verts[i]->color[2] = colors[ (int) mesh1->verts[i]->geoToRef ][2];
	for (size_t i = 0; i < mesh2->verts.size(); i++)
		mesh2->verts[i]->color = new float[3], mesh2->verts[i]->color[0] = colors[ (int) mesh2->verts[i]->geoToRef ][0], mesh2->verts[i]->color[1] = colors[ (int) mesh2->verts[i]->geoToRef ][1], mesh2->verts[i]->color[2] = colors[ (int) mesh2->verts[i]->geoToRef ][2];
	//recapture memory
	for (size_t i = 0; i < nColors; i++)
		delete [] colors[i];
	delete [] colors;
}

void Correspondence::bruteForceMapQuad()
{
	//same as bruteForceMap() except each C(n, 4) source quadruple is compared w/ each C(n, 4) target quadruple here (so output is a 4x4 bijective map)

	clock_t t = clock();
	cout << "Brute-force quadruple map computation..\n";
	int nb1 = 0, nb2 = 0, n1 = (int) mesh1->samples.size(), n2 = (int) mesh2->samples.size(); //n1=n2 in my case but this equality is not important for this function
	if (n1 < 4 || n2 < 4)
		cout << "WARNING: insufficient # samples for quadruples\n\n";
	for (int i1 = 0; i1 < n1; i1++)
		nb1 += (! mesh1->verts[ mesh1->samples[i1] ]->border);
	for (int i2 = 0; i2 < n2; i2++)
		nb2 += (! mesh2->verts[ mesh2->samples[i2] ]->border);
	if (nb1 < 4 || nb2 < 4 && mesh1->borderTreat == 2) //borderTreat=1 already called doubleBorderAreas() that doubled the areas so cannot switch to =0 in =1 mode
		cout << "WARNING: insufficient # non-border samples so switching to borderTreat = 0\n\n", mesh1->borderTreat = mesh2->borderTreat = 0; //cannot switch to borderTreat=1 'cos areas() already called
	//fread permutations from file to perms[][] mtrx
	FILE* fPtr;
	char fName[] = "permutations\\4!.txt"; //obtained using matlab's perms function, e.g., perms(0:3) for 4!.txt
	int factR = factorial(4);
	if (! (fPtr = fopen(fName, "r")))
	{
		cout << fName << " not found\n";
		exit(0);
	}
	//matlab perms(0:10) took 1 sec and fprinting the resulting 10! entries to the 10!.txt file took 2 mins
	int** perms = new int*[factR];
	for (int i = 0; i < factR; i++) //r! reorderings/permutations
	{
		perms[i] = new int[4];
		for (int j = 0; j < 4; j++)
			fscanf(fPtr, "%d\t", &perms[i][j]);
	}
	fclose(fPtr);

	//use permutations to perform the brute-force map computation
	vector< int > optimalMap;
	double minDistortion = INF;
	//for each quadruple in source mesh1
	for (int i1 = 0; i1 < n1; i1++)
		for (int j1 = i1 + 1; j1 < n1; j1++) //each pair
			for (int k1 = j1 + 1; k1 < n1; k1++) //each triplet
				for (int l1 = k1 + 1; l1 < n1; l1++) //each quadruple
			{
				//mesh1.samples[i1], mesh1.samples[j1], mesh1.samples[k1], mesh1.samples[l1] versus the mesh2 quadruple below
				vector< int > mesh1Samples;
				mesh1Samples.push_back(mesh1->samples[i1]), mesh1Samples.push_back(mesh1->samples[j1]), mesh1Samples.push_back(mesh1->samples[k1]), mesh1Samples.push_back(mesh1->samples[l1]);
				for (int i2 = 0; i2 < n2; i2++)
					for (int j2 = i2 + 1; j2 < n2; j2++) //each pair
						for (int k2 = j2 + 1; k2 < n2; k2++) //each triplet
							for (int l2 = k2 + 1; l2 < n2; l2++) //each quadruple
						{
							//compare the current mesh1 quadruple w/ each 4!=24 permutation of this current mesh2 quadruple: mesh2.samples[i2], mesh2.samples[j2], mesh2.samples[k2], mesh2.samples[l2]
							vector< int > mesh2Samples;
							mesh2Samples.push_back(mesh2->samples[i2]), mesh2Samples.push_back(mesh2->samples[j2]), mesh2Samples.push_back(mesh2->samples[k2]), mesh2Samples.push_back(mesh2->samples[l2]);
							for (int i = 0; i < factR; i++) //factR = 24 here for sure
							{
								vector< int > tmpMap; //2 consecutive entries (mesh1-mesh2) gives 1 match of the entire map/correspondence
								int idx1 = 0, idx2 = 0;
								for (int y = 0; y < 8; y++) //tmpMap will always have 4 matches in this function so 2*4 = 8 iterations (4 for mesh1 side 4 for mesh2 side)
									if (y % 2 == 0) //even idxs from the fixed mesh1 side
										tmpMap.push_back(mesh1Samples[idx1++]);
									else			//odd idxs for permutations from mesh2 side
									{
										tmpMap.push_back(mesh2Samples[ perms[i][idx2] ]);
										mesh1->verts[ mesh1Samples[idx1-1] ]->matchIdx = mesh2Samples[ perms[i][idx2++] ]; //mapQuality() needs matchIdx values (-1 'cos idx1 was already incremented)
									}
								if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
									mapQualityPatch(); //global distortion is updated
								else
									mapQualityGeo(); //global distortion is updated
								if (distortion < minDistortion)
									minDistortion = distortion, optimalMap = tmpMap;
								mesh1->verts[ mesh1Samples[0] ]->matchIdx = mesh1->verts[ mesh1Samples[1] ]->matchIdx = mesh1->verts[ mesh1Samples[2] ]->matchIdx = mesh1->verts[ mesh1Samples[3] ]->matchIdx = -1; //reset to unmatched for the next iteration
//if (distortion<6.36e-06+TINY && distortion>6.36e-06-TINY)for (int s = 0; s < 4; s++)cout << tmpMap[2*s] << " -test- " << mesh1->verts[ tmpMap[2*s] ]->matchIdx << endl;
							}
						}
			}
	//optimalMap to matchIdxs
	for (int i1 = 0; i1 < n1; i1++)
		mesh1->verts[ mesh1->samples[i1] ]->matchIdx = -1; //matchIdx have the very last configuration above so reset them to -1
	for (int s = 0; s < 4; s++) //optimalMap will always have 4 matches in this function
		mesh1->verts[ optimalMap[2*s] ]->matchIdx = optimalMap[2*s + 1], mesh2->verts[ optimalMap[2*s + 1] ]->matchIdx = optimalMap[2*s]; //dual operation
	//use the latest/optimal matchIdx values to assign individual v.distortion scores for printing only
	if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
		mapQualityPatch(); //global distortion is updated
	else
		mapQualityGeo(); //global distortion is updated
	for (int s = 0; s < 4; s++) //optimalMap will always have 4 matches in this function
	{
		cout << optimalMap[2*s] << " -bruteQuad- " << optimalMap[2*s + 1] << "\t" << mesh1->verts[ optimalMap[2*s] ]->distortion;
		if (mesh1->verts[ optimalMap[2*s] ]->patchGC == 0.0) //not in use so no print
			cout << endl;
		else
			 cout << "\t" << abs(mesh1->verts[ optimalMap[2*s] ]->patchGC - mesh2->verts[ optimalMap[2*s + 1] ]->patchGC) //difference b/w patch gc also printed after the v.distortion score
				  << "\t" << abs(mesh1->verts[ optimalMap[2*s] ]->area - mesh2->verts[ optimalMap[2*s + 1] ]->area) << "\n";//<< areaVecDiffs[s] << "\t" << areaVecDiffsRing[s] << endl; //area and areas/areasRing[] differences also printed
	}
	//recapture memory
	for (int i = 0; i < factR; i++)
		delete [] perms[i];
	delete [] perms;
	cout << ">>>>>>>> " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs for matching samples (bruteForceQuadruple distortion=" << dType << ": " << distortion << ")\n";	
	mapOK = true;
}

void Correspondence::bruteForceMapQuin()
{
	//same as bruteForceMap() except each C(n, 5) source quintuple is compared w/ each C(n, 5) target quintuple here (so output is a 5x5 bijective map)

	clock_t t = clock();
	cout << "Brute-force quintuple map computation..\n";
	int nb1 = 0, nb2 = 0, n1 = (int) mesh1->samples.size(), n2 = (int) mesh2->samples.size(); //n1=n2 in my case but this equality is not important for this function
	if (n1 < 5 || n2 < 5)
		cout << "WARNING: insufficient # samples for quintuples\n\n";
	for (int i1 = 0; i1 < n1; i1++)
		nb1 += (! mesh1->verts[ mesh1->samples[i1] ]->border);
	for (int i2 = 0; i2 < n2; i2++)
		nb2 += (! mesh2->verts[ mesh2->samples[i2] ]->border);
	if (nb1 < 5 || nb2 < 5 && mesh1->borderTreat == 2) //borderTreat=1 already called doubleBorderAreas() that doubled the areas so cannot switch to =0 in =1 mode
		cout << "WARNING: insufficient # non-border samples so switching to borderTreat = 0\n\n", mesh1->borderTreat = mesh2->borderTreat = 0; //cannot switch to borderTreat=1 'cos areas() already called
	//fread permutations from file to perms[][] mtrx
	FILE* fPtr;
	char fName[] = "permutations\\5!.txt"; //obtained using matlab's perms function, e.g., perms(0:3) for 4!.txt
	int factR = factorial(5);
	if (! (fPtr = fopen(fName, "r")))
	{
		cout << fName << " not found\n";
		exit(0);
	}
	//matlab perms(0:10) took 1 sec and fprinting the resulting 10! entries to the 10!.txt file took 2 mins
	int** perms = new int*[factR];
	for (int i = 0; i < factR; i++) //r! reorderings/permutations
	{
		perms[i] = new int[5];
		for (int j = 0; j < 5; j++)
			fscanf(fPtr, "%d\t", &perms[i][j]);
	}
	fclose(fPtr);

	//use permutations to perform the brute-force map computation
	vector< int > optimalMap;
	double minDistortion = INF;
	//for each quintuple in source mesh1
	for (int i1 = 0; i1 < n1; i1++)
		for (int j1 = i1 + 1; j1 < n1; j1++) //each pair
			for (int k1 = j1 + 1; k1 < n1; k1++) //each triplet
				for (int l1 = k1 + 1; l1 < n1; l1++) //each quadruple
					for (int m1 = l1 + 1; m1 < n1; m1++) //each quintuple
			{
				//mesh1.samples[i1], mesh1.samples[j1], mesh1.samples[k1], mesh1.samples[l1], mesh1.samples[m1] versus the mesh2 quintuple below
				vector< int > mesh1Samples;
				mesh1Samples.push_back(mesh1->samples[i1]), mesh1Samples.push_back(mesh1->samples[j1]), mesh1Samples.push_back(mesh1->samples[k1]), mesh1Samples.push_back(mesh1->samples[l1]), mesh1Samples.push_back(mesh1->samples[m1]);
				for (int i2 = 0; i2 < n2; i2++)
					for (int j2 = i2 + 1; j2 < n2; j2++) //each pair
						for (int k2 = j2 + 1; k2 < n2; k2++) //each triplet
							for (int l2 = k2 + 1; l2 < n2; l2++) //each quadruple
								for (int m2 = l2 + 1; m2 < n2; m2++) //each quintuple
						{
							//compare the current mesh1 quintuple w/ each 5!=120 permutation of this current mesh2 quintuple: mesh2.samples[i2], mesh2.samples[j2], mesh2.samples[k2], mesh2.samples[l2], mesh2.samples[m2]
							vector< int > mesh2Samples;
							mesh2Samples.push_back(mesh2->samples[i2]), mesh2Samples.push_back(mesh2->samples[j2]), mesh2Samples.push_back(mesh2->samples[k2]), mesh2Samples.push_back(mesh2->samples[l2]), mesh2Samples.push_back(mesh2->samples[m2]);
							for (int i = 0; i < factR; i++) //factR = 120 here for sure
							{
								vector< int > tmpMap; //2 consecutive entries (mesh1-mesh2) gives 1 match of the entire map/correspondence
								int idx1 = 0, idx2 = 0;
								for (int y = 0; y < 10; y++) //tmpMap will always have 5 matches in this function so 2*5 = 10 iterations (5 for mesh1 side 5 for mesh2 side)
									if (y % 2 == 0) //even idxs from the fixed mesh1 side
										tmpMap.push_back(mesh1Samples[idx1++]);
									else			//odd idxs for permutations from mesh2 side
									{
										tmpMap.push_back(mesh2Samples[ perms[i][idx2] ]);
										mesh1->verts[ mesh1Samples[idx1-1] ]->matchIdx = mesh2Samples[ perms[i][idx2++] ]; //mapQuality() needs matchIdx values (-1 'cos idx1 was already incremented)
									}
								if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
									mapQualityPatch(); //global distortion is updated
								else
									mapQualityGeo(); //global distortion is updated
								if (distortion < minDistortion)
									minDistortion = distortion, optimalMap = tmpMap;
								mesh1->verts[ mesh1Samples[0] ]->matchIdx = mesh1->verts[ mesh1Samples[1] ]->matchIdx = mesh1->verts[ mesh1Samples[2] ]->matchIdx = mesh1->verts[ mesh1Samples[3] ]->matchIdx = mesh1->verts[ mesh1Samples[4] ]->matchIdx = -1; //reset to unmatched for the next iter
							}
						}
			}
	//optimalMap to matchIdxs
	for (int i1 = 0; i1 < n1; i1++)
		mesh1->verts[ mesh1->samples[i1] ]->matchIdx = -1; //matchIdx have the very last configuration above so reset them to -1
	for (int s = 0; s < 5; s++) //optimalMap will always have 5 matches in this function
		mesh1->verts[ optimalMap[2*s] ]->matchIdx = optimalMap[2*s + 1], mesh2->verts[ optimalMap[2*s + 1] ]->matchIdx = optimalMap[2*s]; //dual operation
	//use the latest/optimal matchIdx values to assign individual v.distortion scores for printing only
	if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
		mapQualityPatch(); //global distortion is updated
	else
		mapQualityGeo(); //global distortion is updated
	for (int s = 0; s < 5; s++) //optimalMap will always have 5 matches in this function
	{
		cout << optimalMap[2*s] << " -bruteQuin- " << optimalMap[2*s + 1] << "\t" << mesh1->verts[ optimalMap[2*s] ]->distortion << "\t";
		if (mesh1->verts[ optimalMap[2*s] ]->patchGC == 0.0) //not in use so no print
			cout << endl;
		else
			 cout << abs(mesh1->verts[ optimalMap[2*s] ]->patchGC - mesh2->verts[ optimalMap[2*s + 1] ]->patchGC) << "\n"; //difference b/w patch gc also printed after the v.distortion score
	}
	//recapture memory
	for (int i = 0; i < factR; i++)
		delete [] perms[i];
	delete [] perms;
	cout << ">>>>>>>> " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs for matching samples (bruteForceQuintuple distortion=" << dType << ": " << distortion << ")\n";
	mapOK = true;
}

bool printOnce = true;
bool Correspondence::bruteForceMapQuinC2F()
{
	//same as bruteForceMapQuin() except 5-size (no C(n, 5) selection) source quintuple is compared w/ 5-size target quintuple here (so output is a 5x5 bijective map)

	if (mesh1->samplesC2F.size() != 5)
	{
		if (printOnce)
			cout << "WARNING: bad # c2f samples for quintuples; falling-back to qapMap() (there may be others; printing stopped)\n\n", printOnce = false; //not the end of the world so this is a light warning hence print once
		return false;
	}

	//use permutations to perform the brute-force map computation
	vector< int > optimalMap;
	double minDistortion = INF;
	int factR = factorial(mesh1->samplesC2F.size()), tmpSize = mesh1->samplesC2F.size() * 2;
	//mesh1.samplesC2F[0], mesh1.samplesC2F[1], mesh1.samplesC2F[2], mesh1.samplesC2F[3], mesh1.samplesC2F[4] versus the mesh2 quintuple below
	//compare the current mesh1 quintuple w/ each 5!=120 permutation of this current mesh2 quintuple: mesh2.samplesC2F[0], mesh2.samplesC2F[1], mesh2.samplesC2F[2], mesh2.samplesC2F[3], mesh2.samplesC2F[4]
	for (int i = 0; i < factR; i++)
	{
		for (size_t s = 0; s < mesh1->samples.size(); s++) //not mesh1->samplesC2F 'cos mapQualityPatch/Geo() goes over all samples[] not just samplesC2F
			mesh1->verts[ mesh1->samples[s] ]->matchIdx = -1; //set by y-loop below and then used by mapQualityPatch/Geo() below
		vector< int > tmpMap; //2 consecutive entries (mesh1-mesh2) gives 1 match of the entire map/correspondence
		int idx1 = 0, idx2 = 0;
		for (int y = 0; y < tmpSize; y++) //tmpMap will always have 5 matches in this function so 2*5 = 10 iterations (5 for mesh1 side 5 for mesh2 side)
			if (y % 2 == 0) //even idxs from the fixed mesh1 side
				tmpMap.push_back(mesh1->samplesC2F[idx1++]);
			else			//odd idxs for permutations from mesh2 side
			{
				tmpMap.push_back(mesh2->samplesC2F[ perms5[i][idx2] ]);
				mesh1->verts[ mesh1->samplesC2F[idx1-1] ]->matchIdx = mesh2->samplesC2F[ perms5[i][idx2++] ]; //mapQuality() needs matchIdx values (-1 'cos idx1 was already incremented)
			}
		if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
			mapQualityPatch(); //global distortion is updated
		else
			mapQualityGeo(); //global distortion is updated
		if (distortion < minDistortion)
			minDistortion = distortion, optimalMap = tmpMap;
	}
	//optimalMap to matchIdxs
	for (size_t i1 = 0; i1 < mesh1->samples.size(); i1++)
		mesh1->verts[ mesh1->samples[i1] ]->matchIdx = -1; //matchIdx have the very last configuration above so reset them to -1
	for (int s = 0; s < tmpSize / 2; s++) //optimalMap will always have 5 matches in this function
		mesh1->verts[ optimalMap[2*s] ]->matchIdx = optimalMap[2*s + 1], mesh2->verts[ optimalMap[2*s + 1] ]->matchIdx = optimalMap[2*s]; //dual operation
	nQAPs++; //counting bruteForce maps (naming sucks)
	return true;
}

void Correspondence::bruteForceMap(size_t m)
{
	//matches mesh1 and mesh2 samples via exhaustive/brute-force search that minimizes a distortion metric; n-m virtual samples are added to the partial mesh2 model

	clock_t t = clock();
	cout << "Brute-force map computation..\n";
	size_t n = mesh1->samples.size(); //force n samples match w/ m samples and n-m virtual samples (virtual matches are then discarded since virtual samples do not even belong to mesh2 surface)
	size_t nb1 = 0, nb2 = 0;
	for (size_t i1 = 0; i1 < n; i1++)
		nb1 += (! mesh1->verts[ mesh1->samples[i1] ]->border);
	for (size_t i2 = 0; i2 < m; i2++)
		nb2 += (! mesh2->verts[ mesh2->samples[i2] ]->border);
	if (nb1 < m || nb2 < m && mesh1->borderTreat == 2) //borderTreat=1 already called areas() that doubled the areas so cannot switch to =0 in =1 mode
		cout << "WARNING: insufficient # non-border samples so switching to borderTreat = 0\n\n", mesh1->borderTreat = mesh2->borderTreat = 0; //cannot switch to borderTreat=1 'cos areas() already called
	//fread permutations from file to perms[][] mtrx
	FILE* fPtr;
	char fName[250];
	int factR = factorial(n); //other permutation less than 10 are available in fName folder but here i just need 10! (r=nBruteSamples)
	sprintf(fName, "permutations\\%d!.txt", n); //obtained using matlab's perms function, e.g., perms(0:5) for 6!.txt
	if (! (fPtr = fopen(fName, "r")))
	{
		cout << fName << " not found\n";
		exit(0);
	}
	//matlab perms(0:10) took 1 sec and fprinting the resulting 10! entries to the 10!.txt file took 2 mins
	int** perms = new int*[factR], winnerPerm = -1;
	for (int i = 0; i < factR; i++) //r! reorderings/permutations
	{
		perms[i] = new int[n];
		for (size_t j = 0; j < n; j++)
			fscanf(fPtr, "%d\t", &perms[i][j]);
	}
	fclose(fPtr);
	cout << ">>>>>>>> " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs for freading permutations\n";

	//use permutations to perform the brute-force map computation
	vector< int > optimalMap;
	double minDistortion = INF;
	t = clock();//mapOK = false; //will be true if at least 1 permutation gives a INF-free distortion
	for (int i = 0; i < factR; i++)
	{
		//this exhaustive loop takes 3 secs for r=10 (10!=3628800 permutations fread above in 5 additional secs so 8 secs in total), which is still slower than my old Genetic Algorithm paper which takes 0.5 secs to compute the same optimal 10x10 map
		//it takes 0.003 secs r=6, 0.03 secs r=8 (of course it depends on the distortion computation time in mapQualityGeo() too)
		vector< int > tmpMap; //2 consecutive entries (mesh1-mesh2) gives 1 match of the entire map/correspondence
		int idx1 = 0, idx2 = 0;
		for (size_t y = 0; y < 2*n; y++)
			if (y % 2 == 0) //even idxs for permutations from mesh1 side
				tmpMap.push_back(mesh1->samples[ perms[i][idx1] ]);
			else			//odd idxs from the fixed mesh2 side w/ virtual samples (n-m indices not covered by the current permutation are assumed to matched to virtual mesh2 samples)
			{
				if (idx2 < (int) m) //fixed 0,1,2,..,m part
				{
					tmpMap.push_back(mesh2->samples[idx2]);
					mesh1->verts[ mesh1->samples[ perms[i][idx1] ] ]->matchIdx = mesh2->samples[idx2]; //mapQuality() needs matchIdx values
				}
				else          //virtual sample part
				{
					tmpMap.push_back(-2);
					mesh1->verts[ mesh1->samples[ perms[i][idx1] ] ]->matchIdx = -2; //special flag for virtual samples
				}
				idx1++, idx2++;
			}
		mapQualityVirtual(); //global distortion is updated
		if (distortion < minDistortion)
		{
			minDistortion = distortion;
			optimalMap = tmpMap;
			winnerPerm = i;
		}
		for (size_t i1 = 0; i1 < n; i1++)
			mesh1->verts[ mesh1->samples[i1] ]->matchIdx = -1; //matchIdx have the very last configuration above so reset them to -1 (optimalMap to matchIdxs below also relies on this resetting)
	}
	if (mapOK) //optimalMap is not empty
	{
		//optimalMap to matchIdxs	
		for (size_t s = 0; s < n; s++)
			if (optimalMap[2*s + 1] != -2) //mesh2 side is not a virtual sample so match is valid
				mesh1->verts[ optimalMap[2*s] ]->matchIdx = optimalMap[2*s + 1], mesh2->verts[ optimalMap[2*s + 1] ]->matchIdx = optimalMap[2*s]; //dual operation
		//use the latest/optimal matchIdx values to assign individual v.distortion scores for printing only
		mapQualityVirtual(); //global distortion is updated
		for (size_t s = 0; s < n; s++)
			if (optimalMap[2*s + 1] != -2) //mesh2 side is not a virtual sample so match is valid
			{
				cout << optimalMap[2*s] << " -brute- " << optimalMap[2*s + 1] << "\t" << mesh1->verts[ optimalMap[2*s] ]->distortion << "\t";
				if (mesh1->verts[ optimalMap[2*s] ]->patchGC == 0.0) //not in use so no print
					cout << endl;
				else
					 cout << abs(mesh1->verts[ optimalMap[2*s] ]->patchGC - mesh2->verts[ optimalMap[2*s + 1] ]->patchGC) //difference b/w patch gc also printed after the v.distortion score
					 << "\t" << abs(mesh1->verts[ optimalMap[2*s] ]->area - mesh2->verts[ optimalMap[2*s + 1] ]->area) << "\n";//<< areaVecDiffs[s] << endl; //area and areas[] differences also printed
			}
		cout << ">>>>>>>> " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs for matching samples (bruteForce virtualSamples perm" << winnerPerm << ", distortion=" << dType << ": " << distortion << ")\n";
	} //else cout << "WARNING: put more samples on mesh2 or increase maxGeoDiff threshold (increasing automatically in main())\n\n";
	//recapture memory
	for (int i = 0; i < factR; i++)
		delete [] perms[i];
	delete [] perms;	
}

void Correspondence::denseMap(int mapType, bool fileFPS, bool trim)
{
	//computes dense map b/w dense samples on source and target ROIs
	
	nQAPs = 0; //to see how long each qapMap() takes: nSamplesPerSubregion=5/10/15/16/17/18/20 --> 0.023/0.4/5.3/10.2/14.5/16.7/65.4 secs; brute-force matching of 5/10/15/16/17/18/20, e.g., 5x5 map, --> 0.0026/209/intractable/intr/intr/intr/intr secs (nPathSamples=nRingAreas=5 case which affects objective complexity); 37x37 qap in matchUnmatcheds() took 3017secs

	if (multipleC2Fs && nDenseSamples > 60) //do multiple c2f w/ new disjoint sample sets each time, e.g., finish c2f matching of teh first 30 of 250; then do the next disjoint 30 and so on
	{
		denseMapLargeN(30);
		return;
	}

	vector< int > coarseMapFull; //coarseMap[] is the trimmed version of this original coarseMapFull formed by main().qapMap() matches (used by ROI() and distortion averaging only)
	if (fileFPS) //rare case (for fair comparison only)
	{
		int v1 = 31674, v2 = 31639; //can also be found automatically by calling qapMap() from main() but that would have complicated code/program too much for a comparison action; setting coarseMap is not that critical anyway, i.e., this block can safely be erased
		mesh1->coarseMap.push_back(v1), mesh1->coarseMap.push_back(v2), mesh1->verts[v1]->matchIdx = v2, mesh2->verts[v2]->matchIdx = v1, mesh1->dijkstraShortestPaths(v1), mesh2->dijkstraShortestPaths(v2), mesh1->verts[v1]->coarseMapMatch = mesh2->verts[v2]->coarseMapMatch = true;
		v1 = 51914, v2 = 51917, mesh1->coarseMap.push_back(v1), mesh1->coarseMap.push_back(v2), mesh1->verts[v1]->matchIdx = v2, mesh2->verts[v2]->matchIdx = v1, mesh1->dijkstraShortestPaths(v1), mesh2->dijkstraShortestPaths(v2), mesh1->verts[v1]->coarseMapMatch = mesh2->verts[v2]->coarseMapMatch = true;
		v1 = 47824, v2 = 47699, mesh1->coarseMap.push_back(v1), mesh1->coarseMap.push_back(v2), mesh1->verts[v1]->matchIdx = v2, mesh2->verts[v2]->matchIdx = v1, mesh1->dijkstraShortestPaths(v1), mesh2->dijkstraShortestPaths(v2), mesh1->verts[v1]->coarseMapMatch = mesh2->verts[v2]->coarseMapMatch = true;
		mesh1->maxGeoDist = mesh1->verts[31674]->geodesic[51914], mesh2->maxGeoDist = mesh2->verts[31639]->geodesic[51917]; //max geo on complete humans is from hand to foot
	}
	else //common case of usual program execution
	{
		for (size_t s = 0; s < mesh1->samples.size(); s++)
			if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
				coarseMapFull.push_back(mesh1->samples[s]), coarseMapFull.push_back(mesh1->verts[ mesh1->samples[s] ]->matchIdx); //each 2 consecutive items make 1 match (even indices belong to mesh1, odd indices mesh2)
		//make coarseMap[] ready for blossom costType > 1 (also for pathSamplingCoarse() and fillQdense() and maybe other things)
		for (size_t s = 0; s < mesh1->samples.size(); s++)
			if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1 && mesh1->verts[ mesh1->samples[s] ]->distortion < 1.1*distortion/(coarseMapFull.size()/2)) //2nd condition to keep coarseMap very small (efficiency) and accurate (accuracy in, e.g., pathSamplingCoarse() and tmp tests) by admitting distortions below the 1.1avg distortion (1.1 'cos some matches must be above average but they may still be very good)
				mesh1->coarseMap.push_back(mesh1->samples[s]), mesh1->coarseMap.push_back(mesh1->verts[ mesh1->samples[s] ]->matchIdx), mesh1->verts[ mesh1->samples[s] ]->coarseMapMatch = mesh2->verts[ mesh1->verts[ mesh1->samples[s] ]->matchIdx ]->coarseMapMatch = true; //even and odd idxs mesh1 and mesh2 side, respectively
			else if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1) //high distortion, i.e., in coarseMapFull[] but not in coarseMap, so reset matchIdx
				mesh2->verts[ mesh1->verts[ mesh1->samples[s] ]->matchIdx ]->matchIdx = -1, mesh1->verts[ mesh1->samples[s] ]->matchIdx = -1; //ROI() uses coarseMapFull[] but not need these matchIdx info so no problem
	}
	mesh2->coarseMap = mesh1->coarseMap; //coarseMap[] will not be updated but if it was updated then a change in mesh2.coarseMap would not affect mesh1.corresp (verified)
	int m = mesh2->samples.size(); //n mesh1 samples are matched w/ m mesh2 samples (m <= n) during main.qapMap() so remember m to adjust areaWeight as follows:
	//3rd dimension of Q may change, i.e., number of ring areas on a sample or number of patch areas between 2 samples, but it does not affect the weighting b/w unary and pairwise terms (each term gets the L2 difference (of areasRing[] for instance) w/ new 3rd dimension so no problem); the change on mesh2 sample size,
	//however, affect the weighting (areaWeight) since 1 unary + m-1 pairwise are added from each Q row to get the overall distortion; new m is nSamplesPerSubregion=5 in C2FQAP mode and nDenseSamples=30 in QAP mode; to equalize unary term's effect (see last line in setAreaWeight()) multiply areaWeight by the newM/m ratio
//nDenseSamples = 100;//250;//150;//50;//30; //nDenseSamples=30&nPathSamples=4 takes 62secs, =50&4 takes 140secs so reduce 4->2 below if nDenseSamples > 30, =50&1 takes 6secs, =50&2 takes 63secs
	mesh1->eucPatching = (mesh1->verts.size() > 20000);//19000); //geodesic-based patching, despite the early break in patchAndDescriptors(), is significantly slow (13secs on 52K-vertex mesh w/ nSamples=nPathSamples=10 (eucPatching takes just 0.3secs), 20secs nSamples=30,nPathSamples=3) so use euclidean distance during patching for a 40x speedup
	mesh2->eucPatching = (mesh2->verts.size() > 20000);//19000); //geodesic-based patching, despite the early break in patchAndDescriptors(), is significantly slow (13secs on 52K-vertex mesh w/ nSamples=nPathSamples=10 (eucPatching takes just 0.3secs), 20secs nSamples=30,nPathSamples=3) so use euclidean distance during patching for a 40x speedup
	if (mesh1->eucPatching || mesh2->eucPatching)
		cout << "WARNING: path patches by euclidean distance for efficiency; geodesic is more accurate but slow\n\n";
	//region of interests based on the correspondence information
	double r2 = mesh2->ROI(coarseMapFull), r1 = mesh1->ROI(coarseMapFull, r2); //mesh2 is the partial model whose roiArea r2 gives an upper bound for the roiArea of the complete mesh1 (for efficiency only as ROI() is costly for dense, e.g., 52K-vertex, models); no problem for complete matching scenario
	if (max(r1/r2, r2/r1) > 2)
		cout << "WARNING: source and target ROI areas too different, implying bad map (recompute the map)\n\n";	
	bool curvFPS = true;//nDenseSamples <= 100
	if (fileFPS) //rare case (for fair comparison only)
		mesh1->FPSfile(nDenseSamples), mesh2->FPSfile(nDenseSamples), mesh1->radius = mesh2->radius = max(mesh1->avgEdgeLen, mesh2->avgEdgeLen) * 6.0 * 2.5, mesh1->nRingAreas = mesh2->nRingAreas = 5; //from main()
	else //common case of usual program execution
		mesh1->FPSroi(nDenseSamples, curvFPS), mesh2->FPSroi(nDenseSamples, curvFPS); //put many FPS samples on ROIs (coarseMap[] samples survive and the rest is erased by FPSroi().clearSamples())
//return;//see rois only (compatibleColoring(false,true)) or dense sampling only

	//mesh1->radius = mesh2->radius = max(mesh1->avgEdgeLen, mesh2->avgEdgeLen) * 2.0; //smaller radius than the one use during coarse map computation (not in use 'cos then patches would be too small; 20secs-->57secs on complete model when not in use)
//	qapSolver = QUCOOPMATCH;
//qapSolver = QUCOOP;
//qapSolver = QMATCH;
//qapSolver = QGMATCH;
	sprintf(qapSolverStr, "%s", (qapSolver == QUCOOP ? "qucoop" : "qucoopmatch"));
	if (qapSolver == QMATCH)
		sprintf(qapSolverStr, "qmatch");

//	dissimType = (nDenseSamples <= 50 /*|| (mesh1->eucPatching && mesh2->eucPatching)*/ ? dissimType : PAIRWISENOPATCH), dissimTypeString(); //PAIRWISENOPATCH (pure geodesic lengths in use) or POINTWISEPAIRWISE (pathing/patching in use: nPathSamples+1 samples/areas between samples (for off-diagonal Q entries) and ring areas (for diagonal Q entries)) or POINTWISE or PAIRWISE
//	dissimType = POINTWISEPAIRWISE, dissimTypeString();
//	dissimType = PAIRWISENOPATCH, dissimTypeString();
//	dissimType = PAIRWISE, dissimTypeString();

	size_t cSize = mesh1->coarseMap.size() / 2; //number of matches in coarseMap[] (2 consecutive entries make 1 match)
	int nPathSamples = (dissimType == POINTWISE || dissimType == PAIRWISENOPATCH ? -1 : 4); //paths from dense samples to coarseMap will mostly be short so use small radius and small nPathSamples (4 (makes total 4+1 samples); too slow if >4)
//nPathSamples = ((nDenseSamples > 50 || mesh1->verts.size() > 20000) && nPathSamples != -1 ? 2/*1*/ : nPathSamples); //disabled 'cos more samples are better as i go to distant coarseMap samples in extraSamples mode; nPathSamples=2: 21secs vs. 4: 33secs
	int nPathSamplesExtra = cSize * nPathSamples + nPathSamples+1; //3*4 (towards 3 coarseMap samples (not counting the common samples[j] end)) + 5 from samples[i] to samples[j] (actual path w/o the extra sampling)
	mesh1->extraSamples = mesh2->extraSamples = (dissimType == PAIRWISE || dissimType == POINTWISEPAIRWISE ? mesh1->extraSamples : false); //extraSamples is ineffective (false) if no pathing/patching in use, i.e., nPathSamples = -1 (regardless of its constructor initialization)	
	bool displayCoarse = true; //ineffective if patched path visualization is active in compatibleColoring()
	char str[10], orgDtype[250]; sprintf(orgDtype, dType); clock_t t; //same as strcpy(orgDtype, dType) or sprintf(orgDtype, "%s", dType);
	if (mapType == QAP || mapType == C2FQAP) //quadratic assignment problem via qucoop (one-shot or coarse-to-fine approach)
	{
		//might be changed to 1 in main(), e.g., david0-10_partial3, so return to default 2 to be safe
		mesh1->nonuniformPathSamplingHeuristic = mesh2->nonuniformPathSamplingHeuristic = 2;//0; //nonuniform sampling (=1 or 2) is 2x slower which may be redundant in dense sampling (not really 'cos extraSamples makes me go long paths towards coarseMap); 0: uniform (2x faster),
		//1: for each patch p (except endpoint patches), replace its centering sample w/ the max-curvature locMaxGC path vertex inside p; 2: replace the centering sample with the max-curvature locMaxGC patch (not path) vertex which makes us go off the path but not that off so sill ok
		//pairwise structure
		if (dissimType == PAIRWISE || dissimType == POINTWISEPAIRWISE) //path sampling requested (nPathSamples!=-1 for sure)
		{
			if (mesh1->extraSamples)
				mesh1->pathSamplingCoarse(nPathSamples), mesh2->pathSamplingCoarse(nPathSamples); //make samplesToCoarse[] ready for the pathSampling() calls below
			mesh1->pathSampling(nPathSamples), mesh2->pathSampling(nPathSamples); //unnormalized euc/geo distance vs. unnormalized radius so geo must be unnormalized (no unitGeoScaling()) prior to pathSampling/pathSamplingCoarse() calls
		}
		//pointwise descriptor
		if (dissimType == POINTWISE || dissimType == POINTWISEPAIRWISE) //pointwise descriptor is set by areas() so dissimType must include pointwise (to save tiny time)
		{
			mesh1->areas(), mesh2->areas(); //radius has not changed since coarseMap[] computation but new dense samples arrived so call this function to fill their areasRing[] and areas[] entries
			if (mesh1->borderTreat == 1)
				mesh1->doubleBorderAreas(), mesh2->doubleBorderAreas(); //avgAreas() called inside doubleBorderAreas()
			else
				mesh1->avgAreas(), mesh2->avgAreas(); //set avgPatch/RingArea values used by areaWeight
		}
//cout << nDenseSamples << " == new" << mesh2->samples.size() << " must be equal\n";//verified
		//newM/m areaWeight adjustments (no barrier adjustment 'cos original main().barrier is based on avgPairwise value which does not change as i use the same radius as main() here; note that barrier is added to EACH pairwise term, via if (pairwise)Q(u,v)+=barrier*.., so it's quite effective)
		if (mapType == C2FQAP)
			areaWeight *= ((double) nSamplesPerSubregion / m), printQAP = false, cout << "C2FQAP map computation via " << qapSolverStr << "..\n"; //no printing during C2FQAP as qapMap() will be called many times in that mode (newM/m adjustment explained above is done here)
		else
			areaWeight *= ((double) nDenseSamples / m); //newM/m adjustment explained above
//areaWeight = 0; //no unary term at all
//barrier = 0; //no barrier at all
#ifndef GEODIFFTHRESHOLDING
		cout << "PARAMS (dense): radius = " << mesh1->radius << ", nRingAreas = " << mesh1->nRingAreas << ", nPathSamples = " << nPathSamples+1 << ", unaryWeight = " << areaWeight << ", barrier = " << (barrier > 0 ? barrier : 0) << (mesh1->weightByGC ? ", wByGC" : "") << ", areasRing, qapEnergy: " << dissimType << ", pathSampler: " << mesh1->nonuniformPathSamplingHeuristic << "\n"; //barrier n/a if negative
#else
		barrier = 0, cout << "PARAMS (dense): radius = " << mesh1->radius << ", nRingAreas = " << mesh1->nRingAreas << ", nPathSamples = " << nPathSamples+1 << ", unaryWeight = " << areaWeight << ", maxGeoDiff = " << maxGeoDiff << (mesh1->weightByGC ? ", wByGC" : "") << ", areasRing, qapEnergy: " << dissimType << ", pathSampler: " << mesh1->nonuniformPathSamplingHeuristic << "\n"; //non-overlapping areasRing[] alternative is overlapping is areas[] for the unary term
#endif
		//although new roi-based maxGeoDist has not changed since FPSroi() above, i need to defer unitGeoScaling() here due to unnormalized euc/geo distance vs. unnormalized radius part above; unitGeoScaling() mandatory due to fillQdense() and also FPSpatch()
		g = mesh1->maxGeoDist/*max(mesh1->maxGeoDist, mesh2->maxGeoDist)*/, mesh1->unitGeoScaling(g), mesh2->unitGeoScaling(g); //mesh1.maxGeoDist = ~mesh2.maxGeoDist expected since ROIs are compatible so max is ~redundant but anyway; actually pfaust partial models may give a bigger maxGeoDist due to lack of cropped regions/paths so remove max()
		cSize = (nPathSamples == -1 ? 0 : cSize); //disable cSize effect to be able to send special flag -1 for Q2D usage (in no patch mode of POINTWISE or PAIRWISENOPATCH) (nPathSamples!=-1 for sure for the else extraSamples part below)
		if (mapType == QAP) //coarse-to-fine qap C2FQAP will work slightly differently so don't do it here
		{
//			double d = max(mesh1->avgRingArea, mesh2->avgRingArea), awObserve = (d > 0 ? max(mesh1->avgPatchArea, mesh2->avgPatchArea) / d : 1.0); not needed 'cos radius has not changed since main() so the same main() weights are valid except newM/m adjustment explained above
//			areaWeight = (dissimType == POINTWISEPAIRWISE ? main().setAreaWeight(mesh1, mesh2, avgPairwise) : awObserve); //automatic weighting that makes average unary and pairwise terms in Q (and in mapQualityPatch/Virtual()) approximately the same (ring/patchAreas affected (if radius changes) but maxGeoDiff is not affected by dense sampling (/=g below in LAP is for normalization only))
			t = clock();
			if (! mesh1->extraSamples) //nPathSamples=-1: means special flag for Q2Ddense usage (diag: lengthsToCoarseMap off-diag: lengthsBetweenSamples); nPathSamples>0: means diag: ringAreas off-diag: patchAreasBetweenSamples (alternative not in use: diag: patchAreasToCoarseMap off-diag: patchAreasBetweenSamples)
				mesh1->fillQdense(nPathSamples + 1, areaWeight), mesh2->fillQdense(nPathSamples + 1, areaWeight), //QAP needs the Q matrix
				qapMap(nPathSamples + 2 + cSize, true); //nPathSamples=-1 sends 1 which is a special flag for Q2D usage; otherwise nPathSamples + 1 gives # path samples and another + 1 gives access to the additional dimension storing pure geodesic length info (for maxGeoDiff filtering or geoDiff-based barrier) and another + cSize for lengths to coarseMap entries
			else //nPathSamples!=-1 for sure 'cos nPathSamples=-1 means no path sampling requested, i.e., POINTWISE or PAIRWISENOPATCH, and in that case extraSamples cannot be true
				mesh1->fillQdense(nPathSamplesExtra, areaWeight), mesh2->fillQdense(nPathSamplesExtra, areaWeight), //+1 added during declaration so no nPathSamplesExtra + 1 here
				qapMap(nPathSamplesExtra + 1 + cSize, true); //+ 1 gives access to the additional dimension storing pure geodesic length info (for maxGeoDiff filtering or geoDiff-based barrier) and another + cSize for lengths to coarseMap entries
		}			
	}
	else if (mapType == LAP) //linear assignment problem via blossom
	{
		cout << "disabled to simplify compilation; not in use anyway as LAP is much less robust\n\n";
/*		mesh1->eucPatching = (mesh1->verts.size() > 20000);//19000); //geodesic-based patching, despite the early break in patchAndDescriptors(), is significantly slow (13secs on 52K-vertex mesh w/ nSamples=nPathSamples=10 (eucPatching takes just 0.3secs), 20secs nSamples=30,nPathSamples=3) so use euclidean distance during patching for a 40x speedup
		mesh2->eucPatching = (mesh2->verts.size() > 20000);//19000); //geodesic-based patching, despite the early break in patchAndDescriptors(), is significantly slow (13secs on 52K-vertex mesh w/ nSamples=nPathSamples=10 (eucPatching takes just 0.3secs), 20secs nSamples=30,nPathSamples=3) so use euclidean distance during patching for a 40x speedup
		if (mesh1->eucPatching || mesh2->eucPatching)
			cout << "WARNING: path patches by euclidean distance for efficiency; geodesic is more accurate but slow\n\n";
		int bloCostType = POINTWISEGEO, //POINTWISEAREAandRING, POINTWISEAREA, POINTWISERING
			nPathSamplesLAP = 4; //same nPathSamples as fillQdense() is not mandatory (but logical) 'cos qap's fillQdense() uses paths b/w samples but this lap case uses paths towards coarseMap endpoints
		//solve linear assignment problem (lap) via the minimum-weight perfect matching provided by the blossom algorithm (a faster hungarian alternative); fillQdense() for np-hard quadratic assignment problem (qap for another project)
		if (bloCostType == POINTWISEAREA || bloCostType == POINTWISEAREAandRING) //patch areas over paths towards coarseMap needed (patches computed by fillQdense().pathSampling() will be cleared inside pathSamplingCoarse())
			mesh1->pathSamplingCoarse(nPathSamplesLAP, true), mesh2->pathSamplingCoarse(nPathSamplesLAP, true);
		if (bloCostType == POINTWISERING || bloCostType == POINTWISEAREAandRING) //area rings from around each sample needed
			mesh1->areas(), mesh2->areas();
		if (mesh1->borderTreat == 1)
			mesh1->doubleBorderAreas(), mesh2->doubleBorderAreas(); //avgAreas() called inside doubleBorderAreas()
		else
			mesh1->avgAreas(), mesh2->avgAreas();
#ifndef GEODIFFTHRESHOLDING
		//although new roi-based maxGeoDist has not changed since FPSroi() above, i need to defer unitGeoScaling() here due to the pathSamplingCoarse() which has unnormalized euc/geo distance vs. unnormalized radius inside patchAndDescriptors() that it calls
		g = mesh1->maxGeoDist, mesh1->unitGeoScaling(g), mesh2->unitGeoScaling(g); //mesh1.maxGeoDist = mesh2.maxGeoDist expected since ROIs are compatible so max is ~redundant but anyway; actually pfaust partial models may give a bigger maxGeoDist due to lack of cropped regions/paths so remove max()
		if (! noBarrier) //no avgPairwise*100 here 'cos removed setAreaWeight2() used pathSamples b/w all dense samples but here i just have pathSamples from coarseMap samples to dense samples; so always barrier = g*10
			barrier = g * 10.0; //make it negative to use geoDiff as a multiplier of pairwise term, not as a barrier term; make it 0 to completely disable geoDiff-based barrier stuff
		maxGeoDiff /= g; //normally maxGeoDiff'll not be used in this GEODIFFTHRESHOLDING-undefined mode so no update was done in constructor but here update is needed 'cos blossom.POINTWISEGEO uses maxGeoDiff filtering (blossom is linear assignment, not quadratic)
#endif
		clock_t t2 = clock();
		(new BloMatch(this, bloCostType, mesh1->coarseMap))->blossomMap(displayCoarse);
		sprintf(str, " -blo- "), sprintf(dType, "puregeo"), mapQualityGeo(), printMap(t2, str); //mapQualityGeo() and hence puregeo distortions will be printed so update dType
disabled to simplify compilation; not in use anyway as LAP is much less robust*/
	}
	if (mapType == C2FQAP) //quadratic assignment problem via qucoop (coarse-to-fine approach)
	{
//qapSolver == QUCOOP or QUCOOPMATCH, matchUnmatcheds(nPathSamples, nPathSamplesExtra); return;//see the pure qapSolver's accuracy and timing, i.e., no C2F heuristic (coarseMap[] (<1.1distortion matches above) are excluded in this function hence expect 27x27 instead of exact 30x30)
		//each matched sample (initially only coarseMap[] samples are matched and hence inside matches1) defines a patch which is resampled via FPSpatch() and matched via qapMap()
		vector< int > matches1, matches2; //nMaxLevels many levels will be used where the number of matches increases as levels increase
//nMaxLevels = 10;//2;//1; //nSamplesPerSubregion = 10;
		double minDisto = INF, minDistoDensest = INF, nMaxMatches = -INF, nMatches, radius = mesh1->maxGeoDist/*max(mesh1->maxGeoDist, mesh2->maxGeoDist)*/ / 2.0; //initial radius here will be halved at each level iteration; mesh1 is the complete model whose maxGeoDist is guaranteed to be 1 after unitGeoScaling()
		vector< int > optimalMap, densestOptMap; //map corresponding to the minimum-distortion map or densest minimum-distortion map through c2f levels
		bool patchCostDistortion = true, //false; //best match update based on patch-based pointwisepairwise (true and also dissimType != PAIRWISENOPATCH) or puregeo (false) distortions
			 printOnce1 = true, printOnce2 = true, backToOptimal = true, doTwice = false, noTrimCoarseMapMatches = true; //matches in the initial coarseMap[], i.e., sparse/coarse matches from main().qapMap() cannot be trimmed, i.e., they may change/refine but never be -1
		cout << "c2f w/ initial c2fRadius: " << radius << ", nSamplesPerSubregion: " << nSamplesPerSubregion << ", " << "back2optimal: " << backToOptimal << ", doTwice: " << doTwice << endl, t = clock();
		//evaluate the matches produced by the very first main().qapMap() using the new parameters shown in PARAMS print above and/or new coarseMap that is trimmed above; v.distortion ready now for level=0 case below (other levels done naturally w/ the mapQualityPatch/Geo() calls inside the loop)
		if (patchCostDistortion && dissimType != PAIRWISENOPATCH) //pathing/patching in use
			nMatches = mapQualityPatch(); //global distortion is updated
		else
			nMatches = mapQualityGeo(); //global distortion is updated (redundant 'cos PARAMS update above does not affect puregeo distortion in mapQualityGeo() that does not use coarseMap either)
		for (size_t s = 0; s < mesh1->samples.size(); s++)
		{
			if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
				mesh1->verts[ mesh1->samples[s] ]->distortion /= nMatches, //distortion normalization as the current map size is coarseMap.size() which is not necessarily compatible w/ the nSamplesPerSubregion map size below (normalize those distortions too via /=nSamplesPerSubregion below)
				matches1.push_back(mesh1->samples[s]), matches2.push_back(mesh1->verts[ mesh1->samples[s] ]->matchIdx); //matches1[] is based on coarseMap[] now (it'll expand below, possibly changing the matches of the previous coarse levels)
			mesh1->verts[ mesh1->samples[s] ]->geoToRef = mesh2->verts[ mesh2->samples[s] ]->geoToRef = -1.0; //to detect the untouched samples during nMaxLevels processing below; hopefully all samples are touched/processed/matched; trimming may leave some samples unmatched which is OK but they must be touched at least once during the nMaxLevels processing; -1 means untouched (for printing, FPSpatch() and distortion restoration)
		}
		for (int level = 0; level < nMaxLevels; level++)
		{
			for (size_t s = 0; s < mesh1->samples.size(); s++)
			{
				mesh1->verts[ mesh1->samples[s] ]->binID = -1, mesh1->verts[ mesh1->samples[s] ]->tmp = INF; //binID is the best matchIdx throughout this level's iterations (naming sucks) and tmp is that match's distortion (naming sucks)
				if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1) //mesh1->samples[s] already has a distortion from coarseMap by qapMap() (level=0 case) so use that distortion as tmp to naturally protect the good coarseMap matches (different from noTrimCoarseMapMatches which protects them from being trimmed; this line protects them from being overwritten inside samplesC2F loop below)
					mesh1->verts[ mesh1->samples[s] ]->tmp = mesh1->verts[ mesh1->samples[s] ]->distortion;//cout << mesh1->samples[s] << " -coarseMap- " << mesh1->verts[ mesh1->samples[s] ]->matchIdx << "\t" << mesh1->verts[ mesh1->samples[s] ]->distortion << dType << "normalized\n";
			}
			vector< int > rd; //verts whose distortions will be restored to the saved normalized distortions in geoToRef (naming sucks)
			int qual = 0; //most recent quality metric (mapQualityGeo/Patch()) for printing only
			for (size_t i = 0; i < matches1.size(); i++) //same as matches2.size()
			{
//if (i > 1) continue; //do matches1[0] and [1] only
//if (i == 0) continue;
//cout << "patch pair based on " << matches1[i] << " - " << matches2[i] << " match\n";
//if (level==1 && i!=4)continue;//check a particular 5x5 qap map
//if (level==1)cout << matches1.size() << " " << matches1[i] << " " << matches2[i] << "sssssss\n";				
				curvFPS = true;//level < 2; level2+ uses pure FPS (curvFPS=false) to increase the chances of touching all samples, e.g., samples w/ small gc'll never be touched (cannot win the pool competition) so give up on curvFPS at some point)
				//do each radius twice, first w/ typical FPSpatch and second w/ modified FPSpatch that selects untouched (geoToRef=-1, naming sucks) samples (to prevent many unmatched matches in the end due to the untouched samples)
				c2fHelper(matches1[i], matches2[i], radius, curvFPS, nPathSamples, nPathSamplesExtra, cSize, patchCostDistortion, qual, rd, false); //first
				if (doTwice) //second processing is optional (makes sense but did not increase the number of matches as much as i expected)
					c2fHelper(matches1[i], matches2[i], radius, curvFPS, nPathSamples, nPathSamplesExtra, cSize, patchCostDistortion, qual, rd, true); //second (did twice)
//break;//see the very first 5x5qap map (make nMaxLevels=1 above and disable matchUnmatcheds() below)
			} //each patch is processed and the winner binID (match) and its tmp (distortion) is saved for this current level (namings suck)
cout << matches1.size() << " --> ";
			vector< int > matches1Cpy = matches1, matches2Cpy = matches2; //keep the copies before clearing them as i'll use the matches here to rematch the ones that are unmatched via matchIdx=-1 below
			//refill matches1 for the next level iteration using the binID values that encode the matches up to this level (naming sucks); set the matchIdx values as well in case this is the last level iteration (other iterations'll reset matchIdx=-1 above)
			matches1.clear(), matches2.clear();
			if (mesh1->samples.size() != mesh2->samples.size())
				cout << "WARNING: FPSroi() must have created same-size mesh1/2 samples\n\n";
			for (size_t s = 0; s < mesh2->samples.size(); s++)
				mesh2->verts[ mesh2->samples[s] ]->matchIdx = -1, mesh2->verts[ mesh2->samples[s] ]->distortion = INF; //mesh2 side must also be reset 'cos 2+ (manyTo1) matches on the same mesh2 sample implies incompatible sampling and hence cancels both 2+ matches (or keep the best of 2+ matches); distortion is for CANCELONMANYTO1 undefined mode only
			for (size_t s = 0; s < mesh1->samples.size(); s++) //same as mesh2->samples.size()
			{
				mesh1->verts[ mesh1->samples[s] ]->matchIdx = -1;
				int v2 = mesh1->verts[ mesh1->samples[s] ]->binID;
				if (v2 != -1)
				{
#ifdef CANCELONMANYTO1 //cancel 2+ matches targeting the same mesh2 sample
					if (mesh2->verts[v2]->matchIdx == -1) //ensure no manyTo1 matches on mesh2 vertices
						matches1.push_back(mesh1->samples[s]), matches2.push_back(v2), mesh1->verts[ mesh1->samples[s] ]->matchIdx = v2, mesh2->verts[ mesh1->verts[ mesh1->samples[s] ]->matchIdx ]->matchIdx = mesh1->samples[s]; //dual operation (v2.distortion assignment is redundant here due to CANCELONMANYTO1 mode)
					else //manyTo1 would have happened so cancel this match (by not entering if above) and also cancel the existing match (by entering this else below)
					{
//if (mesh1->samples[s] == 39450 || v2 == 5234){
//int v1 = mesh2->verts[v2]->matchIdx; cout << "CANCELLLLLL: " << mesh1->samples[s] << " - " << v2 << " w/ " << mesh1->verts[ mesh1->samples[s] ]->tmp << " == " << mesh1->verts[ mesh1->samples[s] ]->distortion 
//										  << "\t\t" << v1 << " - " << v2 << " w/ disto: " << mesh1->verts[v1]->tmp << " == " << mesh1->verts[v1]->distortion << endl;}
						int cancel1 = mesh2->verts[v2]->matchIdx; //matched mesh1 vertex whose match is canceled below
						mesh1->verts[cancel1]->matchIdx = -1, mesh2->verts[v2]->matchIdx = -1; //cancel by resetting matchIdx values on both sides
						for (size_t c = 0; c < matches1.size(); c++) //cancel by erasing canceled entries from both sides
							if (matches1[c] == cancel1)
								matches1.erase(matches1.begin() + c), matches2.erase(matches2.begin() + c), cancel1 = -2; //erase(begin() + c) erases matches1/2[c]
						if (cancel1 != -2)
							cout << "WARNING: entry not found in matches1/2\n\n";
					}
#else //keep only the best/minDistortion of 2+ matches targeting the same mesh2 sample
					if (mesh2->verts[v2]->matchIdx == -1) //save the distortion of the match hitting v2 in case another one hits it later (manyTo1); in such a collision keep the best/minDistortion match and cancel the other
						mesh2->verts[v2]->distortion = mesh1->verts[ mesh1->samples[s] ]->distortion, matches1.push_back(mesh1->samples[s]), matches2.push_back(v2), mesh1->verts[ mesh1->samples[s] ]->matchIdx = v2, mesh2->verts[v2]->matchIdx = mesh1->samples[s]; //dual operation
					else if (mesh2->verts[v2]->distortion > mesh1->verts[ mesh1->samples[s] ]->distortion) //existing match hitting v2 is worse/highDistortion than the current manyTo1 match so cancel the existing match and set the new match
					{
						//typically tmp is same as distortion but rarely/validly tmp < distortion may arise when mesh1 sample already has a good/small tmp value from a previous overlapping patch; use distortion here to discard that previous (currently irrelevant) patch
						int cancel1 = mesh2->verts[v2]->matchIdx; //matched mesh1 vertex whose match is canceled below
//cout << "cancel: " << cancel1 << " - " << v2 << " w/ disto: " << mesh1->verts[cancel1]->tmp << " == " << mesh1->verts[cancel1]->distortion << " > " << mesh1->verts[ mesh1->samples[s] ]->distortion << " hence new match: " << mesh1->samples[s] << " - " << v2 << "\n\n"; //new match has the better/smaller distortion of mesh1->samples[s]
						mesh1->verts[cancel1]->matchIdx = -1, mesh2->verts[v2]->matchIdx = -1; //cancel by resetting matchIdx values on both sides
						for (size_t c = 0; c < matches1.size(); c++) //cancel by erasing canceled entries from both sides
							if (matches1[c] == cancel1)
								matches1.erase(matches1.begin() + c), matches2.erase(matches2.begin() + c), cancel1 = -2; //erase(begin() + c) erases matches1/2[c]
						if (cancel1 != -2)
							cout << "WARNING: entry not found in matches1/2\n\n";
						//canceled the existing match and now set the new match
						mesh2->verts[v2]->distortion = mesh1->verts[ mesh1->samples[s] ]->distortion, matches1.push_back(mesh1->samples[s]), matches2.push_back(v2), mesh1->verts[ mesh1->samples[s] ]->matchIdx = v2, mesh2->verts[v2]->matchIdx = mesh1->samples[s]; //dual operation
					}
#endif
				}
//cout << v2 << "\t\t" << mesh1->samples[s] << " - " << mesh1->verts[ mesh1->samples[s] ]->matchIdx << " set\n";
			} //end of s
			//rematch the ones (in 1to1 manner) that are unmatched above
			for (size_t i = 0; i < matches1Cpy.size(); i++)
				if (mesh1->verts[ matches1Cpy[i] ]->matchIdx == -1 && mesh2->verts[ matches2Cpy[i] ]->matchIdx == -1) //2nd condition prevents manyTo1 assignments
					mesh1->verts[ matches1Cpy[i] ]->matchIdx = matches2Cpy[i], mesh2->verts[ matches2Cpy[i] ]->matchIdx = matches1Cpy[i], //cout << "rematch: " << matches1Cpy[i] << " - " << matches2Cpy[i] << "\n\n";
					matches1.push_back(matches1Cpy[i]), matches2.push_back(matches2Cpy[i]); //matches2Cpy[i].distortion = matches1Cpy[i].distortion not made here 'cos matches1Cpy[i] is unknown; no problem 'cos trimMap() resets mesh1 distortions and mesh2 distortions are unimportant at this point
			if (trim)
			{
				qual = trimMap(/*2.0 1.6 1.2*/1.2, mesh1->samples, noTrimCoarseMapMatches, patchCostDistortion); //trimmedSamples[] is updated so it's ready for the cancellations below
				for (size_t s = 0; s < trimmedSamples.size(); s++) //cancel by erasing canceled entries from both sides
					for (size_t c = 0; c < matches1.size(); c++) //cancel by erasing canceled entries from both sides
						if (matches1[c] == trimmedSamples[s])
							matches1.erase(matches1.begin() + c), matches2.erase(matches2.begin() + c); //erase(begin() + c) erases matches1/2[c]
				for (size_t r = 0; r < rd.size(); r++) //trimMap() changed v.distortion values to unnormalized, and possibly irrelevant (trimType), ones; restore the normalized values for the upcoming tmp = distortion above
					mesh1->verts[ rd[r] ]->distortion = mesh1->verts[ rd[r] ]->geoToRef;
			}
//update coarseMap[] too; but it only improves w/ fillQdenseC2F().Q2Ddense, e.g., pathSamplingCoarse() is already done (i may recompute it if coarseMap changes but not doing that currently; not that critical)
			/*for (size_t i = 0; i < mesh1->coarseMap.size(); i += 2) disabled 'cos there's a final 'original coarseMap changed during ' print below anyway
				if (mesh1->verts[ mesh1->coarseMap[i] ]->matchIdx != mesh1->coarseMap[i+1])
					cout << "original coarseMap changed at level " << level << ": " << mesh1->coarseMap[i] << " - " << mesh1->coarseMap[i+1] << " --> " << mesh1->verts[ mesh1->coarseMap[i] ]->matchIdx << "\n"; //not always a bad thing; 595 - 599 --> 49003 was a good change for instance*/
			cout << matches1.size() << " matches (" << (qual == 1 ? "puregeo" : "pointwisepairwise") << " disto: " << distortion << ") after level " << level << " of radius " << radius << "\n\n"; //global distortion is updated just for the print here
			
			radius /= 2; //smaller radius (after more levels) can always find nSamplesPerSubregion samples to match but those findings will not respect radius at all, i.e., high nOffRadius value, and will be less meaningful (no pool/gc stuff) so don't decrease radius at all (verified that this leads to bad results so decrease radius)
//			radius /= (level < 1 ? 2 : 1); //depends on the loop counter (prevents radius < 0.25 as it divides the initial radius=0.5 only once when level=0)
//			radius /= (level <= 3 ? 2 : 1); //depends on the loop counter (prevents radius < 0.25 as it divides the initial radius=0.5 only once when level=0)

			if (backToOptimal) //update optimalMap, i.e., my final result may be from level 2 despite nMaxLevels=5
			{
				minDisto = (matches1.size() > (optimalMap.size() / 2) ? INF : minDisto); //start optimalMap search again (via =INF) if the current number of matches is bigger than optimalMap size, i.e., get the densest map
				if (distortion < minDisto) //preferring the densest map is also fair 'cos distortion is not an avg it is a total so sparser maps have unfairly low distortion
				{
					minDisto = distortion, optimalMap.clear();
					for (size_t i = 0; i < matches1.size(); i++) //transfer matches1/2 to optimalMap
						optimalMap.push_back(matches1[i]), optimalMap.push_back(matches2[i]); //2 consecutive entries make 1 match
				}
				if (matches1.size() >= nMaxMatches)
				{
					minDistoDensest = (matches1.size() > nMaxMatches ? INF : minDistoDensest); //start densestOptMap search again (via =INF) if the current number of matches is bigger than nMaxMatches so far, i.e., get the densest map
					nMaxMatches = matches1.size();
					if (distortion < minDistoDensest)
					{
						minDistoDensest = distortion, densestOptMap.clear();
						for (size_t i = 0; i < matches1.size(); i++) //transfer matches1/2 to densestOptMap
							densestOptMap.push_back(matches1[i]), densestOptMap.push_back(matches2[i]); //2 consecutive entries make 1 match
					}
				}
			}
//if(level==0)break;
		} //end of level
		for (size_t s = 0; s < mesh1->samples.size(); s++)
		{
			if (mesh1->verts[ mesh1->samples[s] ]->geoToRef == -1.0 && printOnce1 && ! curvFPS) //normal if curvFPS in use: " << mesh1->verts[ mesh1->samples[s] ]->gc << " vs. " << mesh1->avgGC ('cos small-curvature samples can never win the pool competition)
				cout << "WARNING: c2f didn't touch mesh1's " << mesh1->samples[s] << " at all (there may be others; printing stopped)\n\n", printOnce1 = false;; //small-curvature samples can never win the pool competition (disabling curvFPS on level2+ is therefore a logical heuristic)			
			//mesh2 side which has the same number of samples for sure
			if (mesh2->verts[ mesh2->samples[s] ]->geoToRef == -1.0 && printOnce2 && ! curvFPS)
				cout << "WARNING: c2f didn't touch mesh2's " << mesh2->samples[s] << " at all (there may be others; printing stopped)\n\n", printOnce2 = false;
			if (backToOptimal)
				mesh1->verts[ mesh1->samples[s] ]->matchIdx = mesh2->verts[ mesh2->samples[s] ]->matchIdx = -1; //optimalMap below will reset these values below
		}
		if (backToOptimal) //optimalMap to matchIdxs
		{
			if (minDisto < minDistoDensest)
				for (size_t i = 0; i < optimalMap.size(); i += 2)
					mesh1->verts[ optimalMap[i] ]->matchIdx = optimalMap[i + 1], mesh2->verts[ optimalMap[i + 1] ]->matchIdx = optimalMap[i]; //dual operation
			else //unlikely 'cos more matches mean more additions to distortion (no avg) but still rarely observed minDistoDensest <= minDisto
				for (size_t i = 0; i < densestOptMap.size(); i += 2)
					mesh1->verts[ densestOptMap[i] ]->matchIdx = densestOptMap[i + 1], mesh2->verts[ densestOptMap[i + 1] ]->matchIdx = densestOptMap[i]; //dual operation
			cout << "back to optimal distortion: "; if (patchCostDistortion && dissimType != PAIRWISENOPATCH) mapQualityPatch(); else mapQualityGeo(); cout << distortion << endl; //distortion is updated by the mapQualityPatch/Geo() here*/
		}
		matchUnmatcheds(nPathSamples, nPathSamplesExtra); //7 out of 30 remained unmatched in nMaxLevels=5; make 7x7 map and trim bad ones, e.g., match additional 3x3 good subset from the unmatched list
	}	
	if (mapType != LAP) //evaluation of QAP or C2FQAP
	{		
		int nMatches;
		sprintf(str, " -qap- "), sprintf(orgDtype, dType), sprintf(dType, "puregeo"), nMatches = mapQualityGeo(), printMap(t, str), sprintf(dType, orgDtype); //mapQualityGeo() and hence puregeo distortions will be printed so update dType first and then roll back to original (char str[] = " -qap- "; same effect)
		for (size_t v = 0; v < mesh1->verts.size(); v++)
			mesh1->verts[v]->marked = false; //true means display it w/ a sphere
		bool printOnce = true;
		if (displayCoarse) //display coarseMap
			for (size_t i = 0; i < mesh1->coarseMap.size(); i += 2) //display coarseMap matches
			{
				mesh1->verts[ mesh1->coarseMap[i] ]->marked = true;
				if (mesh1->verts[ mesh1->coarseMap[i] ]->matchIdx != mesh1->coarseMap[i+1] && printOnce)
					cout << "original coarseMap changed during dense mapping: " << mesh1->coarseMap[i] << " - " << mesh1->coarseMap[i+1] << " --> " << mesh1->verts[ mesh1->coarseMap[i] ]->matchIdx << " (there may be others; printing stopped)\n\n", printOnce = false;
//				mesh1->verts[ mesh1->coarseMap[i] ]->matchIdx = mesh1->coarseMap[i+1]; //undoes qap's matchIdx for original coarseMap visualization; qucoop.objective and mapQualityVirtual.distortion will print differently so disable this line (for consistent printing only)
			}
		else  //display bad matches
			for (size_t s = 0; s < mesh1->samples.size(); s++)
				if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1 && mesh1->verts[ mesh1->samples[s] ]->distortion >= 1.5*distortion / nMatches) //distortion not divided by number of matches so do it here
					mesh1->verts[ mesh1->samples[s] ]->marked = true; //display bad matches (display=1 means display coarseMap and =2 means display bad matches)
//					mesh2->verts[ mesh1->verts[ mesh1->samples[s] ]->matchIdx ]->matchIdx = -1, mesh1->verts[ mesh1->samples[s] ]->matchIdx = -1, qucoop.objective and mapQualityVirtual.distortion will print differently so disable this trimming (for consistent printing only)
		if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
			mapQualityPatch(), cout << "distortion=" << dType << ": " << distortion; //else mapQualityGeo(); redundant 'cos mapQualityGeo/printMap() is called above to learn nMatches and print geodesic distortions (mapQualityVirtual() uses the same metrics so call either of them here)
		else
			cout << "WARNING: Q is filled w/ lengths to coarseMap but mapQualityGeo() is not so it's OK to get different objective & distortion prints\n\n";
		cout << "; one " << nSamplesPerSubregion << " x " << nSamplesPerSubregion  << (brute ? " brute-force" : " QAP") << " solution took " << (((float) clock()-t) / CLOCKS_PER_SEC) / nQAPs << " secs\n";
	}
}

int Correspondence::trimMap(double f, const vector< int >& samples, bool coarseMapMatchesExempt, bool patchCostDistortion, bool print)
{
	//trims very bad matches based on puregeo and then pointwisepairwise (if applicable) distortions and print the latest distortion; samples is either all mesh1.samples (currentl) or subset mesh1.samplesC2F

	int qual = 0, trimType = (! patchCostDistortion ? 1 : 2/*3*/); //1: only puregeo-based trimming, 2: only pointwisepairwise-based trimming (recommended), 3: puregeo-based followed by pointwisepairwise-based trimming; dissimType!=PAIRWISENOPATCH must hold to be able to use pointwisepairwise-based trimming
	trimmedSamples.clear(); //list of samples whose matches are trimmed (matches1/2 of C2FQAP adjust themselves based on this global)	
	double n = 0, nMatches, disto;	
	if (trimType == 1 || trimType == 3 || (trimType == 2 && (dissimType == PAIRWISENOPATCH || dissimType == POINTWISE))) //last condition to make trimming when trimType=2 but path/patch disabled due to POINTWISE||PAIRWISENOPATCH selection
	{
		qual = 1, nMatches = mapQualityGeo(), disto = distortion / nMatches; //update v.distortion and global distortion values based on the puregeo metric and set disto to make distortion compatible with v.distortion			
if (print) for (size_t s = 0; s < samples.size(); s++) if (mesh1->verts[ samples[s] ]->matchIdx != -1) cout << samples[s] << " -beforeTrim- " << mesh1->verts[ samples[s] ]->matchIdx << "\t" << mesh1->verts[ samples[s] ]->distortion << " vs. avgPureGeo: " << disto << endl;
		for (size_t s = 0; s < samples.size(); s++)
			if (mesh1->verts[ samples[s] ]->matchIdx != -1 && mesh1->verts[ samples[s] ]->distortion > f*disto)
			{
				if (coarseMapMatchesExempt && mesh1->verts[ samples[s] ]->coarseMapMatch)
					continue; //coarseMap matches (from main().qapMap()) are exempt from trimming hence cannot be trimmed
//cout << samples[s] << " -puregeoTrimmed- " << mesh1->verts[ samples[s] ]->matchIdx << "\t" << mesh1->verts[ samples[s] ]->distortion << " vs. avgPureGeo: " << disto << endl,
				mesh2->verts[ mesh1->verts[ samples[s] ]->matchIdx ]->matchIdx = -1, mesh1->verts[ samples[s] ]->matchIdx = -1, trimmedSamples.push_back(samples[s]), n++;
			}
		mapQualityGeo(); //update v.distortion and global distortion using on the trimmed most recent map
	}
	if ((dissimType == POINTWISEPAIRWISE || dissimType == PAIRWISE) && trimType > 1)
	{
		qual = 2, nMatches = mapQualityPatch(), disto = distortion / nMatches; //make distortion compatible with v.distortion (mapQualityVirtual() uses the same metrics so call either of them here)
if (print) for (size_t s = 0; s < samples.size(); s++) if (mesh1->verts[ samples[s] ]->matchIdx != -1) cout << samples[s] << " -beforeTrim- " << mesh1->verts[ samples[s] ]->matchIdx << "\t" << mesh1->verts[ samples[s] ]->distortion << " vs. avgPntwisePairwise: " << disto << endl;
		for (size_t s = 0; s < samples.size(); s++)
			if (mesh1->verts[ samples[s] ]->matchIdx != -1 && mesh1->verts[ samples[s] ]->distortion > f*disto)
			{
				if (coarseMapMatchesExempt && mesh1->verts[ samples[s] ]->coarseMapMatch)
					continue; //coarseMap matches (from main().qapMap()) are exempt from trimming hence cannot be trimmed
//cout << samples[s] << " -pointwisepairwiseTrimmed- " << mesh1->verts[ samples[s] ]->matchIdx << "\t" << mesh1->verts[ samples[s] ]->distortion << "\n",
				mesh2->verts[ mesh1->verts[ samples[s] ]->matchIdx ]->matchIdx = -1, mesh1->verts[ samples[s] ]->matchIdx = -1, trimmedSamples.push_back(samples[s]), n++; //samples whose matches are trimmed
			}
		mapQualityPatch(); //update v.distortion and global distortion using on the trimmed most recent map
	}
	if (n > 0)
		cout << n << " match" << (n > 1 ? "es" : "") << " of this level " << (n > 1 ? "are" : "is") << " trimmed\n";
//qual=1,mapQualityGeo();//print puregeo on return
	return qual;
}

void Correspondence::c2fHelper(int src, int tgt, double radius, bool curvFPS, int nPathSamples, int nPathSamplesExtra, size_t cSize, bool patchCostDistortion, int& qual, vector<int>& rd, bool secondCall) //modifiable references: changes apply to the original int and vector
{
	//helps matching the current level samples during the coarse-to-fine process

	//reset matchIdx values for all samples as qapMap() below will update them for the current FPSpatch() samples and then i'll evaluate matches via mapQualityPatch/Geo() which requires matchIdx=-1 for the non-FPSpatch() samples
//	for (size_t s = 0; s < mesh1->samples.size(); s++) //k in mapQualityPatch().pathSamples[k] requires fixed samples[] indexing so i cannot do the samples1Cpy tactic below; this is the accurate way and aligns well with my distortion normalization too
//		mesh1->verts[ mesh1->samples[s] ]->matchIdx = mesh2->verts[ mesh2->samples[s] ]->matchIdx = -1; //set by qapMap() below and then used by mapQualityPatch/Geo() below; mesh2 size is necessary to prevent 'entry not found' warning
	for (size_t v = 0; v < mesh1->verts.size(); v++) mesh1->verts[v]->matchIdx = -1; for (size_t v = 0; v < mesh2->verts.size(); v++) mesh2->verts[v]->matchIdx = -1; //reset over verts 'cos sample set changes during denseMapLargeN() so be safe (no such problem in denseMap() so reset over samples ok for denseMap())

	//FPSpatch()'ll use verts[ samplesC2F[e] ]->geodesic[s] as samplesC2F[] is populated based on the current k=30 samples[] so normalize all samples[].geodesics as coarseMap[].geodesics are already normalized and hence prevent unnormalized vs normalized incompatibilities
	for (size_t s = 0; s < mesh1->samples.size(); s++) //redundant for denseMap() 'cos samples[] definitely normalized via unitGeoScaling(); necessary if called from denseMapLargeN()
		if (! mesh1->verts[ mesh1->samples[s] ]->geoNormalized)
		{
			for (size_t j = 0; j < mesh1->verts.size(); j++)
				mesh1->verts[ mesh1->samples[s] ]->geodesic[j] /= mesh1->maxGeoDist;
			mesh1->verts[ mesh1->samples[s] ]->geoNormalized = true; //currently normalized
		}
	for (size_t s = 0; s < mesh2->samples.size(); s++) //redundant for denseMap() 'cos samples[] definitely normalized via unitGeoScaling(); necessary if called from denseMapLargeN()
		if (! mesh2->verts[ mesh2->samples[s] ]->geoNormalized)
		{
			for (size_t j = 0; j < mesh2->verts.size(); j++)
				mesh2->verts[ mesh2->samples[s] ]->geodesic[j] /= mesh1->maxGeoDist;
			mesh2->verts[ mesh2->samples[s] ]->geoNormalized = true; //currently normalized
		}
	mesh1->FPSpatch(nSamplesPerSubregion, radius, src, curvFPS, secondCall), mesh2->FPSpatch(nSamplesPerSubregion, radius, tgt, curvFPS, secondCall); //mesh1/2.samplesC2F of size nSamplesPerSubregion ready now
	//upcoming fillQdenseC2F(), mapQualityPatch(), matchUnmatcheds() need normalized geodesics so normalize them if not normalized already (denseMap() definitely normalized all the samples' geodesics via unitGeoScaling() and s-loop above fixed the issue for denseMapLargeN())

	//define the patches by sampling nSamplesPerSubregion samples inside them (via radius restriction and v.roi restriction where the latter is implied as i select a subset of samples from samples[] which are definitely v.roi=true)
	if (! mesh1->extraSamples) //nPathSamples=-1: means special flag for Q2Ddense usage (diag: lengthsToCoarseMap off-diag: lengthsBetweenSamples); nPathSamples>0: means diag: ringAreas off-diag: patchAreasBetweenSamples (alternative not in use: diag: patchAreasToCoarseMap off-diag: patchAreasBetweenSamples)
		mesh1->fillQdenseC2F(nPathSamples + 1, areaWeight), mesh2->fillQdenseC2F(nPathSamples + 1, areaWeight), //QAP needs the Q matrix
		qapMap(nPathSamples + 2 + cSize, true); //nPathSamples=-1 sends 1 which is a special flag for Q2D usage; otherwise nPathSamples + 1 gives # path samples and another + 1 gives access to the additional dimension storing pure geodesic length info (for maxGeoDiff filtering or geoDiff-based barrier) and another + cSize for lengths to coarseMap entries
	else //nPathSamples!=-1 for sure 'cos nPathSamples=-1 means no path sampling requested, i.e., POINTWISE or PAIRWISENOPATCH, and in that case extraSamples cannot be true
		mesh1->fillQdenseC2F(nPathSamplesExtra, areaWeight), mesh2->fillQdenseC2F(nPathSamplesExtra, areaWeight), //+1 added during declaration so no nPathSamplesExtra + 1 here
		qapMap(nPathSamplesExtra + 1 + cSize, true); //+ 1 gives access to the additional dimension storing pure geodesic length info (for maxGeoDiff filtering or geoDiff-based barrier) and another + cSize for lengths to coarseMap entries
	//evaluate the matches produced by qapMap() above; v.distortion ready now for mesh1.samplesC2F				
//	samples1Cpy = mesh1->samples, samples2Cpy = mesh2->samples, mesh1->samples = mesh1->samplesC2F, mesh2->samples = mesh2->samplesC2F; before evaluation reset mesh1/2.samples as mapQualityPatch/Geo() iterates over mesh1/2.samples but here i want to restrict the eval to mesh1/2.samplesC2F (to be compatible w/ the normalized tmp below); no see samples1Cpy tactic comment above
	if (patchCostDistortion && dissimType != PAIRWISENOPATCH) //pathing/patching in use
		mapQualityPatch(), qual = 2; //global distortion is updated
	else
		mapQualityGeo(), qual = 1; //global distortion is updated
//	mesh1->samples = samples1Cpy, mesh2->samples = samples2Cpy; roll back to the original mesh1/2.samples (see samples1Cpy tactic comment above)
	for (size_t j = 0; j < mesh1->samplesC2F.size(); j++) //best match update (binID, naming sucks) if the distortion permits (tmp, naming sucks)
	{
		mesh1->verts[ mesh1->samplesC2F[j] ]->distortion /= nSamplesPerSubregion; //distortion normalization to make it compatible w/ tmp from above
		mesh1->verts[ mesh1->samplesC2F[j] ]->geoToRef = mesh1->verts[ mesh1->samplesC2F[j] ]->distortion, rd.push_back(mesh1->samplesC2F[j]); //save normalized distortion (based on 5x5 qapMap() not main().qapMap()) as trimMap() below will change them (naming sucks)
//cout << mesh1->verts[ mesh1->samplesC2F[j] ]->distortion << " vs. " << mesh1->verts[ mesh1->samplesC2F[j] ]->tmp << endl;
//if (mesh1->samplesC2F[j] == 25979 || mesh1->samplesC2F[j] == 595) cout << level << ": " << mesh1->verts[ mesh1->samplesC2F[j] ]->distortion << " vs. " << mesh1->verts[ mesh1->samplesC2F[j] ]->tmp << " due to " << mesh1->samplesC2F[j] << " - " << mesh1->verts[ mesh1->samplesC2F[j] ]->matchIdx << endl;
		if (mesh1->verts[ mesh1->samplesC2F[j] ]->distortion < mesh1->verts[ mesh1->samplesC2F[j] ]->tmp) //samplesC2F[j] may already have a good/small tmp value from a previous overlapping patch so don't update that good match which has smaller distortion/tmp (naming sucks)
			mesh1->verts[ mesh1->samplesC2F[j] ]->tmp = mesh1->verts[ mesh1->samplesC2F[j] ]->distortion, mesh1->verts[ mesh1->samplesC2F[j] ]->binID = mesh1->verts[ mesh1->samplesC2F[j] ]->matchIdx;
		mesh2->verts[ mesh2->samplesC2F[j] ]->geoToRef = 0.0; //any value different than -1 states that mesh2.samplesC2F[j] is touched now (for printing and FPSpatch() only) mesh1/2.samplesC2F.size() same so no problem here
	}
}

void Correspondence::denseMapLargeN(size_t k)
{
	//same as denseMap() except this ones samples new k=30, matches them, and then repeats until n=250 (large n) samples are matched
	
	vector< int > prevSamples1, prevSamples2, coarseMapFull; //coarseMap[] is the trimmed version of this original coarseMapFull formed by main().qapMap() matches (used by ROI() and distortion averaging only)
	for (size_t s = 0; s < mesh1->samples.size(); s++)
		if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
			coarseMapFull.push_back(mesh1->samples[s]), coarseMapFull.push_back(mesh1->verts[ mesh1->samples[s] ]->matchIdx); //each 2 consecutive items make 1 match (even indices belong to mesh1, odd indices mesh2)
	//make coarseMap ready for blossom costType > 1 (also for pathSamplingCoarse() and fillQdense2/3() and may be other things)
	for (size_t s = 0; s < mesh1->samples.size(); s++)
		if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1 && mesh1->verts[ mesh1->samples[s] ]->distortion < 1.1*distortion/(coarseMapFull.size()/2)) //2nd condition to keep coarseMap very small (efficiency) and accurate (accuracy in, e.g., pathSamplingCoarse() and tmp tests) by admitting distortions below the 1.1avg distortion (1.1 'cos some matches must be above average but they may still be very good)
			mesh1->coarseMap.push_back(mesh1->samples[s]), mesh1->coarseMap.push_back(mesh1->verts[ mesh1->samples[s] ]->matchIdx), mesh1->verts[ mesh1->samples[s] ]->coarseMapMatch = mesh2->verts[ mesh1->verts[ mesh1->samples[s] ]->matchIdx ]->coarseMapMatch = true; //even and odd idxs mesh1 and mesh2 side, respectively
		else if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1) //high distortion, i.e., in coarseMapFull[] but not in coarseMap, so reset matchIdx
			mesh2->verts[ mesh1->verts[ mesh1->samples[s] ]->matchIdx ]->matchIdx = -1, mesh1->verts[ mesh1->samples[s] ]->matchIdx = -1; //ROI() uses coarseMapFull[] but not need these matchIdx info so no problem
	mesh2->coarseMap = mesh1->coarseMap; //coarseMap[] will not be updated but if it was updated then a change in mesh2.coarseMap would not affect mesh1.corresp (verified)
	int iter = 0, nMatches, m = mesh2->samples.size(); //n mesh1 samples are matched w/ m mesh2 samples (m <= n) during main.qapMap() so remember m to adjust areaWeight as follows:
	//3rd dimension of Q may change, i.e., number of ring areas on a sample or number of patch areas between 2 samples, but it does not affect the weighting b/w unary and pairwise terms (each term gets the L2 difference (of areasRing[] for instance) w/ new 3rd dimension so no problem); the change on mesh2 sample size,
	//however, affect the weighting (areaWeight) since 1 unary + m-1 pairwise are added from each Q row to get the overall distortion; new m is nSamplesPerSubregion=5 in C2FQAP mode and nDenseSamples=30 in QAP mode; to equalize unary term's effect (see last line in setAreaWeight()) multiply areaWeight by the newM/m ratio
//nDenseSamples = 100;//250;//150;//50;//30; //nDenseSamples=30&nPathSamples=4 takes 62secs, =50&4 takes 140secs so reduce 4->2 below if nDenseSamples > 30, =50&1 takes 6secs, =50&2 takes 63secs

	mesh1->eucPatching = (mesh1->verts.size() > 20000);//19000); //geodesic-based patching, despite the early break in patchAndDescriptors(), is significantly slow (13secs on 52K-vertex mesh w/ nSamples=nPathSamples=10 (eucPatching takes just 0.3secs), 20secs nSamples=30,nPathSamples=3) so use euclidean distance during patching for a 40x speedup
	mesh2->eucPatching = (mesh2->verts.size() > 20000);//19000); //geodesic-based patching, despite the early break in patchAndDescriptors(), is significantly slow (13secs on 52K-vertex mesh w/ nSamples=nPathSamples=10 (eucPatching takes just 0.3secs), 20secs nSamples=30,nPathSamples=3) so use euclidean distance during patching for a 40x speedup
	if (mesh1->eucPatching || mesh2->eucPatching)
		cout << "WARNING: path patches by euclidean distance for efficiency; geodesic is more accurate but slow\n\n";
	//region of interests based on the correspondence information
	double r2 = mesh2->ROI(coarseMapFull), r1 = mesh1->ROI(coarseMapFull, r2); //mesh2 is the partial model whose roiArea r2 gives an upper bound for the roiArea of the complete mesh1 (for efficiency only as ROI() is costly for dense, e.g., 52K-vertex, models); no problem for complete matching scenario
	if (max(r1/r2, r2/r1) > 2)
		cout << "WARNING: source and target ROI areas too different, implying bad map (recompute the map)\n\n";	
	mesh1->nonuniformPathSamplingHeuristic = mesh2->nonuniformPathSamplingHeuristic = 2;//0; //might be changed to 1 in main(), e.g., david0-10_partial3, so return to default 2 to be safe
	clock_t t = clock(); //includes FPSroi() and pathSampling() timings so the individual timing prints for these functions below are redundant but anyway
	do
	{
		mesh1->FPSroi(k, prevSamples1), mesh2->FPSroi(k, prevSamples2); //put k=30 FPS samples on ROIs while discarding the previous samples (they affect the locations of new samples but they are discarded from being samples[] of this iteration and hence cannot be matched, i.e., retain their existing matchIdx values); samples[] ready now
		size_t cSize = mesh1->coarseMap.size() / 2; //number of matches in coarseMap[] (2 consecutive entries make 1 match)
		int nPathSamples = (dissimType == POINTWISE || dissimType == PAIRWISENOPATCH ? -1 : 4); //paths from dense samples to coarseMap will mostly be short so use small radius and small nPathSamples (4 (makes total 4+1 samples); too slow if >4)
//nPathSamples = ((nDenseSamples > 50 || mesh1->verts.size() > 20000) && nPathSamples != -1 ? 2/*1*/ : nPathSamples); //disabled 'cos more samples are better as i go to distant coarseMap samples in extraSamples mode; nPathSamples=2: 21secs vs. 4: 33secs
		int nPathSamplesExtra = cSize * nPathSamples + nPathSamples+1; //3*4 (towards 3 coarseMap samples (not counting the common samples[j] end)) + 5 from samples[i] to samples[j] (actual path w/o the extra sampling)
		mesh1->extraSamples = mesh2->extraSamples = (dissimType == PAIRWISE || dissimType == POINTWISEPAIRWISE ? mesh1->extraSamples : false); //extraSamples is ineffective (false) if no pathing/patching in use, i.e., nPathSamples = -1 (regardless of its constructor initialization)
		//pairwise structure
		if (dissimType == PAIRWISE || dissimType == POINTWISEPAIRWISE) //path sampling requested (nPathSamples!=-1 for sure)
		{
			if (mesh1->extraSamples)
				mesh1->pathSamplingCoarse(nPathSamples), mesh2->pathSamplingCoarse(nPathSamples); //make samplesToCoarse[] ready for the pathSampling() calls below
			mesh1->pathSampling(nPathSamples), mesh2->pathSampling(nPathSamples); //unnormalized euc/geo distance vs. unnormalized radius so geo must be unnormalized (no unitGeoScaling()) prior to pathSampling/pathSamplingCoarse() calls
		}
		//pointwise descriptor
		if (dissimType == POINTWISE || dissimType == POINTWISEPAIRWISE) //pointwise descriptor is set by areas() so dissimType must include pointwise (to save tiny time)
		{
			mesh1->areas(), mesh2->areas(); //radius has not changed since coarseMap[] computation but new dense samples arrived so call this function to fill their areasRing[] and areas[] entries
			if (mesh1->borderTreat == 1)
				mesh1->doubleBorderAreas(), mesh2->doubleBorderAreas(); //avgAreas() called inside doubleBorderAreas()
			else
				mesh1->avgAreas(), mesh2->avgAreas(); //set avgPatch/RingArea values used by areaWeight
		}
		//newM/m areaWeight adjustments (no barrier adjustment 'cos original main().barrier is based on avgPairwise value which does not change as i use the same radius as main() here; note that barrier is added to EACH pairwise term, via if (pairwise)Q(u,v)+=barrier*.., so it's quite effective)
		areaWeight *= ((double) nSamplesPerSubregion / m), printQAP = false, cout << "C2FQAP map computation via " << qapSolverStr << "..\n"; //no printing during C2FQAP as qapMap() will be called many times in that mode (newM/m adjustment explained above is done here)
#ifndef GEODIFFTHRESHOLDING		
		cout << "PARAMS (dense): radius = " << mesh1->radius << ", nRingAreas = " << mesh1->nRingAreas << ", nPathSamples = " << nPathSamples+1 << ", unaryWeight = " << areaWeight << ", barrier = " << (barrier > 0 ? barrier : 0) << (mesh1->weightByGC ? ", wByGC" : "") << ", areasRing, qapEnergy: " << dissimType << ", pathSampler: " << mesh1->nonuniformPathSamplingHeuristic << "\n"; //barrier n/a if negative
#else
		barrier = 0, cout << "PARAMS (dense): radius = " << mesh1->radius << ", nRingAreas = " << mesh1->nRingAreas << ", nPathSamples = " << nPathSamples+1 << ", unaryWeight = " << areaWeight << ", maxGeoDiff = " << maxGeoDiff << (mesh1->weightByGC ? ", wByGC" : "") << ", areasRing, qapEnergy: " << dissimType << ", pathSampler: " << mesh1->nonuniformPathSamplingHeuristic << "\n"; //non-overlapping areasRing[] alternative is overlapping is areas[] for the unary term
#endif		
		cSize = (nPathSamples == -1 ? 0 : cSize); //disable cSize effect to be able to send special flag -1 for Q2D usage (in no patch mode of POINTWISE or PAIRWISENOPATCH) (nPathSamples!=-1 for sure for the else extraSamples part below)

		//quadratic assignment problem via qucoop (coarse-to-fine approach)
		//each matched sample (initially only coarseMap[] samples are matched and hence inside matches1) defines a patch which is resampled via FPSpatch() and matched via qapMap()
		vector< int > matches1, matches2; //nMaxLevels many levels will be used where the number of matches increases as levels increase
		double minDisto = INF, minDistoDensest = INF, nMaxMatches = -INF, radius = 1.0 / 2.0; //mesh1->maxGeoDist!=1 'cos no unitGeoScaling() in this function so use 1.0 explicitly
		vector< int > optimalMap, densestOptMap; //map corresponding to the minimum-distortion map or densest minimum-distortion map through c2f levels
		bool patchCostDistortion = true, //false; //best match update based on patch-based pointwisepairwise (true and also dissimType != PAIRWISENOPATCH) or puregeo (false) distortions
			 printOnce1 = true, printOnce2 = true, backToOptimal = true, doTwice = false, noTrimCoarseMapMatches = true; //matches in the initial coarseMap[], i.e., sparse/coarse matches from main().qapMap() cannot be trimmed, i.e., they may change/refine but never be -1
		cout << "c2f w/ initial c2fRadius: " << radius << ", nSamplesPerSubregion: " << nSamplesPerSubregion << ", " << "back2optimal: " << backToOptimal << ", doTwice: " << doTwice << endl;
//for (size_t v = 0; v < mesh1->verts.size(); v++) if (mesh1->verts[v]->matchIdx != -1) cout << v << "m1\n"; //27 out of 30 matches of the previous do-while iteration printed
//for (size_t v = 0; v < mesh2->verts.size(); v++) if (mesh2->verts[v]->matchIdx != -1) cout << v << "m2\n"; //27 out of 30 matches of the previous do-while iteration printed
		for (size_t v = 0; v < mesh1->verts.size(); v++)
			mesh1->verts[v]->matchIdx = -1;
		for (size_t v = 0; v < mesh2->verts.size(); v++)
			mesh2->verts[v]->matchIdx = -1;
		//cannot call unitGeoScaling() as it would have divided already-normalized geodesics further into tiny values
		//mapQualityPatch() below (before the nMaxLevels iterations) need normalized geodesics based on mesh1->maxGeoDist set in the very first FPSroi() call above
		for (size_t i = 0; i < mesh1->coarseMap.size(); i += 2) //each 2 consecutive items make 1 match (even indices belong to mesh1, odd indices mesh2)
		{
			int v1 = mesh1->coarseMap[i], v2 = mesh1->coarseMap[i+1];
			mesh1->verts[v1]->matchIdx = v2, mesh2->verts[v2]->matchIdx = v1; //prepare matchIdx values for mapQualityPatch/Geo() below (before the nMaxLevels iterations) based on coarseMap[]; dual operation redundant but anyway
			if (mesh1->verts[v1]->geoNormalized || mesh2->verts[v2]->geoNormalized || mesh1->verts[v1]->geodesic.empty() || mesh2->verts[v2]->geodesic.empty())
				cout << "WARNING: coarseMap[] samples had unexpected normalized values or emptiness issue\n\n";
			for (size_t j = 0; j < mesh1->verts.size(); j++)
				mesh1->verts[v1]->geodesic[j] /= mesh1->maxGeoDist; //very first iteration sets normalizer mesh1->maxGeoDist based on the very first k=30 samples which definitely includes the most extreme/farthest points
			for (size_t j = 0; j < mesh2->verts.size(); j++)
				mesh2->verts[v2]->geodesic[j] /= mesh1->maxGeoDist;
			mesh1->verts[v1]->geoNormalized = mesh2->verts[v2]->geoNormalized = true; //currently normalized
			if (! prevSamples1.empty()) //very first iteration filled matches1/2 below and learned distortionCpy values that are used here later
				matches1.push_back(v1), matches2.push_back(v2), mesh1->verts[v1]->distortion = mesh1->verts[v1]->distortionCpy; //matches1[] is based on coarseMap[] now (it'll expand below, possibly changing the matches of the previous coarse levels)
				




//possibly changing the matches of the previous coarse levels::::::::::::::::: what happens to initial coarseMap[] in the end???????? i mean coarseMap[] never updated but what about the matches it implies?????????
//i mean each do-while iteration is supposed to touch a unique sample but now that coarseMap[] is in matches1/2[] it'll touch coarseMap[] samples multiple times, which one to keep?





		}
		if (prevSamples1.empty()) //very first iteration has definitely coarseMap[] inside current samples[] (due to clearSamples()) so mapQualityPatch/Geo() can be used to evaluate its distortion (there're other samples[] but their matchIdxs are -1 and hence not used by mapQualityPatch/Geo())
		{
			//evaluate the matches produced by the very first main().qapMap() using the new parameters shown in PARAMS print above and/or new coarseMap that is trimmed above; v.distortion ready now for level=0 case below (other levels done naturally w/ the mapQualityPatch/Geo() calls inside the loop)
			if (patchCostDistortion && dissimType != PAIRWISENOPATCH) //pathing/patching in use
				nMatches = mapQualityPatch(); //global distortion is updated
			else
				nMatches = mapQualityGeo(); //global distortion is updated (redundant 'cos PARAMS update above does not affect puregeo distortion in mapQualityGeo() that does not use coarseMap either)
			for (size_t s = 0; s < mesh1->samples.size(); s++)
			{
				if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
					mesh1->verts[ mesh1->samples[s] ]->distortion /= nMatches, //distortion normalization as the current map size is coarseMap.size() which is not necessarily compatible w/ the nSamplesPerSubregion map size below (normalize those distortions too via /=nSamplesPerSubregion below)
					mesh1->verts[ mesh1->samples[s] ]->distortionCpy = mesh1->verts[ mesh1->samples[s] ]->distortion, //to be used by the next do-while iteration
					matches1.push_back(mesh1->samples[s]), matches2.push_back(mesh1->verts[ mesh1->samples[s] ]->matchIdx); //matches1[] is based on coarseMap[] now (it'll expand below, possibly changing the matches of the previous coarse levels)
				mesh1->verts[ mesh1->samples[s] ]->geoToRef = mesh2->verts[ mesh2->samples[s] ]->geoToRef = -1.0; //to detect the untouched samples during nMaxLevels processing below; hopefully all samples are touched/processed/matched; trimming may leave some samples unmatched which is OK but they must be touched at least once during the nMaxLevels processing; -1 means untouched (for printing, FPSpatch() and distortion restoration)
			}
		}
		else //for the other iterations, matches1/2[] and s.distortion are already filled by the i-loop above so just init the geoToRef here
			for (size_t s = 0; s < mesh1->samples.size(); s++)
				mesh1->verts[ mesh1->samples[s] ]->geoToRef = mesh2->verts[ mesh2->samples[s] ]->geoToRef = -1.0; //to detect the untouched samples during nMaxLevels processing below; hopefully all samples are touched/processed/matched; trimming may leave some samples unmatched which is OK but they must be touched at least once during the nMaxLevels processing; -1 means untouched (for printing, FPSpatch() and distortion restoration)
		for (int level = 0; level < nMaxLevels; level++)
		{
			for (size_t s = 0; s < mesh1->samples.size(); s++)
			{
				mesh1->verts[ mesh1->samples[s] ]->binID = -1, mesh1->verts[ mesh1->samples[s] ]->tmp = INF; //binID is the best matchIdx throughout this level's iterations (naming sucks) and tmp is that match's distortion (naming sucks)
				if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1) //mesh1->samples[s] already has a distortion from coarseMap by qapMap() (level=0 case) so use that distortion as tmp to naturally protect the good coarseMap matches (different from noTrimCoarseMapMatches which protects them from being trimmed; this line protects them from being overwritten inside samplesC2F loop below)
					mesh1->verts[ mesh1->samples[s] ]->tmp = mesh1->verts[ mesh1->samples[s] ]->distortion;//cout << mesh1->samples[s] << " -coarseMap- " << mesh1->verts[ mesh1->samples[s] ]->matchIdx << "\t" << mesh1->verts[ mesh1->samples[s] ]->distortion << dType << "normalized\n";
			}
			vector< int > rd; //verts whose distortions will be restored to the saved normalized distortions in geoToRef (naming sucks)
			int qual = 0; //most recent quality metric (mapQualityGeo/Patch()) for printing only
			for (size_t i = 0; i < matches1.size(); i++) //same as matches2.size()
			{
				//do each radius twice, first w/ typical FPSpatch and second w/ modified FPSpatch that selects untouched (geoToRef=-1, naming sucks) samples (to prevent many unmatched matches in the end due to the untouched samples)
				c2fHelper(matches1[i], matches2[i], radius, true, nPathSamples, nPathSamplesExtra, cSize, patchCostDistortion, qual, rd, false); //first
				if (doTwice) //second processing is optional (makes sense but did not increase the number of matches as much as i expected)
					c2fHelper(matches1[i], matches2[i], radius, true, nPathSamples, nPathSamplesExtra, cSize, patchCostDistortion, qual, rd, true); //second (did twice)
			} //each patch is processed and the winner binID (match) and its tmp (distortion) is saved for this current level (namings suck)
cout << matches1.size() << " --> ";
			vector< int > matches1Cpy = matches1, matches2Cpy = matches2; //keep the copies before clearing them as i'll use the matches here to rematch the ones that are unmatched via matchIdx=-1 below
			//refill matches1 for the next level iteration using the binID values that encode the matches up to this level (naming sucks); set the matchIdx values as well in case this is the last level iteration (other iterations'll reset matchIdx=-1 above)
			matches1.clear(), matches2.clear();
			if (mesh1->samples.size() != mesh2->samples.size())
				cout << "WARNING: FPSroi() must have created same-size mesh1/2 samples\n\n";
//			for (size_t s = 0; s < mesh2->samples.size(); s++)
//				mesh2->verts[ mesh2->samples[s] ]->matchIdx = -1, mesh2->verts[ mesh2->samples[s] ]->distortion = INF; //mesh2 side must also be reset 'cos 2+ (manyTo1) matches on the same mesh2 sample implies incompatible sampling and hence keep the best of 2+ matches based on distortion
			for (size_t v = 0; v < mesh1->verts.size(); v++) mesh1->verts[v]->matchIdx = -1; for (size_t v = 0; v < mesh2->verts.size(); v++) mesh2->verts[v]->matchIdx = -1, mesh2->verts[v]->distortion = INF; //reset over verts 'cos sample set changes during denseMapLargeN() so be safe (no such problem in denseMap() so reset over samples ok for denseMap())
			for (size_t s = 0; s < mesh1->samples.size(); s++) //same as mesh2->samples.size()
			{
				//mesh1->verts[ mesh1->samples[s] ]->matchIdx = -1; done above via reset over verts
				int v2 = mesh1->verts[ mesh1->samples[s] ]->binID;
//cout << endl << mesh1->samples[s] << " ssssssssssss " << v2 << ":\n";
				if (v2 != -1)
				{
					//keep only the best/minDistortion of 2+ matches targeting the same mesh2 sample
					if (mesh2->verts[v2]->matchIdx == -1) //save the distortion of the match hitting v2 in case another one hits it later (manyTo1); in such a collision keep the best/minDistortion match and cancel the other
					{
//cout << "in case later: " << mesh1->samples[s] << " - " << v2 << "\t" << mesh1->verts[ mesh1->samples[s] ]->distortion << "; "  << matches1.size() << " " << mesh2->verts[v2]->matchIdx << endl;
						mesh2->verts[v2]->distortion = mesh1->verts[ mesh1->samples[s] ]->distortion, matches1.push_back(mesh1->samples[s]), matches2.push_back(v2), mesh1->verts[ mesh1->samples[s] ]->matchIdx = v2, mesh2->verts[v2]->matchIdx = mesh1->samples[s]; //dual operation
					}
					else if (mesh2->verts[v2]->distortion > mesh1->verts[ mesh1->samples[s] ]->distortion) //existing match hitting v2 is worse/highDistortion than the current manyTo1 match so cancel the existing match and set the new match
					{
						//typically tmp is same as distortion but rarely/validly tmp < distortion may arise when mesh1 sample already has a good/small tmp value from a previous overlapping patch; use distortion here to discard that previous (currently irrelevant) patch
						int cancel1 = mesh2->verts[v2]->matchIdx; //matched mesh1 vertex whose match is canceled below
//cout << "cancel: " << cancel1 << " - " << v2 << " w/ disto: " << mesh1->verts[cancel1]->tmp << " == " << mesh1->verts[cancel1]->distortion << " > " << mesh1->verts[ mesh1->samples[s] ]->distortion << " hence new match: " << mesh1->samples[s] << " - " << v2 << "\n\n"; //new match has the better/smaller distortion of mesh1->samples[s]
						mesh1->verts[cancel1]->matchIdx = -1, mesh2->verts[v2]->matchIdx = -1; //cancel by resetting matchIdx values on both sides
						for (size_t c = 0; c < matches1.size(); c++) //cancel by erasing canceled entries from both sides
							if (matches1[c] == cancel1)
								matches1.erase(matches1.begin() + c), matches2.erase(matches2.begin() + c), cancel1 = -2; //erase(begin() + c) erases matches1/2[c]
						if (cancel1 != -2)
							cout << "WARNING: entry not found in matches1/2\n\n";
						//canceled the existing match and now set the new match
						mesh2->verts[v2]->distortion = mesh1->verts[ mesh1->samples[s] ]->distortion, matches1.push_back(mesh1->samples[s]), matches2.push_back(v2), mesh1->verts[ mesh1->samples[s] ]->matchIdx = v2, mesh2->verts[v2]->matchIdx = mesh1->samples[s]; //dual operation
					}
				}
//cout << v2 << "\t\t" << mesh1->samples[s] << " - " << mesh1->verts[ mesh1->samples[s] ]->matchIdx << " set\n";
			} //end of s
			//rematch the ones (in 1to1 manner) that are unmatched above
			for (size_t i = 0; i < matches1Cpy.size(); i++)
				if (mesh1->verts[ matches1Cpy[i] ]->matchIdx == -1 && mesh2->verts[ matches2Cpy[i] ]->matchIdx == -1) //2nd condition prevents manyTo1 assignments
					mesh1->verts[ matches1Cpy[i] ]->matchIdx = matches2Cpy[i], mesh2->verts[ matches2Cpy[i] ]->matchIdx = matches1Cpy[i], //cout << "rematch: " << matches1Cpy[i] << " - " << matches2Cpy[i] << "\n\n";
					matches1.push_back(matches1Cpy[i]), matches2.push_back(matches2Cpy[i]); //matches2Cpy[i].distortion = matches1Cpy[i].distortion not made here 'cos matches1Cpy[i] is unknown; no problem 'cos trimMap() resets mesh1 distortions and mesh2 distortions are unimportant at this point
//			if (trim)
//			{
				qual = trimMap(/*2.0 1.6 1.2*/1.2, mesh1->samples, noTrimCoarseMapMatches, patchCostDistortion); //trimmedSamples[] is updated so it's ready for the cancellations below
				for (size_t s = 0; s < trimmedSamples.size(); s++) //cancel by erasing canceled entries from both sides
					for (size_t c = 0; c < matches1.size(); c++) //cancel by erasing canceled entries from both sides
						if (matches1[c] == trimmedSamples[s])
							matches1.erase(matches1.begin() + c), matches2.erase(matches2.begin() + c); //erase(begin() + c) erases matches1/2[c]
				for (size_t r = 0; r < rd.size(); r++) //trimMap() changed v.distortion values to unnormalized, and possibly irrelevant (trimType), ones; restore the normalized values for the upcoming tmp = distortion above
					mesh1->verts[ rd[r] ]->distortion = mesh1->verts[ rd[r] ]->geoToRef;
//			}
			cout << matches1.size() << " matches (" << (qual == 1 ? "puregeo" : "pointwisepairwise") << " disto: " << distortion << ") after level " << level << " of radius " << radius << "\n\n"; //global distortion is updated just for the print here	
			radius /= 2; //smaller radius (after more levels) can always find nSamplesPerSubregion samples to match but those findings will not respect radius at all, i.e., high nOffRadius value, and will be less meaningful (no pool/gc stuff) so don't decrease radius at all (verified that this leads to bad results so decrease radius)
			if (backToOptimal) //update optimalMap, i.e., my final result may be from level 2 despite nMaxLevels=5
			{
				minDisto = (matches1.size() > (optimalMap.size() / 2) ? INF : minDisto); //start optimalMap search again (via =INF) if the current number of matches is bigger than optimalMap size, i.e., get the densest map
				if (distortion < minDisto) //preferring the densest map is also fair 'cos distortion is not an avg it is a total so sparser maps have unfairly low distortion
				{
					minDisto = distortion, optimalMap.clear();
					for (size_t i = 0; i < matches1.size(); i++) //transfer matches1/2 to optimalMap
						optimalMap.push_back(matches1[i]), optimalMap.push_back(matches2[i]); //2 consecutive entries make 1 match
				}
				if (matches1.size() >= nMaxMatches)
				{
					minDistoDensest = (matches1.size() > nMaxMatches ? INF : minDistoDensest); //start densestOptMap search again (via =INF) if the current number of matches is bigger than nMaxMatches so far, i.e., get the densest map
					nMaxMatches = matches1.size();
					if (distortion < minDistoDensest)
					{
						minDistoDensest = distortion, densestOptMap.clear();
						for (size_t i = 0; i < matches1.size(); i++) //transfer matches1/2 to densestOptMap
							densestOptMap.push_back(matches1[i]), densestOptMap.push_back(matches2[i]); //2 consecutive entries make 1 match
					}
				}
			}
//if(level==0)break;
		} //end of level
		if (backToOptimal) //optimalMap to matchIdxs
		{
			for (size_t s = 0; s < mesh1->samples.size(); s++)
				mesh1->verts[ mesh1->samples[s] ]->matchIdx = mesh2->verts[ mesh2->samples[s] ]->matchIdx = -1; //optimalMap below will reset these values below
			if (minDisto < minDistoDensest)
				for (size_t i = 0; i < optimalMap.size(); i += 2)
					mesh1->verts[ optimalMap[i] ]->matchIdx = optimalMap[i + 1], mesh2->verts[ optimalMap[i + 1] ]->matchIdx = optimalMap[i]; //dual operation
			else //unlikely 'cos more matches mean more additions to distortion (no avg) but still rarely observed minDistoDensest <= minDisto
				for (size_t i = 0; i < densestOptMap.size(); i += 2)
					mesh1->verts[ densestOptMap[i] ]->matchIdx = densestOptMap[i + 1], mesh2->verts[ densestOptMap[i + 1] ]->matchIdx = densestOptMap[i]; //dual operation
			cout << "back to optimal distortion: "; if (patchCostDistortion && dissimType != PAIRWISENOPATCH) mapQualityPatch(); else mapQualityGeo(); cout << distortion << endl; //distortion is updated by the mapQualityPatch/Geo() here
		}
//		matchUnmatcheds(nPathSamples, nPathSamplesExtra); //7 out of 30 remained unmatched in nMaxLevels=5; make 7x7 map and trim bad ones, e.g., match additional 3x3 good subset from the unmatched list
		
		//prevSamples1 = mesh1->samples, prevSamples2 = mesh2->samples; prevSamples1/2 prevents using the current samples as new samples in the next iteration (append below is a better idea)
		prevSamples1.insert(prevSamples1.end(), mesh1->samples.begin(), mesh1->samples.end()), prevSamples2.insert(prevSamples2.end(), mesh2->samples.begin(), mesh2->samples.end()); //append samples1/2 to prevSamples1/2 to prevent using all prev iters' samples as new samples
		nMatches = 0; //break condition is about number of matches which is related to the desired nDenseSamples
		for (size_t v = 0; v < mesh1->verts.size(); v++)
			if (mesh1->verts[v]->matchIdx != -1)
				mesh1->verts[v]->tmpIdx = mesh1->verts[v]->matchIdx; //matchIdx values will be reset to -1 above in the next iteration but tmpIdx are frozen to indicate all valid matches
		for (size_t v = 0; v < mesh1->verts.size(); v++)
			nMatches += (mesh1->verts[v]->tmpIdx != -1);
		cout << nMatches << " / " << nDenseSamples << " matched overall at iter " << iter++ << "\n\n";
//if(iter==3)break;
	} while (nMatches < nDenseSamples);
	//refill samples[] w/ the verts w/ valid matches; samples[] above were always of size k=30 but now samples[] will be of size ~nDenseSamples
	mesh1->samples.clear(), mesh2->samples.clear();
	for (size_t v = 0; v < mesh2->verts.size(); v++)
		mesh2->verts[v]->sample = false;
	for (size_t v = 0; v < mesh1->verts.size(); v++)
		if (mesh1->verts[v]->tmpIdx != -1)
			mesh1->verts[v]->matchIdx = mesh1->verts[v]->tmpIdx, mesh2->verts[ mesh1->verts[v]->tmpIdx ]->matchIdx = v, //matchIdx update w/ dual operation
			mesh1->verts[v]->sample = true, mesh1->samples.push_back(v), mesh2->verts[ mesh1->verts[v]->matchIdx ]->sample = true, mesh2->samples.push_back(mesh1->verts[v]->matchIdx);
//			cout << v << "\t" << mesh1->verts[v]->geodesic.size() << " " << mesh1->verts[v]->geoNormalized << " " << mesh1->verts[v]->geodesic[10] << "; " << mesh1->verts[v]->coarseMapMatch << endl,
//			cout << mesh1->verts[v]->matchIdx << "\t" << mesh2->verts[mesh1->verts[v]->matchIdx]->geodesic.size() << " " << mesh2->verts[mesh1->verts[v]->matchIdx]->geoNormalized << " " << mesh2->verts[mesh1->verts[v]->matchIdx]->geodesic[10] << "; " << mesh2->verts[mesh1->verts[v]->matchIdx]->coarseMapMatch << "\n\n";
		else
			mesh1->verts[v]->sample = false;

	//evaluation of C2FQAP
	mesh1->unitGeoScaling(mesh1->maxGeoDist, false), mesh2->unitGeoScaling(mesh1->maxGeoDist, false);
	char str[10], orgDtype[250]; sprintf(orgDtype, dType); //same as strcpy(orgDtype, dType) or sprintf(orgDtype, "%s", dType);
	sprintf(str, " -qap- "), sprintf(orgDtype, dType), sprintf(dType, "puregeo"), nMatches = mapQualityGeo(), printMap(t, str), sprintf(dType, orgDtype); //mapQualityGeo() and hence puregeo distortions will be printed so update dType first and then roll back to original (char str[] = " -qap- "; same effect)
	cout << "one " << nSamplesPerSubregion << " x " << nSamplesPerSubregion  << " QAP solution took " << (((float) clock()-t) / CLOCKS_PER_SEC) / nQAPs << " secs\n";
	for (size_t v = 0; v < mesh1->verts.size(); v++)
		mesh1->verts[v]->marked = false; //true means display it w/ a sphere
	bool printOnce = true;
	for (size_t i = 0; i < mesh1->coarseMap.size(); i += 2) //display coarseMap matches
	{
		mesh1->verts[ mesh1->coarseMap[i] ]->marked = true;
		if (mesh1->verts[ mesh1->coarseMap[i] ]->matchIdx != mesh1->coarseMap[i+1] && printOnce)
			cout << "original coarseMap changed during dense mapping: " << mesh1->coarseMap[i] << " - " << mesh1->coarseMap[i+1] << " --> " << mesh1->verts[ mesh1->coarseMap[i] ]->matchIdx << " (there may be others; printing stopped)\n\n", printOnce = false;
	}
	/*if (dissimType != PAIRWISENOPATCH) disabled 'cos mapQualityPatch().pathSamples[] are already computed for k=30 samples (through pathSampling() above) but now i have nDenseSamples=250 samples so it's inconsistent; hence no evaluation here (recall samples on the path going from the ith sample to the jth sample are given by samples[i].pathSamples[j])
		mapQualityPatch(), cout << "distortion=" << dType << ": " << distortion << "\n"; //else mapQualityGeo(); redundant 'cos mapQualityGeo/printMap() is called above to learn nMatches and print geodesic distortions (mapQualityVirtual() uses the same metrics so call either of them here)
	else
		cout << "WARNING: Q is filled w/ lengths to coarseMap but mapQualityGeo() is not so it's OK to get different objective & distortion prints\n\n";*/
}


void Correspondence::fullMap(bool roiReady, bool forceSample)
{
	//computes the fully dense map b/w mesh1 and mesh2 vertices using feature vector comparisons (if mesh2 is partial then some mesh1 vertices will remain unmatched as desired/expected)

	vector< int > anchors1, anchors2; //landmark/anchor matches obtained in the qap part, i.e., after denseMap()
	int nMatches = 0; //interpolate nMatches matches below
	for (size_t b = 0; b < mesh1->samples.size(); b++)		
		if (mesh1->verts[ mesh1->samples[b] ]->matchIdx != -1)
			nMatches++;
	if (nMatches == 0)
		cout << "WARNING: nothing to interpolate\n\n";
	distortion /= nMatches; //make it compatible w/ b.distortions below by taking the avg of the total sum in distortion
	for (size_t b = 0; b < mesh1->samples.size(); b++)
		if (mesh1->verts[ mesh1->samples[b] ]->matchIdx != -1)
			if (mesh1->verts[ mesh1->samples[b] ]->distortion < distortion * 1.1)
				anchors1.push_back(mesh1->samples[b]), anchors2.push_back(mesh1->verts[ mesh1->samples[b] ]->matchIdx);
	clock_t t = clock();
	cout << "\n(optional) full matching via " << anchors1.size() << " landmarks..";
	bool oneToOne = false; //true: ensure 1-to-1 matching by preserving for each mesh2.v the match w/ the closest descriptor distance (may leave some mesh1 verts unmatched)
	int nCloseAnchors = -10; //use distance to this many close robust matches in order to eliminate tooFar/erroneous robusts in the construction of v.geoAll
	mesh1->geoAll(anchors1, nCloseAnchors), mesh2->geoAll(anchors2, nCloseAnchors);
	//clear existing matches so feature-based matching below starts from fresh (i might keep the anchor matches but they'll be naturally preserved anyway thanks to geoAll construction)
	for (size_t i = 0; i < mesh1->verts.size(); i++)
		mesh1->verts[i]->matchIdx = -1;
	for (size_t i = 0; i < mesh2->verts.size(); i++)
	{
		mesh2->verts[i]->matchIdx = -1;
		mesh2->verts[i]->tmp = INF; //preserve for each mesh2.v the match w/ the closest descriptor distance tmp (naming sucks)
	}
	for (size_t v1 = 0; v1 < mesh1->verts.size(); v1++)
	{
		if (! roiReady || (roiReady && mesh1->verts[v1]->roi)) //process only the ROI vertices; if roi not available, i.e., NODENSEMAP mode, then all verts pass this condition trivially
		{
			double minDescDiff = INF, val;
			for (size_t v2 = 0; v2 < mesh2->verts.size(); v2++)
				if (! roiReady || (roiReady && mesh2->verts[v2]->roi)) //process only the ROI vertices; if roi not available, i.e., NODENSEMAP mode, then all verts pass this condition trivially
				{
					//match like many-to-one (may give many matches to one mesh2 vert) and then convert to one-to-one if requested via oneToOne=true
					double maxVal = -INF, currDescDiff = 0.0; //for the match i-j
					for (size_t d = 0; d < mesh1->verts[v1]->geoAll.size(); d++) //for each dimension of the descriptor vector
					{
						val = abs(mesh1->verts[v1]->geoAll[d] - mesh2->verts[v2]->geoAll[d]); //geoAll might have represented diffusion/biharmonic distances or gps coordinates etc as well
						currDescDiff += (val * val); //L2-norm
					}
					currDescDiff = sqrt(currDescDiff); //L2-norm
					if (currDescDiff < minDescDiff)
					{
						minDescDiff = currDescDiff;
						mesh1->verts[v1]->matchIdx = v2;
						mesh1->verts[v1]->tmp = minDescDiff; //score of this match
					}
				} //end of v2			
		}
		if (v1 % 10000 == 0) //show progress
			cout << v1 << "/" << mesh1->verts.size() << " ";
	}
	int nm = 0;//, nTrimmed = 0;
	//set mesh2.matchIdx/tmp based on the best matches of mesh1.verts set above
	for (size_t v1 = 0; v1 < mesh1->verts.size(); v1++)
	{
		int v22 = mesh1->verts[v1]->matchIdx;
		if (mesh1->verts[v1]->matchIdx != -1)
		{
			if (mesh1->verts[v1]->tmp < mesh2->verts[v22]->tmp) //current match v1-v22 is better than the existing one another1-v22 (less dissimilarity); so update mesh2 vert's dissim and matchIdx
			{
				mesh2->verts[v22]->tmp = mesh1->verts[v1]->tmp;
				mesh2->verts[v22]->matchIdx = v1;
			}
			nm++;
		}
	}
	cout << endl << nm << " matched manyTo1" << (oneToOne ? " (not in use) " : " "), nm=0;
	if (oneToOne) //many-to-one case may give many matches to one mesh2 vert; in one-to-one case remove all but the best of these matches hence giving 1 (best/recovered) match to that mesh2 vert
	{
		for (size_t i = 0; i < mesh1->verts.size(); i++)
			mesh1->verts[i]->matchIdx = -1; //remove all matches
		for (size_t v2 = 0; v2 < mesh2->verts.size(); v2++)
		{
			int v1 = mesh2->verts[v2]->matchIdx; //best match for this v2 (set above this oneToOne block) is v1
			if (v1 != -1)
				mesh1->verts[v1]->matchIdx = v2; //recover the best match
		}
		for (size_t i = 0; i < mesh1->verts.size(); i++)
			nm += (mesh1->verts[i]->matchIdx != -1);
		cout << nm << " matched 1to1 in " << ((float) clock()-t) / CLOCKS_PER_SEC;
	}
	else
		cout << "in " << ((float) clock()-t) / CLOCKS_PER_SEC;
	if (forceSample) //force mesh1.samples to match w/ mesh2.samples by changing their current match (that may be a non-sample) w/ the closest non-sample; changed/new matches guaranteed to be sample only if oneToOne=false below
	{
		oneToOne = true;//false; //forced samples must respect 1-to-1 matching (update existing variable oneToOne 'cos this block is different and 1-to-1 is in fact encouraged here)
		if (oneToOne) //marked heuristic may violate 1-to-1 property but it is worth it; simply use 'oneToOne && mesh2->verts[ mesh2->samples[j] ]->matchIdx != -1' condition below instead to disable this heuristic and get exact 1-to-1
			for (size_t j = 0; j < mesh2->samples.size(); j++)
				if (mesh2->verts[ mesh2->samples[j] ]->matchIdx != -1 && mesh1->verts[ mesh2->verts[ mesh2->samples[j] ]->matchIdx ]->sample)
					mesh2->verts[ mesh2->samples[j] ]->marked = true; //matched and the match is a sample so cannot be the winner closestSample below; note that this old match and a potential closestSample match below will violate 1-to-1 but prevents selecting a ~bad (but available) sample within radius
				else
					mesh2->verts[ mesh2->samples[j] ]->marked = false; //unmatched or matched by a non-sample so can be the winner closestSample below (true means unavailable to be a closestSample below)
		for (size_t i = 0; i < mesh1->samples.size(); i++)
			if (mesh1->verts[ mesh1->samples[i] ]->matchIdx != -1 && ! mesh2->verts[ mesh1->verts[ mesh1->samples[i] ]->matchIdx ]->sample) //a non-sample match for this mesh1.sample
			{
				int v2 = mesh1->verts[ mesh1->samples[i] ]->matchIdx, closestSample = -1; //v2 is an roi vertex for sure in roiReady=true mode; closestSample below will be an roi vertex as it is sampled inside roi (via FPSroi()) for sure
				double minDist = INF;
				for (size_t j = 0; j < mesh2->samples.size(); j++) //find the closest mesh2 sample to non-sample v2
					if (mesh2->verts[ mesh2->samples[j] ]->geodesic[v2] < minDist)
					{
						//if (oneToOne && mesh2->verts[ mesh2->samples[j] ]->matchIdx != -1) //1-to-1 and j is already matched so cannot be the winner closestSample now (enabled here and disable below for an exact 1-to-1 behavior)
						if (oneToOne && mesh2->verts[ mesh2->samples[j] ]->marked) //1-to-1 may be violated but better heuristic that prevents j only if it's matched to a sample, i.e., it's ok if it's matched to a non-sample (see marked = false comment above)
							continue;
						if (oneToOne)
						{
							if (mesh2->verts[ mesh2->samples[j] ]->geoNormalized) //may be true (verified) so obtain unnormalized geodesics via a new dijkstra
								mesh2->verts[ mesh2->samples[j] ]->geodesic.clear(), mesh2->dijkstraShortestPaths(mesh2->samples[j]);//cout << "unnormalized geo obtained\n";
							if (mesh2->verts[ mesh2->samples[j] ]->geodesic[v2] > 2 * mesh1->radius) //j is too far away so cannot be the winner closestSample now (same as mesh2.radius)
								continue;
						}
						minDist = mesh2->verts[ mesh2->samples[j] ]->geodesic[v2], closestSample = mesh2->samples[j]; //many-to-1 is possible if oneToOne=false
					}
				if (closestSample != -1) //always true in manyTo1=true mode
					mesh1->verts[ mesh1->samples[i] ]->matchIdx = closestSample, mesh2->verts[closestSample]->matchIdx = mesh1->samples[i], mesh2->verts[closestSample]->marked = true; //dual operation
			}
		if (oneToOne) //undo dijkstra's above to be ready for the upcoming mapQualityGT() that requires normalized geodesics
			mesh2->unitGeoScaling(g, false); //use g instead of mesh1->maxGeoDist which may be 1 due to a previous unitGeoScaling() call
		cout << " secs\n";
	}
	else
		cout << " secs\n";
}

void Correspondence::transferColors(bool roiReady)
{
	//transfers predefined mesh2 colors to mesh1 using the matches; if multiple mesh1 verts are matched to a single mesh2 vert, i color all those verts in mesh1 by transferring predefined mesh2 colors to mesh1

	//geodesic distances to the reference vertex for visualization only (to transfer smooth colors from mesh2 to mesh1)
	for (size_t j = 0; j < mesh2->verts.size(); j++)
	{
		size_t v = 3643;//0; //reference geo/seo source vertex
		v = (v < mesh2->verts.size() ? v : 0); //stay in bounds
		mesh2->dijkstraShortestPaths(v); //geodesic[] ready (returns immediately if v.geodesic is already not empty)
		mesh2->verts[j]->geoToRef = mesh2->verts[v]->geodesic[j]; //make sure vertColors() below uses geoToRef too
	}
	mesh2->vertColors(0, roiReady); //predefined mesh2 colors are set
	bool print1 = true;
	int nMatched = 0;
	for (size_t i = 0; i < mesh1->verts.size(); i++) //all or some (some in partial matching case) mesh1 vert has a match and it'll use that matchIdx to learn its color from the fully colored mesh2
	{
		mesh1->verts[i]->color = new float[3]; //alloc memo for v.color and init it to mesh color (gray)
		mesh1->verts[i]->color[0] = mesh1->color[0]; //unmatched mesh2 verts will have this default color (gray set in Mesh's constructor)
		mesh1->verts[i]->color[1] = mesh1->color[1];
		mesh1->verts[i]->color[2] = mesh1->color[2];
		if (mesh1->verts[i]->matchIdx != -1)
			mesh1->verts[i]->color = mesh2->verts[ mesh1->verts[i]->matchIdx ]->color, nMatched++;
		else if (print1)
			cout << "WARNING: unmatched mesh1 vert (there may be other unmatched ones; printing stopped)\n\n", print1 = false;
	}
	cout << nMatched << " colors pulled back from mesh2\n";	
}

void Correspondence::bestWorstMatches()
{
	//sets the best and worst matches according to current (or below) distortion values

/*	if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
		mapQualityPatch();
	else//
		mapQualityGeo();//*/
	double minDisto = INF, maxDisto = -INF;
	for (size_t i = 0; i < mesh1->samples.size(); i++)
	{
		if (mesh1->verts[ mesh1->samples[i] ]->matchIdx != -1 && mesh1->verts[ mesh1->samples[i] ]->distortion < minDisto)
			minDisto = mesh1->verts[ mesh1->samples[i] ]->distortion, best = mesh1->samples[i];
		if (mesh1->verts[ mesh1->samples[i] ]->matchIdx != -1 && mesh1->verts[ mesh1->samples[i] ]->distortion > maxDisto)
			maxDisto = mesh1->verts[ mesh1->samples[i] ]->distortion, worst = mesh1->samples[i];
	}
	cout << "worst (by current distortion): " << worst << "-" << mesh1->verts[worst]->matchIdx << " " << mesh1->verts[worst]->distortion << endl;
	cout << "best (by current distortion): " << best << "-" << mesh1->verts[best]->matchIdx << " " << mesh1->verts[best]->distortion << endl;
}

void Correspondence::matchUnmatcheds(int nPathSamples, int nPathSamplesExtra)
{
	//matches the unmatched samples after C2FQAP and trims unsafe ones, e.g. 7 of 30 may be unmatched and qap below finds 7x7 (no c2f) and then trimming returns the 3x3 safe subset (if 4 matches are trimmed)

	mesh1->samplesC2F.clear(), mesh2->samplesC2F.clear(); //put the unmatcheds inside samplesC2F (naming sucks) as fillQdenseC2F() below uses samplesC2F[]
	vector< int > matched1; //used in trimByCurrentQAP=false mode only
	for (size_t i = 0; i < mesh1->samples.size(); i++)
		if (mesh1->verts[ mesh1->samples[i] ]->matchIdx == -1)
			mesh1->samplesC2F.push_back(mesh1->samples[i]);//cout << mesh1->samples[i] << " "; cout << "unmatched1\n";
		else
			matched1.push_back(mesh1->samples[i]);
	for (size_t i = 0; i < mesh2->samples.size(); i++)
		if (mesh2->verts[ mesh2->samples[i] ]->matchIdx == -1)
			mesh2->samplesC2F.push_back(mesh2->samples[i]);//cout << mesh2->samples[i] << " "; cout << "unmatched2\n";
	if (mesh1->samplesC2F.size() != mesh2->samplesC2F.size())
		cout << "WARNING: bad unmatched size\n\n"; //sizes must have been equal and each size should be <50 to make qap feasible	
	if (mesh1->samplesC2F.empty() || mesh1->samplesC2F.size() > 40 || dissimType != POINTWISEPAIRWISE) //PAIRWISENOPATCH may have -ve nPathSamplesExtra here which may cause problems; not that important so just return early via 3rd condition
		return; //nothing to match so return (no problem if this occurs) or too big to match so memory problems arise so return
	else if (mesh1->samplesC2F.size() > 20)
		cout << "WARNING: prepare to wait as 27x27 (37x37) QAP takes 343 (3017) secs w/ the fastest QuCM solver (current QAP is " << mesh1->samplesC2F.size() << " x " << mesh1->samplesC2F.size() << ")\n\n";
	clock_t t = clock();
	bool trimByCurrentQAP = false; //true: trim using the current qap map b/w newly matched 7x7 samples; false: use the previous 23x23 map which is supposed to be safer (as the matches there have survived the previous trimmings)
	int nMatches = mapQualityGeo(); //global distortion ready
	double safeDistortion = distortion / nMatches, n = 0.0;
	//////////////////////////// copied from denseMap() ////////////////////////////
	cout << mesh1->samplesC2F.size() << " x " << mesh2->samplesC2F.size() << " qap to match the unmatcheds..\n";
	size_t cSize = mesh1->coarseMap.size() / 2; //number of matches in coarseMap[] (2 consecutive entries make 1 match)
	//define the patches by sampling nSamplesPerSubregion samples inside them (via radius restriction and v.roi restriction where the latter is implied as i select a subset of samples from samples[] which are definitely v.roi=true)
	if (! mesh1->extraSamples) //nPathSamples=-1: means special flag for Q2Ddense usage (diag: lengthsToCoarseMap off-diag: lengthsBetweenSamples); nPathSamples>0: means diag: ringAreas off-diag: patchAreasBetweenSamples (alternative not in use: diag: patchAreasToCoarseMap off-diag: patchAreasBetweenSamples)
		mesh1->fillQdenseC2F(nPathSamples + 1, areaWeight), mesh2->fillQdenseC2F(nPathSamples + 1, areaWeight), //QAP needs the Q matrix
		qapMap(nPathSamples + 2 + cSize, true, true); //nPathSamples=-1 sends 1 which is a special flag for Q2D usage; otherwise nPathSamples + 1 gives # path samples and another + 1 gives access to the additional dimension storing pure geodesic length info (for maxGeoDiff filtering or geoDiff-based barrier) and another + cSize for lengths to coarseMap entries
	else //nPathSamples!=-1 for sure 'cos nPathSamples=-1 means no path sampling requested, i.e., POINTWISE or PAIRWISENOPATCH, and in that case extraSamples cannot be true
		mesh1->fillQdenseC2F(nPathSamplesExtra, areaWeight), mesh2->fillQdenseC2F(nPathSamplesExtra, areaWeight), //+1 added during declaration so no nPathSamplesExtra + 1 here
		qapMap(nPathSamplesExtra + 1 + cSize, true, true); //+ 1 gives access to the additional dimension storing pure geodesic length info (for maxGeoDiff filtering or geoDiff-based barrier) and another + cSize for lengths to coarseMap entries
	if (trimByCurrentQAP)
		trimMap(1.0, mesh1->samplesC2F, true, true), n = trimmedSamples.size(); //be more aggressive than denseMap() via 1.0 'cos possibly inconsistent unmatcheds are matched in this function
	else
	{
		//assign distortion to each new match by traversing them over the previous safer 23x23 map
		double avgDisto = 0.0;
		for (size_t i = 0; i < mesh1->samplesC2F.size(); i++)
		{
			int v1 = mesh1->samplesC2F[i], v2 = mesh1->verts[ mesh1->samplesC2F[i] ]->matchIdx;
			mesh1->verts[v1]->distortion = 0.0;
			for (size_t t = 0; t < matched1.size(); t++) //traverse v1-v2 match over all matches1 to get its individual distortion
				mesh1->verts[v1]->distortion += abs( mesh1->verts[v1]->geodesic[ matched1[t] ] - mesh2->verts[v2]->geodesic[ mesh1->verts[ matched1[t] ]->matchIdx ] );
			avgDisto += mesh1->verts[v1]->distortion;
		}
		avgDisto /= mesh1->samplesC2F.size(); //use this avg based on current matches that traversed safe map or safeDistortion avg that is direct distortion of the safe map; typically avgDisto > safeDistortion, e.g., 1.12957 > 0.574992
//cout << avgDisto << " > " << safeDistortion << endl; 1.12957 > 0.574992 typically
		//trimming
		for (size_t s = 0; s < mesh1->samplesC2F.size(); s++)
		{
//			cout << mesh1->samplesC2F[s] << " -puregeoTrim?- " << mesh1->verts[ mesh1->samplesC2F[s] ]->matchIdx << "\t" << mesh1->verts[ mesh1->samplesC2F[s] ]->distortion << " vs. avg:" << safeDistortion << ": ";
			if (mesh1->verts[ mesh1->samplesC2F[s] ]->distortion > safeDistortion)//avgDisto) //be more aggressive than denseMap() via 1.0 'cos possibly inconsistent unmatcheds are matched in this function
				mesh2->verts[ mesh1->verts[ mesh1->samplesC2F[s] ]->matchIdx ]->matchIdx = -1, mesh1->verts[ mesh1->samplesC2F[s] ]->matchIdx = -1, n++;//cout << "yes\n"; //else cout << "no\n";
		}
//		if (n > 0) //no mapQualityGeo() 'cos it'll be called on return and then puregeo disto will be printed; see the nMatches = mapQualityGeo(), printMap(t, str) part
//			cout << n << " match" << (n > 1 ? "es" : "") << (n > 1 ? " are" : " is") << " trimmed\n";
	}//disable trim (put comment opener before if(trimByCurrentQAP)) to see all samples for visualization only*/
	if (mesh1->samplesC2F.size() > n)
		cout << mesh1->samplesC2F.size() - n << " unmatcheds are matched in " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs\n";
	else
		cout << "all trimmed so no new match after " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs\n";
	//////////////////////////// copied from denseMap() ends ////////////////////////////
}

void Correspondence::dissimTypeString()
{
	//sets dType which stores the name of the current dissimilarity type as a string

	sprintf(dType, "%s", (dissimType == PAIRWISE ? "pairwise" : "pointwisepairwise")); //for printing only
	if (dissimType == POINTWISE)
		sprintf(dType, "pointwise");
	else if (dissimType == PAIRWISENOPATCH)
		sprintf(dType, "puregeo");
}

double Correspondence::distortionByPermutation(const int* perm, const int& n, const int& m)
{
	//return the distortion of the map implied by the permutation perm

	for (int j = 0; j < n; ++j)
		mesh1->verts[ mesh1->samples[j] ]->matchIdx = -1; //reset matchIdx values that have been updated during the previous call to this function
	vector< int > tmpMap; //2 consecutive entries (mesh1-mesh2) gives 1 match of the entire map/correspondence
	int idx1 = 0, idx2 = 0;
	for (int y = 0; y < 2*m; y++)
		if (y % 2 == 0) //even idxs for permutations from mesh1 side
			tmpMap.push_back(mesh1->samples[ perm[idx1] ]);
		else			//odd idxs from the fixed mesh2 side w/ virtual samples (n-m indices not covered by the current permutation are assumed to matched to virtual mesh2 samples)
		{
			if (idx2 < m) //fixed 0,1,2,..,m part
			{
				tmpMap.push_back(mesh2->samples[idx2]);
				mesh1->verts[ mesh1->samples[ perm[idx1] ] ]->matchIdx = mesh2->samples[idx2]; //mapQuality() needs matchIdx values
			}
			else          //virtual sample part
			{
				tmpMap.push_back(-2);
				mesh1->verts[ mesh1->samples[ perm[idx1] ] ]->matchIdx = -2; //special flag for virtual samples
			}
			idx1++, idx2++;
		}
	mapQualityVirtual(); //global distortion is updated
	return distortion;
}
void Correspondence::qapMap2(int nMaxIters, int nSweeps, double beta_min, double beta_max)
{
	//same as qapMap() except this one uses pure simulated annealing (no qucoop stuff) (adapted from simulatedAnnealing())

	cout << "QAP map computation via simulated annealing..\n";

//	std::mt19937 rng(42); //fixed seed (10 calls to realDist(rng) always generates 0.796543 0.183435 0.779691 0.59685 0.445833 0.0999749 0.459249 0.333709 0.142867 0.650888; similarly 7 calls to rng() % 2 always generates 0 1 0 0 0 1 0)
	std::mt19937 rng( (unsigned) time(NULL) ); //non-fixed seed that changes with every execution (recommended)
	std::uniform_real_distribution<double> realDist(0.0, 1.0);

	//m-size permutations from mesh1 side will match w/ the fixed m-size mesh2 side (n-m indices not covered by the current permutation are assumed to matched to virtual mesh2 samples)
	int n = mesh1->samples.size(), m = mesh2->samples.size(), * bestPerm = new int[m], * perm = new int[m], * inUse = new int[n], swap1, swap2, tmp;
	double E, newE, delta, bestE = INF, finalE;
	for (int i = 0; i < nMaxIters; i++) //multistart: each of nMaxIters iterations initialized by a random perm array tries to update bestPerm below
	{
		for (int j = 0; j < n; ++j)
			inUse[j] = 0; //prevent duplications in perm[]; inUse[5]=1 means number 5 is inside perm[]
		for (int j = 0; j < m; ++j)
		{
			do
			{
				perm[j] = rng() % n; //random initialization w/ an integer in [0, n)
			} while (inUse[ perm[j] ] == 1);
			inUse[ perm[j] ] = 1;
		}
		
		E = distortionByPermutation(perm, n, m);

		for (int sweep = 0; sweep < nSweeps; ++sweep) //cooling-down: beta-based
		{
			double frac = double(sweep) / nSweeps;
			double beta = beta_min * pow(beta_max / beta_min, frac);

			//swap (2 matches are affected per swap) and check the new energy (revert in reject state)
			for (int b = 0; b < m; b++)
			{
				swap1 = rng() % m, swap2 = rng() % m;
				while (swap2 == swap1) //swap1 and swap2 must be different indices so they can swap validly
					swap2 = rng() % m;
				tmp = perm[swap1], perm[swap1] = perm[swap2], perm[swap2] = tmp; //permutation is still bijective, i.e., no duplication after swap

				newE = distortionByPermutation(perm, n, m), delta = newE - E;
				if (delta < 0 || realDist(rng) < exp(-beta*delta)) //newE is better/smaller or temparature-based random access is granted (accept state)
					E = newE;
				else
					tmp = perm[swap2], perm[swap2] = perm[swap1], perm[swap1] = tmp; //revert
			}
			if (i % 10 == 0 && sweep % 25 == 0)
				cout << "iterInner: " << i << "\tcost: " << E << "\tbeta: " << beta << endl;
		}
		finalE = distortionByPermutation(perm, n, m); //energy based on the swapped perm
		if (finalE < bestE)
		{
			bestE = finalE; //bestPerm = perm; not ok 'cos a change in perm[] affects bestPerm[]
			for (int j = 0; j < m; j++)
				bestPerm[j] = perm[j];
		}
		if (i % 10 == 0)
			cout << "iterationOuter: " << i << "\tcost: " << finalE << "\tbestCost: " << bestE << "\n\n";
	}
	//bestPerm to matchIdxs
	//..

	delete [] perm;
	delete [] inUse;
	delete [] bestPerm;
	mapOK = true;
}

bool Correspondence::printMap(clock_t t, char* mapType)
{
	//prints map produced by qapMap() or blossomMap()

	bool infPrinted = false;
	int nMatches = 0;
/*	char subStr[250] = "";
	if (strlen(mapType) > 4) //either " -qap- " or " -blo- " so take qap or blo for the printing below (strlen(mapType) is 7 for " -qap- ")
		sprintf(subStr, "%c%c%c", mapType[2], mapType[3], mapType[4]); this info already printed so no subStr now*/
	for (size_t s = 0; s < mesh1->samples.size(); s++)
		if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
		{
			nMatches++;
			cout << mesh1->samples[s] << mapType << mesh1->verts[ mesh1->samples[s] ]->matchIdx << "\t" << mesh1->verts[ mesh1->samples[s] ]->distortion << dType << endl; //subStr
			if (! infPrinted)
				infPrinted = (mesh1->verts[ mesh1->samples[s] ]->distortion == INF);
			if (mesh1->samples.size() > 10 && s > 4)//&& strcmp(mapType, " -blo- ") == 0 || strcmp(mapType, " -qap- ")) //print only a subset of dense blossom or qap matches
				cout << ".. not printing the remaining matches\n", s = mesh1->samples.size();
		}
	cout << ">>>>>>>> " << ((float) clock()-t) / CLOCKS_PER_SEC << " secs for matching samples (distortion=" << dType << ": " << distortion << ", avg: " << distortion / nMatches << ")\n";
	/*//print indices in samples[] arrays
	for (size_t s = 0; s < mesh1->samples.size(); s++)
		if (mesh1->verts[ mesh1->samples[s] ]->matchIdx != -1)
		{
			size_t s2 = 0;
			for (; s2 < mesh2->samples.size(); s2++)
				if (mesh2->samples[s2] == mesh1->verts[ mesh1->samples[s] ]->matchIdx)
					break;
			cout << s << "\t" << s2 << endl;			
		}//*/
	return infPrinted;
}

MatrixXd padMatrix(const MatrixXd& matrix, int pad)
{
	//pads the matrix w/ pad many 0-rows and 0-cols

    int n = matrix.rows();
    MatrixXd padded = MatrixXd::Zero(n + pad, n + pad);
    padded.block(0, 0, n, n) = matrix;
    return padded;
}

double avgUnary = 0.0, nu = 0, avgPairwise = 0.0, np = 0; //for statistic only
MatrixXd computeQ(const MatrixXd& X, const MatrixXd& Y, int p)
{
	//computes the joint n^2 x n^2 Q matrix based on mesh1 and mesh2 content encoded in nxn X and Y separately (last p rows and cols of Y is padded w/ 0s so set related Q entries to 0)

    int n = X.rows(), u, v, q = n - p; //rows and cols after index q (inclusive) are 0-padded so force the related dummy Q entries to be 0
    MatrixXd Q = MatrixXd::Zero(n * n, n * n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				for (int l = 0; l < n; ++l)
				{
					bool pairwise = false;
					u = n * i + k, v = n * j + l; //to match pairs (i, j) vs. (k, l) which also corresponds to the Kronecker product indexing
					if (i == j && k == l)
						Q(u, v) = (k < q ? abs(X(i, i) - Y(k, k)) : 0.0); //unary descriptor l1 (=l2) difference (l2 difference reduces to l1 difference since descriptors are 1-D)
					else if ((i == j && k != l) || (i != j && k == l))
						Q(u, v) = 0.0; //skip 'cos i don't want to take differences of unary and pairwise terms (nonsensical units)
					else
						Q(u, v) = (k < q && l < q ? abs(X(i, j) - Y(k, l)) : 0.0), pairwise = (k < q && l < q); //pairwise geodesic l1 (=l2) difference (l2 difference reduces to l1 difference since geodesics are 1-D)
					//statistics to see how the unary weight (areaWeight) combines unary and pairwise terms (stdDev should be low for a successful combination; use barrier = 0 to see barrier-free compatibility for a better idea)
					if (i == j && k == l && k < q) //unary term that is not padded
						avgUnary += Q(u, v), nu++; //unaryTerms.push_back(Q(u, v));
					else if (pairwise) //pairwise term that is neither padded nor dummy (set this inside else above: pairwise = (k < q && l < q);)
						avgPairwise += Q(u, v), np++; //pairwiseTerms.push_back(Q(u, v));*/
				}
    return Q;
}

bool infCost; //there's at least 1 INF cost in Q due to maxGeoDiff thresholding trick
MatrixXd computeQ2(const vector< MatrixXd >& Xs, const vector< MatrixXd >& Ys, int p, double maxGeoDiff, double barrier, size_t cSize)
{
	//same as computeQ() except accumulates Q(i, j) value using each Xs[dim] and Ys[dim], i.e., have >1 X and Y matrices here (last p rows and cols of Y is padded w/ 0s so set related Q entries to 0)

    int n = Xs[0].rows(), nDimensions = Xs.size(), u, v, q = n - p; //rows and cols after index q (inclusive) are 0-padded so force the related dummy Q entries to be 0
    MatrixXd Q = MatrixXd::Zero(n * n, n * n);
	infCost = false; //becomes true if 1+ Q entry is set to INF	
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				for (int l = 0; l < n; ++l)
				{
					u = n * i + k, v = n * j + l; //to match pairs (i, j) vs. (k, l) which also corresponds to the Kronecker product indexing
					double val = 0.0, geoDiff = 1.0; //from source i to j vs. from target k to l
					vector< double > geoDiffs; //from source j to coarseMap source sides vs. from target l to coarseMap target sides; cSize many entries to be filled below where cSize is the number of matches in coarseMap[]
					bool pairwise = false;
					if (i == j && k == l)
						for (int dim = 0; dim < nDimensions; ++dim) //no dim<nDimensions-1 unlike the else part 'cos geodesic length b/w i & i and k & k is always 0 so no effect on val
							val += (k < q ? pow(Xs[dim](i, i) - Ys[dim](k, k), 2) : 0.0);
					else if ((i == j && k != l) || (i != j && k == l))
						; //skip 'cos i don't want to take differences of unary and pairwise terms (nonsensical units)
					else
					{
						geoDiff = abs(Xs[nDimensions - cSize - 1](i, j) - Ys[nDimensions - cSize - 1](k, l)); //source i to j vs target k to l geodesic distance difference; cSize coarseMap entries and then comes the desired entry hence -cSize-1
						for (size_t e = 0; e < cSize; e++) //extra info that does not appear in paper's pseudocode (Alg. 1); omitted in pseudocode to prevent overcomplication (needed to define coarseMap and cSize then)
							geoDiffs.push_back( abs(Xs[nDimensions - cSize + e](i, j) - Ys[nDimensions - cSize + e](k, l)) ); //source j to coarseMap entry vs target l to coarseMap entry geodesic distance difference
#ifdef GEODIFFTHRESHOLDING //adapted from mapQualityVirtual()
						if (geoDiff > maxGeoDiff) //last entry in Xs and Ys is special, holding the length info; very different length paths imply high overall dissimilarity
							val = INF * INF, infCost = true; //will be sqrted below to get INF that is consistent w/ mapQualityVirtual() (no early return here 'cos nIters may quit the caller loop in which case i need a valid Q (ok despite INF entries in it))
						else
#endif
							//for (int dim = 0; dim < nDimensions - 1; ++dim) //-1 'cos last dimension stores geodesic length b/w i & j (for X) and k & l (for Y) so it's treated separately below during maxGeoDiff filtering (or geoDiff-based barrier)
							for (int dim = 0; dim < nDimensions - (int) cSize - 1; ++dim) //-1 for extra i & j (for X) and k & l (for Y) geoDiff and -cSize for extra j & coarseMap (for X) and l & coarseMap (for Y) geoDiffs
								val += (k < q && l < q ? pow(Xs[dim](i, j) - Ys[dim](k, l), 2) : 0.0);
						pairwise = (k < q && l < q);
					}
					Q(u, v) = sqrt(val); //unary nDimensions-dimensional descriptor l2 difference (i == j && k == l case) or pairwise (nDimensions-1)-dimensional augmented geodesic l2 difference (else case)
#ifndef GEODIFFTHRESHOLDING //geoDiff=1 for unary and skip parts (from declaration); also infCost=false for sure which causes the caller while loop break after the first iteration (so no loop)
					if (barrier < 0) //geoDiff=1 of unary always makes *= 1 effect (as desired) except the tanh case so never use tanh below
						//Q(u, v) *= geoDiff; //linear weighting y=x line where x is geoDiff
						Q(u, v) *= (geoDiff * (geoDiff*geoDiff - 3*geoDiff + 3)); //non-linear cubic weighting that stays above y=x line where x is geoDiff (endpoints 0 to 1)
						//Q(u, v) *= sqrt(geoDiff); //non-linear weighting that stays above y=x line for x < 1 (endpoints 0 to 1)
						//Q(u, v) *= tanh(geoDiff); don't use tanh weight 'cos geoDiff=1 of unary doesn't give *= 1 effect in tanh(1)=0.76 case //non-linear weighting that stays below y=x line for x < 1 (endpoints 0 to 0.76) (very close to 1-e^-x which ends at 0.63)
						//Q(u, v) *= (geoDiff*geoDiff); //non-linear weighting that stays below y=x line for x < 1 (endpoints 0 to 1)
						//Q(u, v) *= 1.0; //no weighting
					//else if (! ((i == j && k == l) || (i == j && k != l) || (i != j && k == l)) && k < q && l < q) //same effect w/ a simpler pairwise boolean below
					else if (pairwise) //geoDiff=1 for unary and skip parts so prevent their nonzero += to Q below (by admitting only pairwise=true terms below); similarly only k<q && l<q are valid entries which is encoded in pairwise=true
					{
						Q(u, v) += (barrier * geoDiff * (geoDiff*geoDiff - 3*geoDiff + 3)); //geoDiff is added to the unweighted pairwise term (geoDiff=1 for unary and skip parts from declaration)
						for (size_t e = 0; e < geoDiffs.size(); e++)
							Q(u, v) += (barrier * geoDiffs[e] * (geoDiffs[e]*geoDiffs[e] - 3*geoDiffs[e] + 3)); //geoDiffs[e] is added to the unweighted pairwise term (geoDiffs is empty for unary and skip parts from declaration)
					}
#endif
					//statistics to see how the unary weight (areaWeight) combines unary and pairwise terms (stdDev should be low for a successful combination; use barrier = 0 to see barrier-free compatibility for a better idea)
					if (i == j && k == l && k < q) //unary term that is not padded
						avgUnary += Q(u, v), nu++; //unaryTerms.push_back(Q(u, v));
					else if (pairwise && Q(u, v) < INF) //pairwise term that is neither padded nor dummy (set this inside else above: pairwise = (k < q && l < q)); 2nd condition ensures that maxGeoDiff-thresholded INF values are not used in statistics
						avgPairwise += Q(u, v), np++; //pairwiseTerms.push_back(Q(u, v));*/
				}
	return Q;
}

void Correspondence::qapMap(int k, bool dense, bool unmatchedCall)
{
	//matches n points on source to m points on target where m <= n by solving a quadratic assignment problem (qap) via simulated annealing based qucoop (k is the dimension of unary and pairwise features that guide the QAP)
	
	if (brute && nSamplesPerSubregion == 5 && ! mesh1->samplesC2F.empty() && ! unmatchedCall) //3rd condition means coarse-to-fine mode, i.e., n=m=30+ typically; 4th condition true if called from matchUnmatcheds() in which case bruteForceMapQuinC2F() cannot be done (as # unmatcheds is not necessarily 5)
	{		
		if (bruteForceMapQuinC2F()) //match 5x5 patches via brute-force combinatorial search
			return; //no qap solution
	}

//clock_t tAll = clock();
	if (printQAP)
		cout << "QAP map computation via " << qapSolverStr << "..\n";
	//compute Q from 2D mesh1.Q2D and mesh2.Q2D (k == 1) or 3D mesh1.Q and mesh2.Q (k > 1)
	int n = mesh1->samples.size(), m = mesh2->samples.size(), rows = n * n, cols = n * n, nIters = 0;
	if (! mesh1->samplesC2F.empty()) //coarse-to-fine mode
		n = mesh1->samplesC2F.size(), m = mesh2->samplesC2F.size(), rows = n * n, cols = n * n; //n=m for sure in this mode
	avgUnary = nu = avgPairwise = np = 0.0; //for statistics only
	MatrixXd Q(rows, cols);
	if (m > n)
		cout << "WARNING: partial mesh2 sample size <= complete model mesh1 inequality not satisfied\n\n"; //= 'cos partialMatching=false mode makes mesh2 complete as well (hence mesh2 has the same sample size now)
	if (k == 1) //2D Q2D in use 'cos entries are 1-D (k=1), where off-diagonals store pure geodesic lengths (no pathing/patching) and diagonals store gaussian curvatures
	{
		MatrixXd X(n, n), Y(n, n); //mesh1 and mesh2 content stored in X and Y, respectively
		if (n == m) //no padding
			for (int s1 = 0; s1 < n; s1++) //Q2D is an s x s matrix where s is the number of samples (same for both models for sure)
				for (int s2 = 0; s2 < n; s2++)
					if (dense)
						X(s1, s2) = mesh1->Q2Ddense[s1][s2], Y(s1, s2) = mesh2->Q2Ddense[s1][s2]; //unary descriptor for diagonals, geodesic length for off-diagonals
					else
						X(s1, s2) = mesh1->Q2D[s1][s2], Y(s1, s2) = mesh2->Q2D[s1][s2]; //unary descriptor for diagonals, geodesic length for off-diagonals
		else //padding is needed for the partial model Y
		{
			for (int s1 = 0; s1 < n; s1++) //Q2D is an s x s matrix where s is the number of samples
				for (int s2 = 0; s2 < n; s2++)
					if (dense)
						X(s1, s2) = mesh1->Q2Ddense[s1][s2]; //unary descriptor for diagonals, geodesic length for off-diagonals
					else
						X(s1, s2) = mesh1->Q2D[s1][s2]; //unary descriptor for diagonals, geodesic length for off-diagonals
			Y.resize(m, m); //n - m entries will be padded below so start w/ m size here
			for (int s1 = 0; s1 < m; s1++) //Q2D is an s x s matrix where s is the number of samples (padding needed here)
				for (int s2 = 0; s2 < m; s2++)
					if (dense)
						Y(s1, s2) = mesh2->Q2Ddense[s1][s2]; //unary descriptor for diagonals, geodesic length for off-diagonals
					else
						Y(s1, s2) = mesh2->Q2D[s1][s2]; //unary descriptor for diagonals, geodesic length for off-diagonals
			Y = padMatrix(Y, n - m); //n > m for sure
		}
//cout << X << "xxxxxx\n\n" << Y << "yyyyyyy\n\n";
		Q = computeQ(X, Y, n - m);
	}
	else //3D Q in use w/ k-d entries, where off-diagonals store a set of patch areas and diagonals store a set of areasRing (padded w/ 0s in fillQ() when sets have size less than k)
	{
		vector< MatrixXd > Xs(k), Ys(k); //one X and one Y for each of the k dimensions (these matrices will be added for the same i, j index in computeQ2())
		for (int dim = 0; dim < k; dim++)
			if (n == m) //no padding
			{
				Xs[dim].resize(n, n); //necessary to allocate memory for this particular entry (not done during declaration above)
				Ys[dim].resize(m, m); //necessary to allocate memory for this particular entry (not done during declaration above)
				for (int  s1 = 0; s1 < n; s1++) //Q is an s x s x k matrix where s is the number of samples (same for both models for sure) and k is the dimension of unary and pairwise features
					for (int s2 = 0; s2 < n; s2++)
						if (dense)
							Xs[dim](s1, s2) = mesh1->Qdense[s1][s2][dim], Ys[dim](s1, s2) = mesh2->Qdense[s1][s2][dim]; //areaRing at dimension dim for diagonals, patch area at dimension dim for off-diagonals
						else
							Xs[dim](s1, s2) = mesh1->Q[s1][s2][dim], Ys[dim](s1, s2) = mesh2->Q[s1][s2][dim]; //areaRing at dimension dim for diagonals, patch area at dimension dim for off-diagonals
			}
			else //padding is needed for the partial model Y 
			{
				Xs[dim].resize(n, n); //necessary to allocate memory for this particular entry (not done during declaration above)
				for (int s1 = 0; s1 < n; s1++) //Q is an s x s x k matrix where s is the number of samples and k is the dimension of unary and pairwise features
					for (int s2 = 0; s2 < n; s2++)
						if (dense)
							Xs[dim](s1, s2) = mesh1->Qdense[s1][s2][dim]; //areaRing at dimension dim for diagonals, patch area at dimension dim for off-diagonals
						else
							Xs[dim](s1, s2) = mesh1->Q[s1][s2][dim]; //areaRing at dimension dim for diagonals, patch area at dimension dim for off-diagonals
				Ys[dim].resize(m, m); //n - m entries will be padded below so start w/ m size here
				for (int s1 = 0; s1 < m; s1++) //Q2D is an s x s matrix where s is the number of samples (padding needed here)
					for (int s2 = 0; s2 < m; s2++)
						if (dense)
							Ys[dim](s1, s2) = mesh2->Qdense[s1][s2][dim]; //areaRing at dimension dim for diagonals, patch area at dimension dim for off-diagonals
						else
							Ys[dim](s1, s2) = mesh2->Q[s1][s2][dim]; //areaRing at dimension dim for diagonals, patch area at dimension dim for off-diagonals
				Ys[dim] = padMatrix(Ys[dim], n - m); //n > m for sure
			}
//cout << Xs[0] << "firstxxxxxx\n\n" << Ys[0] << "firstyyyyyyy\n\n";
//cout << Xs[k-1] << "lastxxxxxx\n\n" << Ys[k-1] << "lastyyyyyyy\n\n";
		do
		{
			Q = computeQ2(Xs, Ys, n - m, maxGeoDiff, barrier, mesh1->coarseMap.size() / 2); //last parameter is the same cSize value used in denseMap(), i.e., number of matches in coarseMap[]




//break;
//consider: multistart with new x0 in qucoop




			if (! infCost || nIters++ > 100)
				break;
			else
				maxGeoDiff *= 2, cout << "maxGeoDiff doubled (qap)\n"; //try mapping w/ a relaxed maxGeoDiff in case corresp->mapOK stayed false in this iteration
		} while (true);
	}
	if (printQAP)
		cout << "Q's avgUnary: " << avgUnary / nu << " vs. avgPairwise: " << avgPairwise / np << " before " << m << "x" << m << " map " //" w/ " << nu << " & " << np << " counts //avgUnary and avgPairwise'll be close if removed setAreaWeight2() used w/o *m-1 (*m-1 added on purpose to make unary more powerful)
			 << "(1 unary + " << m-1 << " pairwise terms" << (noBarrier ? "" : " (optional barrier term is added to each pairwise term)") <<  " will be added for each match cost)\n";
//cout << Q << " q matrix\n\n";

	clock_t t = clock();
	MatrixXd P_sol; //final solution
	if (qapSolver == QUCOOP)
	{
		//QuCOOP solver(Q, {{2, 7}});
		QuCOOP solver(Q, printQAP); //no fixed matches provided	
		solver.startTemp = mesh1->maxGeoDist * mesh1->maxGeoDist; //required by SA() which is replaced by better simulatedAnnealing() so this is redundant
		solver.solve(), P_sol = solver.P_sol; //solver.P_sol ready after solve()
	}
	else if (qapSolver == QMATCH)
	{
		QMatch solver(Q, printQAP);
		solver.solve(), P_sol = solver.P_sol; //solver.P_sol ready after solve()
	}
	else if (qapSolver == QGMATCH)
	{
		QGM solver(Q, printQAP);
		solver.solve(), P_sol = solver.P_sol; //solver.P_sol ready after solve()
/*		VectorXd idx = VectorXd::LinSpaced(n, 0, n-1);
		VectorXd perm_action = solver.P_sol * idx; //gives 1-valued column ids for each row
		cout << "\nResults:            " << perm_action.transpose().cast<int>() 
				<< ", objective: " << objective(Map<VectorXd>(solver.P_sol.data(), n*n), Q) << endl;*/
	}
	else //our novel QA-ready QAP-solver
	{
		QuCOOPMatch solver(Q, printQAP); //alternative to QuCOOP solver(Q, printQAP);
		solver.solve(), P_sol = solver.P_sol; //solver.P_sol ready after solve()
	}
//solver.identityMap(); //source.samples[i] goes to target.samples[i] for a quick debugging (disable solve() above)

	mapOK = true;
	//P_sol to matchIdxs
	for (int i = 0; i < P_sol.rows(); i++)
	{
		int bitAt = -1; //bit 1 is at this index which implies mesh1.samples[i] matches with mesh2.samples[bitAt]
		for (int j = 0; bitAt == -1 && j < P_sol.cols(); j++)
			if (P_sol(i, j) == 1)
				bitAt = j;
		if (bitAt == -1)
			cout << "WARNING: illegal P_sol matrix\n\n";
		int v1 = mesh1->samples[i], v2 = mesh2->samples[bitAt];
		if (! mesh1->samplesC2F.empty()) //coarse-to-fine mode
			v1 = mesh1->samplesC2F[i], v2 = mesh2->samplesC2F[bitAt];
		if (bitAt >= m)
			mesh1->verts[v1]->matchIdx = -1; //matched to a virtual sample, i.e., no match
		else
			mesh1->verts[v1]->matchIdx = v2, mesh2->verts[v2]->matchIdx = v1; //dual operation
	}
	nQAPs++;//cout << ((float) clock()-tAll) / CLOCKS_PER_SEC << " secs for 1 qapMap() all, " << ((float) clock()-t) / CLOCKS_PER_SEC << " for solver only\n";
	if (dense)
		return; //no further evaluation as i'll print puregeo distortions later in dense mode
	//use the latest/optimal matchIdx values to assign individual v.distortion scores for printing only
	if (dissimType != PAIRWISENOPATCH) //pathing/patching in use
		mapQualityVirtual(); //global distortion is updated (mapQualityPatch() uses the same metrics so call either of them here)
	else
		mapQualityGeo(); //global distortion is updated
	char str[] = " -qap- "; double tmp;
	if (printMap(t, str) && maxGeoDiff != INF && dissimType != PAIRWISENOPATCH) //if INF not printed at all then printMap() returns false and i don't reprint the same INF-free map below
		tmp = maxGeoDiff, maxGeoDiff = INF, mapQualityVirtual(), printMap(t, str), maxGeoDiff = tmp; //maxGeoDiff=INF made to see the winner map w/o maxGeoDiff filtering (roll back to the original for denseMap())
}
