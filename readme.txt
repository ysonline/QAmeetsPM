Code and executable for the paper Quantum Annealing Meets Partial Matching, EG, 2027 (submitted).

Executable command:
PartialCorrespQAP.exe <sourceID> <targetID> <nInitialSourceSamples> <mInitialTargetSamples> <nDenseSamples> <qapEnergy> <nPathSamples> <qapSolver> <nSamplesPerSubregion> <nMaxC2FLevels> <fullMap>

Examples:
PartialCorrespQAP.exe 100 101 10 5 -1 1 9 1 10 10 0 //horse complete-to-partial no dense C2F part (-1)
PartialCorrespQAP.exe 100 101 10 5 250 1 9 1 10 10 1 //horse complete-to-partial w/ dense C2F part and fullMap part (last 1)

More Examples:
see sampleRuns.png

This folder also includes drop-in solver functions of our quantum-classical competitors QuCOOP and Q-Match as qucoop.h and qmatch.h, respectively. Our solver is qucoopmatch.h and mesh processing/coarse-to-fine stuff are mainly at Mesh.cpp and Correspondence.cpp.

QAPSolversDWave folder, on the other hand, provides DWave-ready python implementations (that could hav been tested due to their canceled policy on 1-minute-per month free access. Still valuable to be used as a references for those who have access.

You may get visuals by defining #define GRAPHICS_STUFF in Mesh.h although this requires Coin3D installation (dll's provided).

Please cite if you are using this package:
Quantum Annealing Meets Partial Matching, EG, 2027 (submitted).
