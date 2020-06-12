#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/*
	Random number library function
	http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
*/
#include "mt19937ar.h"

#ifdef _MSC_VER
#define 	C_DIR_DELIMITER '\\'
#include 	<direct.h>
#include 	<process.h>
#else
#define 	C_DIR_DELIMITER '/'
#include 	<sys/types.h>
#include 	<sys/stat.h>
#include 	<unistd.h>
#endif

/*
	Stop Visual C++ from warning about thread safety when asked to compile idiomatic ANSI
*/
#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

#define		_MAX_STATIC_BUFF_LEN	1024
#define		MY_REALLOC_BLOCK_SIZE	1024
#define		BRK_ENDLINE				"\r\n"
#define		TOK_WHITESPACE			" \t,"
#define		_MY_TINY_EPS			1e-10
#define		_SCREEN_PRINT_STEP		100

#define		OUT_DIR				"samplingPattern"
#define		HOST_INFO_FILE			"activeLandscape.txt"		/* This file is created by the landscape scale model when it runs */
#define		DEBUG_DUMP_INFO			0							/* Whether (1) or not (0) to dump information at the end to check the calculation */
#define		PARAM_OBJ_FUNC_TYPE		0							/*
																	= 0 means average of detection probabilities over runs
																	= 1 means simulate and count up number of successes
																	= 2 means maximise average of expected detections/run
																*/

/*
	Global variables storing global configuration options
*/
char	INPUT_DIR[_MAX_STATIC_BUFF_LEN];
char	SIM_OUTPUT_STUB[_MAX_STATIC_BUFF_LEN];
char	OBJ_FUNC_OUT[_MAX_STATIC_BUFF_LEN];
int		NUM_RUNS;
int		PARAM_N;
int		PARAM_n;
double	PARAM_WITHIN_CELL_R;
double	PARAM_WITHIN_CELL_S0;
int		PARAM_TRUE_MIN_FLAG;
double	TEST_SENS;
double	DET_LAG;
double	PARAM_DELTA;
double	PARAM_COOL;
double	PARAM_ALPHA;
int		SIMANN_N;
int		B_ALLOW_DUPLICATES;

typedef struct
{
	int				hostID;
	int				hostX;
	int				hostY;
	double			hostDensity;
} t_HostInfo;

typedef struct
{
	int				hostID;
	int				hostPos;
} t_HostLookup;

typedef struct
{
	double			maxTimeInf;
	int				numInf;
	t_HostLookup	*aHostLookup;
	double			*aTimeInf;
	double			*aHostDensity;
	double			*aPDetect;
} t_RunInfo;

typedef struct
{
	int				hostID;
	int				numSims;
} t_InfInfo;

typedef struct
{
	int				numRuns;
	t_RunInfo		*aRunInfo;	/* this stores times of infection and p(detect) for infected hosts in the individual runs */
	int				numHosts;
	t_HostInfo		*aHostInfo;	/* this stores the information on location of hosts */
	int				numInf;
	t_InfInfo		*aInfInfo;	/* this stores information on how frequently hosts infected (used to avoid ever choosing non-infected hosts) */
} t_SSAInfo;

/*
	Utility functions for reading configuration options
*/
int	getCfgFileName(char *szProgName, char *szCfgFile)
{
	/*
		Work out configuration file name from that of the executable
		and check whether it exists by attempting to read it
	*/

	char	*pPtr;
	FILE 	*fp;

	szCfgFile[0] = '\0';
	{
		if ((pPtr = strrchr(szProgName, C_DIR_DELIMITER)) != NULL)
		{
			strcpy(szCfgFile, pPtr + 1);
		}
		else
		{
			strcpy(szCfgFile, szProgName);
		}
		if ((pPtr = strstr(szCfgFile, ".exe")) != NULL)
		{
			*pPtr = '\0';
		}
		strcat(szCfgFile, ".cfg");
	}
	/* check file exists */
	fp = fopen(szCfgFile, "rb");
	if (fp)
	{
		fclose(fp);
		return 1;
	}
	return 0;
}

/*
	Following set of routines find values of parameters from the command line options,
	or, failing that, from the cfg file
*/
int findKey(int argc, char **argv, char*szCfgFile, char *szKey, char *szValue)
{
	char *pVal;
	int	 bRet, i;
	FILE *fp;
	char *pThisPair;
	char *szArgvCopy;

	bRet = 0;
	i = 0;
	/* try to find the relevant key on the command line */
	while (bRet == 0 && i<argc)
	{
		szArgvCopy = strdup(argv[i]);
		if (szArgvCopy)
		{
			pThisPair = strtok(szArgvCopy, " \t");
			while (pThisPair)
			{
				if (strncmp(pThisPair, szKey, strlen(szKey)) == 0)
				{
					pVal = strchr(pThisPair, '=');
					if (pVal)
					{
						/* make sure isn't just start of string matching the key */
						if (pThisPair[strlen(szKey)] == '=')
						{
							strcpy(szValue, pVal + 1);
							fprintf(stdout, "extracted %s->%s from command line\n", szKey, szValue);
							bRet = 1;
						}
					}
				}
				pThisPair = strtok(NULL, " \t");
			}
			free(szArgvCopy);
		}
		i++;
	}
	/* otherwise, look in the cfg file */
	if (bRet == 0)
	{
		fp = fopen(szCfgFile, "rb");
		if (fp)
		{
			char szLine[_MAX_STATIC_BUFF_LEN];

			while (!bRet && fgets(szLine, _MAX_STATIC_BUFF_LEN, fp))
			{
				char *pPtr;
				if ((pPtr = strchr(szLine, '=')) != NULL)
				{
					*pPtr = '\0';
					if (strcmp(szKey, szLine) == 0)
					{
						strcpy(szValue, pPtr + 1);
						/* strip off newline (if any) */
						if ((pPtr = strpbrk(szValue, "\r\n")) != NULL)
							*pPtr = '\0';
						bRet = 1;
					}
				}
			}
			fclose(fp);
		}
	}
	return bRet;
}

int readStringFromCfg(int argc, char **argv, char *szCfgFile, char *szKey, char *szValue)
{
	return(findKey(argc, argv, szCfgFile, szKey, szValue));
}

int readDoubleFromCfg(int argc, char **argv, char *szCfgFile, char *szKey, double *pdValue)
{
	char szValue[_MAX_STATIC_BUFF_LEN];

	if (findKey(argc, argv, szCfgFile, szKey, szValue))
	{
		*pdValue = atof(szValue);
		return 1;
	}
	return 0;
}

int readIntFromCfg(int argc, char **argv, char *szCfgFile, char *szKey, int *pnValue)
{
	char szValue[_MAX_STATIC_BUFF_LEN];

	if (findKey(argc, argv, szCfgFile, szKey, szValue))
	{
		*pnValue = atoi(szValue);
		return 1;
	}
	return 0;
}

/*
	Return uniform random number between 0 and 1

	Encapsulated to allow easy replacement if necessary
*/
double	uniformRandom()
{
#if 0
	int nRet;
	nRet = RAND_MAX;
	while(nRet == RAND_MAX || nRet == 0)
	{
		nRet = rand();
	}
	return ((double)nRet/(double)(RAND_MAX));
#else
	return genrand_real3();
#endif
}

/*
	seed random number generator

	encapsulated to allow easy replacement if necessary
*/
void	seedRandom()
{
	unsigned long		ulnSeed;
	unsigned long		myPID;

	ulnSeed =(unsigned long) time(NULL);
#ifndef _WIN32	/* make sure different processes started at same time have different seeds */
	myPID = (unsigned long) getpid();
#else
	myPID = (unsigned long) _getpid();
#endif
	ulnSeed += myPID;
#if 0
	srand((unsigned int)ulnSeed);
#else
	init_genrand(ulnSeed);
#endif
}

int readHostInfo(t_SSAInfo *pSSAInfo)
{
	int		bRet;
	FILE	*fIn;
	char	szBuff[_MAX_STATIC_BUFF_LEN];
	char	szHostInf[_MAX_STATIC_BUFF_LEN];
	char	*pPtr;
	int		thisTok;
	int		numAlloc;

	bRet = 0;
	sprintf(szHostInf, "%s%s", INPUT_DIR, HOST_INFO_FILE);
	fIn = fopen(szHostInf, "rb");
	if(fIn)
	{
		bRet = 1;
		numAlloc = 0;
		while(bRet && fgets(szBuff, _MAX_STATIC_BUFF_LEN, fIn))
		{
			if(pSSAInfo->numHosts == numAlloc)
			{
				numAlloc += MY_REALLOC_BLOCK_SIZE;
				pSSAInfo->aHostInfo = realloc(pSSAInfo->aHostInfo, sizeof(*pSSAInfo->aHostInfo) * numAlloc);
				if(!pSSAInfo->aHostInfo)
				{
					bRet = 0;
				}
			}
			if(bRet)
			{
				pSSAInfo->aHostInfo[pSSAInfo->numHosts].hostID = pSSAInfo->numHosts;
				thisTok = 0;
				pPtr = strpbrk(szBuff, BRK_ENDLINE);
				if(*pPtr)
				{
					*pPtr = '\0';
				}
				pPtr = strtok(szBuff, TOK_WHITESPACE);
				while(bRet && pPtr)
				{
					switch(thisTok)
					{
					case 0:
						/* x */
						pSSAInfo->aHostInfo[pSSAInfo->numHosts].hostX = atoi(pPtr);
						break;
					case 1:
						/* y */
						pSSAInfo->aHostInfo[pSSAInfo->numHosts].hostY = atoi(pPtr);
						break;
					case 2:
						/* density */
						pSSAInfo->aHostInfo[pSSAInfo->numHosts].hostDensity = atof(pPtr);
						break;
					case 3:
						if(pSSAInfo->numHosts != atoi(pPtr))
						{
							fprintf(stderr, "host mismatch when parsing landscape\n");
							bRet = 0;
						}
						break;
					default:
						break;
					}
					pPtr = strtok(NULL, TOK_WHITESPACE);
					thisTok++;
				}
			}
			pSSAInfo->numHosts++;
		}
		fclose(fIn);
		if(bRet && pSSAInfo->numHosts)
		{
			fprintf(stdout, "readHostInfo():\n\tread %d hosts\n", pSSAInfo->numHosts);
		}
		else
		{
			bRet = 0;
			fprintf(stderr, "couldn't find hosts\n");
		}
	}
	else
	{
		fprintf(stderr, "couldn't read host info file (%s)\n", szHostInf);
	}
	return bRet;
}

static int cmpHostLookup(const void *p1, const void *p2)
{
	t_HostLookup *pHL1 = (t_HostLookup *)p1;
	t_HostLookup *pHL2 = (t_HostLookup *)p2;

	return pHL1->hostID - pHL2->hostID;
}

static int cmpInfInfo(const void *p1, const void *p2)
{
	t_InfInfo *pHL1 = (t_InfInfo *)p1;
	t_InfInfo *pHL2 = (t_InfInfo *)p2;

	return pHL1->hostID - pHL2->hostID;
}

int	readSims(t_SSAInfo *pSSAInfo)
{
	int		bRet;
	int		i;
	char	szInputFile[_MAX_STATIC_BUFF_LEN];
	char	szEndTimeFile[_MAX_STATIC_BUFF_LEN];
	FILE	*fIn,*fEnd;
	char	szBuffer[_MAX_STATIC_BUFF_LEN];
	char	*pPtr;
	int		thisTok;
	int		numAlloc;
	int		hostID;
	double	readMaxTime;
	int		everInfAlloc;
	int		infFrom,infTo;

	fprintf(stdout, "readSims()\n");
	everInfAlloc = 0;
	bRet = 1;
	pSSAInfo->numRuns = NUM_RUNS;

	/* Infer number of runs from the input files themselves */
	if (pSSAInfo->numRuns < 0)
	{
		int foundGoodFile,lastRunNumber;

		fprintf(stdout, "Inferring number of runs\n");
		pSSAInfo->numRuns = 0;
		do
		{
			foundGoodFile = 0;
			sprintf(szInputFile, "%s%s_%d.txt", INPUT_DIR, SIM_OUTPUT_STUB, pSSAInfo->numRuns);
			if (fIn = fopen(szInputFile, "rb"))
			{
				fprintf(stdout, "\tFound input %s\n", szInputFile);
				foundGoodFile = 1;
				fclose(fIn);
				pSSAInfo->numRuns++;
			}
		} while (foundGoodFile);
		fprintf(stdout, "\tInferred %d input files\n", pSSAInfo->numRuns);
		/* Check that this is correct based on lastRunNumber file */
		sprintf(szInputFile, "%slastRunNumber.txt", INPUT_DIR);
		fIn = fopen(szInputFile, "rb");
		lastRunNumber = -1;
		if (fIn)
		{
			fgets(szBuffer, _MAX_STATIC_BUFF_LEN, fIn);
			lastRunNumber = atoi(szBuffer);
			fclose(fIn);
		}
		if (pSSAInfo->numRuns == lastRunNumber)
		{
			fprintf(stderr, "\t\tCorrectly...");
		}
		else
		{
			fprintf(stderr, "\t\tIncorrect number inferred...exiting");
			pSSAInfo->numRuns = 0;
		}
	}
	if (pSSAInfo->numRuns)
	{
		pSSAInfo->aRunInfo = malloc(sizeof(*pSSAInfo->aRunInfo)*pSSAInfo->numRuns);
		if (pSSAInfo->aRunInfo)
		{
			for (i = 0; bRet && i < pSSAInfo->numRuns; i++)
			{
				sprintf(szInputFile, "%s%s_%d.txt", INPUT_DIR, SIM_OUTPUT_STUB, i);
				fprintf(stdout, "\t%s_%d\n", SIM_OUTPUT_STUB, i);
				fIn = fopen(szInputFile, "rb");
				if (fIn)
				{
					pSSAInfo->aRunInfo[i].numInf = 0;
					pSSAInfo->aRunInfo[i].aHostLookup = NULL;
					pSSAInfo->aRunInfo[i].aPDetect = pSSAInfo->aRunInfo[i].aTimeInf = pSSAInfo->aRunInfo[i].aHostDensity = NULL;
					numAlloc = 0;
					while (bRet && fgets(szBuffer, _MAX_STATIC_BUFF_LEN, fIn))
					{
						if (pSSAInfo->aRunInfo[i].numInf == numAlloc)
						{
							numAlloc += MY_REALLOC_BLOCK_SIZE;
							pSSAInfo->aRunInfo[i].aHostLookup = realloc(pSSAInfo->aRunInfo[i].aHostLookup, sizeof(*pSSAInfo->aRunInfo[i].aHostLookup) * numAlloc);
							pSSAInfo->aRunInfo[i].aPDetect = realloc(pSSAInfo->aRunInfo[i].aPDetect, sizeof(*pSSAInfo->aRunInfo[i].aPDetect) * numAlloc);
							pSSAInfo->aRunInfo[i].aTimeInf = realloc(pSSAInfo->aRunInfo[i].aTimeInf, sizeof(*pSSAInfo->aRunInfo[i].aTimeInf) * numAlloc);
							pSSAInfo->aRunInfo[i].aHostDensity = realloc(pSSAInfo->aRunInfo[i].aHostDensity, sizeof(*pSSAInfo->aRunInfo[i].aHostDensity) * numAlloc);
							if (!(pSSAInfo->aHostInfo && pSSAInfo->aRunInfo[i].aPDetect && pSSAInfo->aRunInfo[i].aTimeInf && pSSAInfo->aRunInfo[i].aHostDensity))
							{
								bRet = 0;
							}
						}
						if (bRet)
						{
							thisTok = 0;
							pPtr = strpbrk(szBuffer, BRK_ENDLINE);
							if (*pPtr)
							{
								*pPtr = '\0';
							}
							pPtr = strtok(szBuffer, TOK_WHITESPACE);
							while (pPtr)
							{
								switch (thisTok)
								{
								case 12:
									hostID = atoi(pPtr);
									pSSAInfo->aRunInfo[i].aHostLookup[pSSAInfo->aRunInfo[i].numInf].hostID = hostID;
									pSSAInfo->aRunInfo[i].aHostLookup[pSSAInfo->aRunInfo[i].numInf].hostPos = pSSAInfo->aRunInfo[i].numInf;
									if (pSSAInfo->numInf == everInfAlloc)
									{
										everInfAlloc += MY_REALLOC_BLOCK_SIZE;
										pSSAInfo->aInfInfo = realloc(pSSAInfo->aInfInfo, sizeof(*pSSAInfo->aInfInfo)*everInfAlloc);
										if (!pSSAInfo->aInfInfo)
										{
											bRet = 0;
										}
									}
									if (bRet)
									{
										pSSAInfo->aInfInfo[pSSAInfo->numInf].hostID = hostID;
										pSSAInfo->numInf++;
									}
									break;
								case 2:
									/* time */
									pSSAInfo->aRunInfo[i].aTimeInf[pSSAInfo->aRunInfo[i].numInf] = atof(pPtr);
									break;
								case 6:
									/* host density */
									pSSAInfo->aRunInfo[i].aHostDensity[pSSAInfo->aRunInfo[i].numInf] = atof(pPtr);
									break;
								default:
									break;
								}
								thisTok++;
								pPtr = strtok(NULL, TOK_WHITESPACE);
							}
						}
						pSSAInfo->aRunInfo[i].numInf++;
					}
					fclose(fIn);
					sprintf(szEndTimeFile, "%sendTime_%d.txt", INPUT_DIR, i);
					readMaxTime = 0.0;
					fEnd = fopen(szEndTimeFile, "rb");
					if (fEnd)
					{
						fgets(szBuffer, _MAX_STATIC_BUFF_LEN, fEnd);
						readMaxTime = atof(szBuffer);
						pSSAInfo->aRunInfo[i].maxTimeInf = readMaxTime;
						fprintf(stdout, "\t\tread %d infections (maxTime=%.4f)\n", pSSAInfo->aRunInfo[i].numInf, pSSAInfo->aRunInfo[i].maxTimeInf);
						qsort(pSSAInfo->aRunInfo[i].aHostLookup, pSSAInfo->aRunInfo[i].numInf, sizeof(*pSSAInfo->aRunInfo[i].aHostLookup), cmpHostLookup);
						fclose(fEnd);
					}
					else
					{
						bRet = 0;
						fprintf(stderr, "couldn't read %s\n", szEndTimeFile);
					}
				}
				else
				{
					bRet = 0;
					fprintf(stderr, "couldn't open %s\n", szInputFile);
				}
			}
			if (bRet)
			{
				fprintf(stdout, "\t%d hosts infected in total\n", pSSAInfo->numInf);
				qsort(pSSAInfo->aInfInfo, pSSAInfo->numInf, sizeof(*pSSAInfo->aInfInfo), cmpInfInfo);
				pSSAInfo->aInfInfo[0].numSims = 1;
				infTo = 0;
				for (infFrom = 1; infFrom < pSSAInfo->numInf; infFrom++)
				{
					if (pSSAInfo->aInfInfo[infTo].hostID == pSSAInfo->aInfInfo[infFrom].hostID)
					{
						pSSAInfo->aInfInfo[infTo].numSims++;
					}
					else
					{
						infTo++;
						pSSAInfo->aInfInfo[infTo].hostID = pSSAInfo->aInfInfo[infFrom].hostID;
						pSSAInfo->aInfInfo[infTo].numSims = 1;
					}
				}
				pSSAInfo->numInf = infTo + 1;
				fprintf(stdout, "\t%d unique hosts infected\n", pSSAInfo->numInf);
			}
		}
		else
		{
			bRet = 0;
		}
	}
	else
	{
		bRet = 0;
	}
	return bRet;
}

/*
 Handle build up of detectability
	sigma(t) = 1./(1+J.*exp(-r.*t));

 where
		- t = time since first infection of the cell
		- sigma(t) = detectability of cell at time t
		- r = (logistic growth) rate at which infectivity increases
		- J = (1-w0)/w0;
		- detLag = period within which detection is not possible following infection
*/
double detectProbSingleSurvey(double timeSinceInf, int numSamples, double dHostDensity)
{
	/*
		use 1-p(do not detect)
	*/
	double pDetectSingleSample;
	double J;
	double thisWCM;

	if (PARAM_TRUE_MIN_FLAG)
	{
		if (PARAM_WITHIN_CELL_S0 >= dHostDensity)
		{
			J = 0.0;
		}
		else
		{
			thisWCM = PARAM_WITHIN_CELL_S0 / dHostDensity;
			J = (1.0 - thisWCM) / thisWCM;
		}
	}
	else
	{
		J = (1.0 - PARAM_WITHIN_CELL_S0) / PARAM_WITHIN_CELL_S0;
	}
	if (timeSinceInf < DET_LAG)
	{
		pDetectSingleSample = 0.0;
	}
	else
	{
		pDetectSingleSample = TEST_SENS * (1.0 / (1.0 + J * exp(-PARAM_WITHIN_CELL_R * (timeSinceInf - DET_LAG))));
	}
	/*
		Probability of one or more detections (given numSamples taken)
	*/
	return 1.0 - pow(1.0 - pDetectSingleSample, numSamples);
}

int calcProbDetect(t_SSAInfo *pSSAInfo)
{
	int				numToAverage,bRet,i,j,k,numSurveys;
	double			timeSurvey,timeInf,pDetect,pDontDetect,pDetectThisTime,firstOffset,hostDensity;

	fprintf(stdout, "calcProbDetect()\n");
	bRet = 1;
	for(i=0;i<pSSAInfo->numRuns;i++)
	{
		fprintf(stdout, "\tdoing simulation %d\n", i);
		numSurveys = (int)((pSSAInfo->aRunInfo[i].maxTimeInf+_MY_TINY_EPS)/PARAM_DELTA) + 1;
		for(j=0;j<pSSAInfo->aRunInfo[i].numInf;j++)
		{
			timeInf = pSSAInfo->aRunInfo[i].aTimeInf[j];
			hostDensity = pSSAInfo->aRunInfo[i].aHostDensity[j];
			pDetect = 0.0;
			numToAverage = 0;
			firstOffset = 0.0;
			/*
				This loops over the different times surveying could start relative to time of first infection (step of one day)
			*/
			while(firstOffset < PARAM_DELTA)
			{
#ifdef _EXHAUSTIVE_SURVEY_DIAGNOSTICS
				fprintf(stdout, "****\n");
				fprintf(stdout, "surveyOffset=%.4f for host number %d to be infected (hostID=%d,firstInf=%.4f) in simulation %d (maxTimeInf=%.4f)\n",firstOffset,j, pSSAInfo->aRunInfo[i].aHostLookup[j].hostID, timeInf, i, pSSAInfo->aRunInfo[i].maxTimeInf);
				fprintf(stdout, "****\n");
#endif
				pDontDetect = 1.0;
				for(k=0;k<numSurveys;k++)
				{
					timeSurvey = firstOffset + k * PARAM_DELTA;
					/*
						Avoid doing any sample that is actually too late
					*/
					if(timeSurvey < timeInf)
					{
						pDetectThisTime = 0.0;
#ifdef _EXHAUSTIVE_SURVEY_DIAGNOSTICS
						fprintf(stdout, "do survey at %f (ignored in calculation since before first infection of this cell at %.4f)\n", timeSurvey, timeInf);
#endif
					}
					else
					{
						if(timeSurvey > pSSAInfo->aRunInfo[i].maxTimeInf)
						{
							pDetectThisTime = 0.0;
#ifdef _EXHAUSTIVE_SURVEY_DIAGNOSTICS
							fprintf(stdout, "do survey at %f (ignored in calculation since after maximum incidence reached at %.4f)\n", timeSurvey, pSSAInfo->aRunInfo[i].maxTimeInf);
#endif
						}
						else
						{
							pDetectThisTime = detectProbSingleSurvey(timeSurvey - timeInf, PARAM_n, hostDensity);
#ifdef _EXHAUSTIVE_SURVEY_DIAGNOSTICS
							fprintf(stdout, "do survey at %f\n", timeSurvey);
#endif
						}
					}
					pDontDetect *= (1.0-pDetectThisTime);
				}
				pDetect += (1.0-pDontDetect);
				firstOffset+=(1.0/365.0);
				numToAverage++;
			}
			pDetect /= numToAverage;
			pSSAInfo->aRunInfo[i].aPDetect[j] = pDetect;
		}
	}
	return bRet;
}

double calcObjFunction(t_SSAInfo *pSSAInfo, int numToSurvey, int *anHostID)
{
	int				thisHost,i,j;
	double			randomDraw,objFunc,pNotDetectOverall,pDetectOnThisHost,pDetectFromThisPatternInThisRun;
	double			expectedFindsThisRun; /* Calculating the expected maximum number of finds before the disease reaches a certain incidence */
	t_HostLookup	sLookup,*pFound;

	objFunc = 0.0;
	for(i=0; i < pSSAInfo->numRuns;i++)
	{
		expectedFindsThisRun = 0.0;
		pNotDetectOverall = 1.0;
		for(j=0;j<numToSurvey;j++)
		{
			thisHost = anHostID[j];
			sLookup.hostID = thisHost;
			pFound = bsearch(&sLookup,pSSAInfo->aRunInfo[i].aHostLookup, pSSAInfo->aRunInfo[i].numInf, sizeof(*pSSAInfo->aRunInfo[i].aHostLookup), cmpHostLookup);
			if(pFound)
			{
				pDetectOnThisHost = pSSAInfo->aRunInfo[i].aPDetect[pFound->hostPos];
			}
			else
			{
				pDetectOnThisHost = 0.0;
			}
			pNotDetectOverall *= (1.0 - pDetectOnThisHost); /*PRODUCT OF NO DETECTION PROBABILITIES*/
			expectedFindsThisRun += pDetectOnThisHost; /* ADDS UP EXPECTED NUMBER OF FINDS */
		}

		pDetectFromThisPatternInThisRun = 1 - pNotDetectOverall;
		switch(PARAM_OBJ_FUNC_TYPE)
		{
		case 0:
			objFunc += pDetectFromThisPatternInThisRun; /*THIS CALCULATES THE SUM OF THE PROBS OF DETECTION WHEN ALL J SAMPLED, SUMMED OVER ALL RUNS*/
			break;
		case 1:
			randomDraw = uniformRandom();
			if(randomDraw < pDetectFromThisPatternInThisRun)
			{
				objFunc++; /*THIS CALCULATES SUM OF DETECTIONS OR NOT WHEN ALL J SAMPLED, SUMMED OVER ALL RUNS*/
			}
			break;
		case 2:
			objFunc += expectedFindsThisRun;
			break;
		default:
			fprintf(stderr, "Not implemented\n");
		}
	} /*DOES THIS FOR ALL I RUNS*/
	objFunc /= (double)pSSAInfo->numRuns; /*THIS TAKES THE MEAN OF WHATEVER FORM THE OF TAKES AS SPECIFIED ABOVE*/
	return objFunc;
}

/*
	Will only ever find hosts not already being sampled
*/
int	randomValidHost(t_SSAInfo *pSSAInfo, int bAllowDups, int numToCheck, int *anPatternToCheck, int newPos)
{
	int hostIndex,bHostOK,i;

	bHostOK = 0;
	do
	{
		hostIndex = (int)(uniformRandom()*(double)pSSAInfo->numInf);
		bHostOK = 1;
		/* if host is already the pattern, flag it */
		for (i = 0; i < numToCheck; i++)
		{
			if (pSSAInfo->aInfInfo[hostIndex].hostID == anPatternToCheck[i])
			{
/*
				fprintf(stdout, "\tChoice of %d for element %d of new pattern clashes with element %d of current pattern [%d]...rechoosing\n", pSSAInfo->aInfInfo[hostIndex].hostID, newPos, i, anPatternToCheck[i]);
*/
				bHostOK = 0;
			}
		}
	} while (bHostOK == 0 && bAllowDups == 0);
	return pSSAInfo->aInfInfo[hostIndex].hostID;
}

void debugDumpInfo(t_SSAInfo *pSSAInfo)
{
	int		i,j;
	FILE	*fOut;
	char szOutFile[_MAX_STATIC_BUFF_LEN];

	fprintf(stdout, "debugDumpInfo():\n");
	sprintf(szOutFile, "%s//debug_HostInfo.txt", OUT_DIR);
	fOut = fopen(szOutFile, "wb");
	if(!fOut)
	{
		fprintf(stderr, "couldn't open %s", szOutFile);
		return;
	}
	for(i=0;i<pSSAInfo->numHosts;i++)
	{
		fprintf(fOut, "%d %d %d %f\n", pSSAInfo->aHostInfo[i].hostID, pSSAInfo->aHostInfo[i].hostX, pSSAInfo->aHostInfo[i].hostY, pSSAInfo->aHostInfo[i].hostDensity);
	}
	fclose(fOut);

	sprintf(szOutFile, "%s//debug_UniqInf.txt", OUT_DIR);
	fOut = fopen(szOutFile, "wb");
	if(!fOut)
	{
		fprintf(stderr, "couldn't open %s", szOutFile);
		return;
	}
	for(i=0;i<pSSAInfo->numInf;i++)
	{
		fprintf(fOut, "%d %d\n", pSSAInfo->aInfInfo[i].hostID, pSSAInfo->aInfInfo[i].numSims);
	}
	fclose(fOut);

	/* Debug files: jused to check that simulation data is being read in processed correctly */
	/* One of these is created for each simulation, and the number of lines is equal to the total number of infections */
	for(j=0;j<pSSAInfo->numRuns;j++)
	{
		sprintf(szOutFile, "%s//debug_PDetect_%d.txt", OUT_DIR, j);
		fOut = fopen(szOutFile, "wb");
		if(!fOut)
		{
			fprintf(stderr, "couldn't open %s", szOutFile);
			return;
		}
		for(i=0;i<pSSAInfo->aRunInfo[j].numInf;i++)
		{
			fprintf(fOut, "%d %d %d %d %f %f %f\n", /* Output is "in a sort of random order" */
				pSSAInfo->aRunInfo[j].aHostLookup[i].hostPos,  /* Order of infection in simulation */
				pSSAInfo->aRunInfo[j].aHostLookup[i].hostID, /* Host ID */
				pSSAInfo->aHostInfo[pSSAInfo->aRunInfo[j].aHostLookup[i].hostID].hostX, /* Raster column ID */
				pSSAInfo->aHostInfo[pSSAInfo->aRunInfo[j].aHostLookup[i].hostID].hostY, /* Raster row ID */
				pSSAInfo->aRunInfo[j].aTimeInf[pSSAInfo->aRunInfo[j].aHostLookup[i].hostPos], /* Time of first infection */
				pSSAInfo->aRunInfo[j].aHostDensity[pSSAInfo->aRunInfo[j].aHostLookup[i].hostPos], /* Host density */
				pSSAInfo->aRunInfo[j].aPDetect[pSSAInfo->aRunInfo[j].aHostLookup[i].hostPos]); /* Probability of detection */
		}
		fclose(fOut);
	}
}

int readParams(int argc, char **argv)
{
	int		bRet = 1;
	char	szCfgFile[_MAX_STATIC_BUFF_LEN];

	fprintf(stdout, "readParams()\n");
	if (!getCfgFileName(argv[0], szCfgFile))
	{
		fprintf(stderr, "Couldn't find cfg file for program name '%s'\n", argv[0]);
		return 0;
	}
	if (!readStringFromCfg(argc, argv, szCfgFile, "inputDirectory", INPUT_DIR))
	{
		fprintf(stdout, "Couldn't read inputDirectory\n");
		return 0;
	}
	if (!readStringFromCfg(argc, argv, szCfgFile, "outStub", SIM_OUTPUT_STUB))
	{
		fprintf(stdout, "Couldn't read outStub\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "numIts", &NUM_RUNS))
	{
		fprintf(stdout, "Couldn't read numIts\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "allowDuplicates", &B_ALLOW_DUPLICATES))
	{
		fprintf(stdout, "Couldn't read allowDuplicates\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "numSites", &PARAM_N))
	{
		fprintf(stdout, "Couldn't read numSites\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "samplesPerSite", &PARAM_n))
	{
		fprintf(stdout, "Couldn't read samplesPerSite\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "withinCellBulkUp", &PARAM_WITHIN_CELL_R))
	{
		fprintf(stdout, "Couldn't read withinCellBulkUp\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "withinCellMin", &PARAM_WITHIN_CELL_S0))
	{
		fprintf(stdout, "Couldn't read withinCellMin\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "testSens", &TEST_SENS))
	{
		fprintf(stdout, "Couldn't read testSens\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "detLag", &DET_LAG))
	{
		fprintf(stdout, "Couldn't read detLag\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "delta", &PARAM_DELTA))
	{
		fprintf(stdout, "Couldn't read delta\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "cool", &PARAM_COOL))
	{
		fprintf(stdout, "Couldn't read cool\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "alpha", &PARAM_ALPHA))
	{
		fprintf(stdout, "Couldn't read alpha\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "simann_n", &SIMANN_N))
	{
		fprintf(stdout, "Couldn't read simann_n\n");
		return 0;
	}
	if (!readStringFromCfg(argc, argv, szCfgFile, "objFuncOut", OBJ_FUNC_OUT))
	{
		fprintf(stdout, "Couldn't read objFuncOut\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "trueMinFlag", &PARAM_TRUE_MIN_FLAG))
	{
		fprintf(stdout, "Couldn't read trueMinFlag\n");
		return 0;
	}
	if (bRet)
	{
		/* dump out parameters as read and understood by the programme */
		FILE* paramsOut;
		char	outFile[_MAX_STATIC_BUFF_LEN];

		sprintf(outFile, "%s%c%s_paramsOut.txt", OUT_DIR, C_DIR_DELIMITER, SIM_OUTPUT_STUB);
		fprintf(stdout, "\twriting copy of parameters to %s\n", outFile);
		paramsOut = fopen(outFile, "wb");
		if (paramsOut)
		{
			fprintf(paramsOut, "inputDirectory=%s\n", INPUT_DIR);
			fprintf(paramsOut, "outStub=%s\n", SIM_OUTPUT_STUB);
			fprintf(paramsOut, "numIts=%d\n", NUM_RUNS);
			fprintf(paramsOut, "allowDuplicates=%d\n", B_ALLOW_DUPLICATES);
			fprintf(paramsOut, "numSites=%d\n", PARAM_N);
			fprintf(paramsOut, "samplesPerSite=%d\n", PARAM_n);
			fprintf(paramsOut, "withinCellBulkUp=%f\n", PARAM_WITHIN_CELL_R);
			fprintf(paramsOut, "withinCellMin=%f\n", PARAM_WITHIN_CELL_S0);
			fprintf(paramsOut, "testSens=%f\n", TEST_SENS);
			fprintf(paramsOut, "detLag=%f\n", DET_LAG);
			fprintf(paramsOut, "delta=%f\n", PARAM_DELTA);
			fprintf(paramsOut, "cool=%f\n", PARAM_COOL);
			fprintf(paramsOut, "alpha=%f\n", PARAM_ALPHA);
			fprintf(paramsOut, "simann_n=%d\n", SIMANN_N);
			fprintf(paramsOut, "objFuncOut=%s\n", OBJ_FUNC_OUT);
			fprintf(paramsOut, "trueMinFlag=%d\n", PARAM_TRUE_MIN_FLAG);
			fclose(paramsOut);
		}
		else
		{
			fprintf(stderr, "couldn't dump parameters...exiting\n");
			return 0;
		}
	}
	return bRet;
}

int main(int argc, char **argv)
{
	t_SSAInfo	sSSAInfo;
	int			bContinue;

	seedRandom();
	memset(&sSSAInfo,0,sizeof(sSSAInfo));
#ifdef _MSC_VER
	if(!mkdir(OUT_DIR))
#else
	if(!mkdir(OUT_DIR, 0777))
#endif
	{
		fprintf(stdout, "created directory %s for output\n", OUT_DIR);
	}
	else
	{
		fprintf(stdout, "directory %s already exists\n", OUT_DIR);
	}
	bContinue = 1;
	/* read in parameters (stored globally) */
	if (!readParams(argc,argv))
	{
		fprintf(stderr, "couldn't read parameters\n");
		bContinue = 0;
	}
	/* read in host info */
	if(bContinue && !readHostInfo(&sSSAInfo))
	{
		fprintf(stderr, "couldn't read host info\n");
		bContinue = 0;
	}
	/* read in inf info */
	if(bContinue && !readSims(&sSSAInfo))
	{
		fprintf(stderr, "couldn't read simulation info\n");
		bContinue = 0;
	}
	/* calculate detection probabilities */
	if(bContinue && !calcProbDetect(&sSSAInfo))
	{
		fprintf(stderr, "couldn't calculate detection probabilities\n");
		bContinue = 0;
	}
	if(bContinue && DEBUG_DUMP_INFO)
	{
		debugDumpInfo(&sSSAInfo);
	}
	/* do annealing */
	if(bContinue)
	{
		int				j,i,thisHost,randToChange,oldVal;
		int				*anPattern, simann_n;
		double			oldObj,newObj,objFunc, probAcc, dice, cool, alpha;
		double			pDetect,pNotDetectThisRun,timeInf;
		t_HostLookup	sLookup,*pFound;
		FILE			*fpObjOut,*fpFakeOutput;

		cool = PARAM_COOL;
		alpha = PARAM_ALPHA;
		simann_n = SIMANN_N;

		anPattern = malloc(sizeof(int) * PARAM_N);
		if(anPattern)
		{
			char szObjOutFile[_MAX_STATIC_BUFF_LEN];

			sprintf(szObjOutFile, "%s//%s", OUT_DIR, OBJ_FUNC_OUT);
			fpObjOut = fopen(szObjOutFile, "wb");
			if (fpObjOut)
			{
				/* create an initial pattern */
				for (i = 0; i < PARAM_N; i++)
				{
					/* Random host from those that were ever infected */
					anPattern[i] = randomValidHost(&sSSAInfo, B_ALLOW_DUPLICATES, i, anPattern, i);
				}
				/* find its objective function */
				oldObj = calcObjFunction(&sSSAInfo, PARAM_N, anPattern);
				/* do SIMANN_N iterations of the spatial annealing */
				for (j = 0; j <= SIMANN_N; j++)
				{
					/* randomly choose a position in the pattern to change */
					randToChange = (int)(uniformRandom()*PARAM_N);
					/* change it (storing old value) */
					oldVal = anPattern[randToChange];
					anPattern[randToChange] = randomValidHost(&sSSAInfo, B_ALLOW_DUPLICATES, PARAM_N, anPattern, randToChange);
					/* calculate new objective function */
					newObj = calcObjFunction(&sSSAInfo, PARAM_N, anPattern);
					/* always accept change if change increases objective function */
					if(newObj > oldObj)
					{
						probAcc = 1.;
					}
					/* if ratio of difference to cooling parameter is very negative (lower than -99), do not swap */
					else if ((newObj - oldObj)/cool < -99)
					{
						probAcc = 0.;
					}
					/* if objective function increases accept swap with this prob */
					else
					{
						probAcc = exp((newObj - oldObj)/cool);
					}
					/* generate a random number between 0 and 1 */
					dice = uniformRandom();

					/* update obj function */
					if(dice<probAcc)
					{
						oldObj = newObj;
					}
					/* revert change */
					else
					{
						anPattern[randToChange] = oldVal;
					}

					/* dump out information on the objective function */
					if((j%_SCREEN_PRINT_STEP)==0)
					{
						/* print to screen every so often */
						fprintf(stdout, "sample %d (%f)...\n\t", j, oldObj);
						for (i = 0; i < PARAM_N; i++)
						{
							fprintf(stdout, " %d", anPattern[i]);
						}
						fprintf(stdout, "\n");
					}
						fprintf(fpObjOut, "%d %.4f", j, oldObj);
						for (i = 0; i < PARAM_N; i++)
						{
							fprintf(fpObjOut, " %d", anPattern[i]);
						}
						fprintf(fpObjOut, "\n");

					/*
						Reduce the cooling parameter
					*/
					cool = cool*alpha;
				}
				/*
					Only used in debugging
				*/
				if (DEBUG_DUMP_INFO)
				{
					char szOutFile[_MAX_STATIC_BUFF_LEN];

					sprintf(szOutFile, "%s//debugOutput.txt", OUT_DIR);
					/* print out information on the actual info on last pattern */
					fpFakeOutput = fopen(szOutFile, "wb");
					if (fpFakeOutput)
					{
						fprintf(fpFakeOutput, "pattern:\n");
						for (i = 0; i < PARAM_N; i++)
						{
							fprintf(fpFakeOutput, "\tentry=%d, id=%d, x=%d, y=%d, density=%.4f\n",
								i,												/* position in pattern */
								anPattern[i],									/* hostID */
								sSSAInfo.aHostInfo[anPattern[i]].hostX,			/* hostX */
								sSSAInfo.aHostInfo[anPattern[i]].hostY,			/* hostY */
								sSSAInfo.aHostInfo[anPattern[i]].hostDensity);	/* hostDensity */
						}
						fprintf(fpFakeOutput, "pattern performance by run:\n");
						/* print out information on times and probabilities for this pattern from each sim... */
						for (j = 0; j < sSSAInfo.numRuns; j++)
						{
							fprintf(fpFakeOutput, "\trun %d (maxTime=%.4f):\n", j, sSSAInfo.aRunInfo[j].maxTimeInf);
							pNotDetectThisRun = 1.0;
							for (i = 0; i < PARAM_N; i++)
							{
								thisHost = sSSAInfo.aHostInfo[anPattern[i]].hostID;
								sLookup.hostID = thisHost;
								pFound = bsearch(&sLookup, sSSAInfo.aRunInfo[j].aHostLookup, sSSAInfo.aRunInfo[j].numInf, sizeof(*sSSAInfo.aRunInfo[j].aHostLookup), cmpHostLookup);
								if (pFound)
								{
									timeInf = sSSAInfo.aRunInfo[j].aTimeInf[pFound->hostPos];
									pDetect = sSSAInfo.aRunInfo[j].aPDetect[pFound->hostPos];
								}
								else
								{
									timeInf = -1.0;
									pDetect = 0.0;
								}
								pNotDetectThisRun *= (1.0 - pDetect);
								fprintf(fpFakeOutput, "\t\tentry=%d, id=%d, x=%d, y=%d, density=%.4f, t=%.4f (p=%.4f)\n", i, anPattern[i], sSSAInfo.aHostInfo[anPattern[i]].hostX, sSSAInfo.aHostInfo[anPattern[i]].hostY, sSSAInfo.aHostInfo[anPattern[i]].hostDensity, timeInf, pDetect);
							}
							fprintf(fpFakeOutput, "\t\t\t=> pDetectThisRun=%.4f\n", 1.0 - pNotDetectThisRun);
						}
						objFunc = calcObjFunction(&sSSAInfo, PARAM_N, anPattern);
						fprintf(fpFakeOutput, "calculated objective function:\n");
						fprintf(fpFakeOutput, "\t%.4f\n", objFunc);
						fclose(fpFakeOutput);
					}
				}
				fclose(fpObjOut);
			}
			free(anPattern);
		}
		else
		{
			fprintf(stderr, "Memory error; exiting");
			exit(EXIT_FAILURE);
		}

		/* free up the memory that has been allocated */
		if(sSSAInfo.numHosts)
		{
			free(sSSAInfo.aHostInfo);
		}
		if(sSSAInfo.numInf)
		{
			free(sSSAInfo.aInfInfo);
		}
		if(sSSAInfo.numRuns)
		{
			for(i=0;i<sSSAInfo.numRuns;i++)
			{
				if(sSSAInfo.aRunInfo[i].numInf)
				{
					free(sSSAInfo.aRunInfo[i].aHostLookup);
					free(sSSAInfo.aRunInfo[i].aPDetect);
					free(sSSAInfo.aRunInfo[i].aTimeInf);
					free(sSSAInfo.aRunInfo[i].aHostDensity);
				}
			}
			free(sSSAInfo.aRunInfo);
		}
	}
	return EXIT_SUCCESS;
}
