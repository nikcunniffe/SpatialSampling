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

/*
	Stop Visual C++ from warning about thread safety when asked to compile idiomatic ANSI
*/
#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

/*
	Constants used throughout program
*/
#define		_EMPTY_CELL						-1
#define		_UNDEF_TIME						-1.0
#define		_MAX_STATIC_BUFF_LEN			1024
#define		_MAX_DYNAMIC_BUFF_LEN			(1024*1024)
#define		_GIS_HEADER_LENGTH				6
#define		_GIS_WHITESPACE					",\t "
#define		_GIS_NEWLINES					"\r\n"
#define		_GIS_NODATA						"-9999"
#define		_LANDSCAPE_BLOCK_SIZE			128
#define		_PI								3.1415926535897932384626433
#define		_PRI_INF_TYPE					1
#define		_SEC_INF_TYPE					2
#define		_SETUP_DISPERSAL_PRINT_DOT		10000
#define		_VERY_LONG_TIME					10000000.0

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
	Information on a single cell
*/
typedef struct
{
	int		xPos;			/* column in gis raster */
	int		yPos;			/* row in gis raster */
	double	propFull;		/* proportion of cell with host */
	double	relInf;			/* relative infectivity */
	double	relSus;			/* relative susceptibility */
	double	relPri;			/* relative force of primary infection */
	double	tInf;			/* time of first infection of this cell */
	double	tNext;			/* time of next possible secondary infection caused by this cell */
	int		infType;		/* whether this cell became infected via primary or secondary infection */
	int		infBy;			/* which host infected (=_EMPTY_CELL for primary) */
} t_Cell;

/*
	Cells are collated in a landscape
*/
typedef struct
{
	int		numRows;
	int		numCols;
	t_Cell	*aCells;
	int		numCells;
	int		numCellsSpace;
	int		*aCellLookup;
	double	totalFull;		/* this stores the total number of cells that are full, accounting for fractions */
} t_Landscape;

/*
	Store the rate of primary infection on each cell
	Keep track of the cumulative pressure to make it easy to find which cell is going to be (primary) infected next
*/
typedef struct
{
	double	*aCumPressure;
	double	totalPressure;
	double	ratePri;
	double	nextT;
} t_PriInf;

/*
	Progress of a single epidemic is stored by keeping track of which cells have become infected
	And a priority queue for the next infection each generates
*/
typedef struct
{
	int		*aQueueCells;
	int		queueLen;
	int		*aInfCells;
	int		totalInf;
} t_Epidemic;

/*
	Parameters which control behaviour of simulation
*/
typedef struct
{
	double	cellThresh;		/* minimum proportion covered for cell to be included */
	char	filePropFull[_MAX_STATIC_BUFF_LEN];
	char	fileRelInf[_MAX_STATIC_BUFF_LEN];
	char	fileRelSus[_MAX_STATIC_BUFF_LEN];
	char	fileRelPri[_MAX_STATIC_BUFF_LEN];
	int		numIts;			/* number of iterations of the simulation to run */
	char	outStub[_MAX_STATIC_BUFF_LEN];
	double	ratePriInf;		/* this is the max rate at which expect primary infections over entire landscape */
	double	rateSecInf;		/* this is the secondary infection rate */
	double	dispScale;		/* average dispersal scale (measured in cells) */
	double	reportTime;		/* how frequently to report to the screen */
	double	maxTime;		/* maximum time to run epidemics up to (only used if don't set incidence) */
	double	maxIncidence;	/* this gives an alternate stopping condition */
	double	withinCellBulkUp;	/* logistic rate of increase of within-cell infection */
	double	withinCellMin;		/* minimum fraction of cell that can be infected (note if trueMinFlag=0 this is relative to carrying capacity of cell) */
	int		trueMinFlag;	/* whether to make withinCellMin relative to amount of hosts in cell (trueMinFlag=0) or a raw proportion (trueMinFlag=1) */
} t_Params;

/*
	This structure stores the probabilities of dispersing to cells in a single quadrant, from (0,0) -> (nCols, nRows) in a flattened list
	note that when doing dispersal, program decides whether to multiply both x and y by -1
*/
typedef struct
{
	int		numProbs;		/* number of probabilities that are stored */
	double  *aProbs;		/* array of dispersal probabilities */
	double	inCell;			/* probability of dispersing back to original cell (=aProbs[0]) */
	double	onLandscape;	/* probability of dispersing on the landscape */
} t_Dispersal;

/*
	Keep track of how many of each type of event were attempted
*/
typedef struct
{
	long	numSecondaryAttempts;
	long	numNonEmpty;
	long	numNonInfected;
	long	numSuccessful;
	long	numFindNextSecondary;
} t_RunStats;

/*
	Utility functions for the priority queue
*/
int getLeftChildIndex(int nodeIndex)
{
	return 2 * nodeIndex + 1;
}

int getRightChildIndex(int nodeIndex)
{
	return 2 * nodeIndex + 2;
}

int getParentIndex(int nodeIndex)
{
	return (nodeIndex - 1) / 2;
}

void siftUp(int nodeIndex, t_Epidemic *pEpidemic, t_Landscape *pLandscape)
{
      int parentIndex, tmp;

      if (nodeIndex != 0)
	  {
            parentIndex = getParentIndex(nodeIndex);
            if (pLandscape->aCells[pEpidemic->aQueueCells[parentIndex]].tNext >  pLandscape->aCells[pEpidemic->aQueueCells[nodeIndex]].tNext)
			{
                  tmp = pEpidemic->aQueueCells[parentIndex];
                  pEpidemic->aQueueCells[parentIndex] = pEpidemic->aQueueCells[nodeIndex];
                  pEpidemic->aQueueCells[nodeIndex] = tmp;
                  siftUp(parentIndex, pEpidemic, pLandscape);
            }
      }
}

void siftDown(int nodeIndex, t_Epidemic *pEpidemic, t_Landscape *pLandscape)
{
      int leftChildIndex, rightChildIndex, minIndex, tmp;

      leftChildIndex = getLeftChildIndex(nodeIndex);
      rightChildIndex = getRightChildIndex(nodeIndex);
      if (rightChildIndex >= pEpidemic->queueLen)
	  {
            if (leftChildIndex >= pEpidemic->queueLen)
			{
                  return;
			}
            else
			{
                  minIndex = leftChildIndex;
			}
      }
	  else
	  {
            if (pLandscape->aCells[pEpidemic->aQueueCells[leftChildIndex]].tNext <=  pLandscape->aCells[pEpidemic->aQueueCells[rightChildIndex]].tNext)
			{
                  minIndex = leftChildIndex;
			}
            else
			{
                  minIndex = rightChildIndex;
			}
      }
	  if (pLandscape->aCells[pEpidemic->aQueueCells[nodeIndex]].tNext >  pLandscape->aCells[pEpidemic->aQueueCells[minIndex]].tNext)
	  {
            tmp = pEpidemic->aQueueCells[minIndex];
            pEpidemic->aQueueCells[minIndex] = pEpidemic->aQueueCells[nodeIndex];
            pEpidemic->aQueueCells[nodeIndex] = tmp;
            siftDown(minIndex, pEpidemic, pLandscape);
      }
}

double removeMinElement(t_Epidemic *pEpidemic, t_Landscape *pLandscape)
{
	double min;

	min = -1;
	if (pEpidemic->queueLen == 0)
	{
		fprintf(stderr, "Heap is empty");
	}
	else
	{
		min = pLandscape->aCells[pEpidemic->aQueueCells[0]].tNext;
		pEpidemic->aQueueCells[0] = pEpidemic->aQueueCells[pEpidemic->queueLen - 1];
		pEpidemic->queueLen--;
		if (pEpidemic->queueLen > 0)
		{
				siftDown(0, pEpidemic, pLandscape);
		}
	}
	return min;
}

double removeArbitaryElement(int index, t_Epidemic *pEpidemic, t_Landscape *pLandscape)
{
	double thisElement;

	thisElement = -1.0;
	if(index >= pEpidemic->queueLen)
	{
		fprintf(stderr, "Trying to remove element that doesn't exist");
	}
	else
	{
		thisElement = pLandscape->aCells[pEpidemic->aQueueCells[index]].tNext;
		pEpidemic->aQueueCells[index] = pEpidemic->aQueueCells[pEpidemic->queueLen - 1];
		pEpidemic->queueLen--;
		if(pEpidemic->queueLen > 0)
		{
			siftDown(index, pEpidemic, pLandscape);
			siftUp(index, pEpidemic, pLandscape);
		}
	}
	return thisElement;
}

void insertElement(int cellIndex, t_Epidemic *pEpidemic, t_Landscape *pLandscape)
{
      if (pEpidemic->queueLen == pLandscape->numCells)
	  {
            fprintf(stderr, "Heap's storage has overflowed");
	  }
      else
	  {
            pEpidemic->queueLen++;
			pEpidemic->aQueueCells[pEpidemic->queueLen - 1] = cellIndex;
            siftUp(pEpidemic->queueLen - 1, pEpidemic, pLandscape);
      }
}

void checkHeap(FILE *heapOut, t_Epidemic *pEpidemic, t_Landscape *pLandscape)
{
	int i,c;

	/*
		In a heap, a node's value should be smaller than (or equal to) its children
	*/
	for(i=0;i<pEpidemic->queueLen;i++)
	{
		fprintf(heapOut, "%d -> %d -> %.5f\n", i, pEpidemic->aQueueCells[i], pLandscape->aCells[pEpidemic->aQueueCells[i]].tNext);
		c = getRightChildIndex(i);
		if(c >= pEpidemic->queueLen)
		{
			fprintf(heapOut, "\tright child: empty\n");
		}
		else
		{
			fprintf(heapOut, "\tright child: (%d %d %.5f) ", c, pEpidemic->aQueueCells[c], pLandscape->aCells[pEpidemic->aQueueCells[c]].tNext);
			if(pLandscape->aCells[pEpidemic->aQueueCells[c]].tNext >= pLandscape->aCells[pEpidemic->aQueueCells[i]].tNext)
			{
				fprintf(heapOut, "ok\n");
			}
			else
			{
				fprintf(heapOut, "failed\n");
			}
		}
		c = getLeftChildIndex(i);
		if(c >= pEpidemic->queueLen)
		{
			fprintf(heapOut, "\tleft child: empty\n");
		}
		else
		{
			fprintf(heapOut, "\tleft child: (%d %d %.5f) ", c, pEpidemic->aQueueCells[c], pLandscape->aCells[pEpidemic->aQueueCells[c]].tNext);
			if(pLandscape->aCells[pEpidemic->aQueueCells[c]].tNext >= pLandscape->aCells[pEpidemic->aQueueCells[i]].tNext)
			{
				fprintf(heapOut, "ok\n");
			}
			else
			{
				fprintf(heapOut, "failed\n");
			}
		}
	}
}

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
		if((pPtr = strrchr(szProgName,C_DIR_DELIMITER))!=NULL)
		{
			strcpy(szCfgFile,pPtr+1);
		}
		else
		{
			strcpy(szCfgFile,szProgName);
		}
		if((pPtr = strstr(szCfgFile,".exe"))!=NULL)
		{
			*pPtr = '\0';
		}
		strcat(szCfgFile,".cfg");
	}
	/* check file exists */
	fp = fopen(szCfgFile, "rb");
	if(fp)
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
	int	 bRet,i;
	FILE *fp;
	char *pThisPair;
	char *szArgvCopy;

	bRet = 0;
	i=0;
	/* try to find the relevant key on the command line */
	while(bRet == 0 && i<argc)
	{
		szArgvCopy = strdup(argv[i]);
		if(szArgvCopy)
		{
			pThisPair = strtok(szArgvCopy, " \t");
			while(pThisPair)
			{
				if(strncmp(pThisPair,szKey,strlen(szKey))==0)
				{
					pVal = strchr(pThisPair,'=');
					if(pVal)
					{
						/* make sure isn't just start of string matching the key */
						if(pThisPair[strlen(szKey)] == '=')
						{
							strcpy(szValue,pVal+1);
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
	if(bRet == 0)
	{
		fp = fopen(szCfgFile,"rb");
		if(fp)
		{
			char szLine[_MAX_STATIC_BUFF_LEN];

			while(!bRet && fgets(szLine,_MAX_STATIC_BUFF_LEN,fp))
			{
				char *pPtr;
				if((pPtr = strchr(szLine,'='))!=NULL)
				{
					*pPtr = '\0';
					if(strcmp(szKey,szLine)==0)
					{
						strcpy(szValue,pPtr+1);
						/* strip off newline (if any) */
						if((pPtr = strpbrk(szValue,"\r\n"))!=NULL)
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
	return(findKey(argc, argv, szCfgFile,szKey,szValue));
}

int readDoubleFromCfg(int argc, char **argv, char *szCfgFile, char *szKey, double *pdValue)
{
	char szValue[_MAX_STATIC_BUFF_LEN];

	if(findKey(argc, argv, szCfgFile,szKey,szValue))
	{
		*pdValue = atof(szValue);
		return 1;
	}
	return 0;
}

int readIntFromCfg(int argc, char **argv, char *szCfgFile, char *szKey, int *pnValue)
{
	char szValue[_MAX_STATIC_BUFF_LEN];

	if(findKey(argc, argv, szCfgFile,szKey,szValue))
	{
		*pnValue = atoi(szValue);
		return 1;
	}
	return 0;
}

/*
	Internally the two dimensional landscape is stored in a flattened structure
	These routines convert between the two
*/

int	gridToPos(int x, int y, int numCols)
{
	return x + y*numCols;
}

void posToGrid(int pos, int numCols, int *pX, int *pY)
{
	*pY = pos/numCols;
	*pX = pos%numCols;
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
	Seed random number generator

	Encapsulated to allow easy replacement if necessary
*/
void	seedRandom()
{
	unsigned long		ulnSeed;
	unsigned long		myPID;

	ulnSeed =(unsigned long) time(NULL);
	/* Make sure that multiple different processes started at same time have different seeds */
#ifndef _MSC_VER
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

/*
	Read all parameters from the configuration file
	(and dump to a new file so can tell which values program used)
*/

int readParams(t_Params *pParams, int argc, char **argv)
{
	char szCfgFile[_MAX_STATIC_BUFF_LEN];

	fprintf(stdout, "readParams()\n");
	memset(pParams,0,sizeof(t_Params));
	if(!getCfgFileName(argv[0], szCfgFile))
	{
		fprintf(stderr, "Couldn't find cfg file for program name '%s'\n", argv[0]);
		return 0;
	}
	if(!readDoubleFromCfg(argc, argv, szCfgFile, "cellThresh", &pParams->cellThresh))
	{
		fprintf(stdout, "Couldn't read cellThresh\n");
		return 0;
	}
	if(!readIntFromCfg(argc, argv, szCfgFile, "numIts", &pParams->numIts))
	{
		fprintf(stdout, "Couldn't read numIts\n");
		return 0;
	}
	if(!readStringFromCfg(argc, argv, szCfgFile, "filePropFull", pParams->filePropFull))
	{
		fprintf(stdout, "Couldn't read filePropFull\n");
		return 0;
	}
	if(!readStringFromCfg(argc, argv, szCfgFile, "fileRelInf", pParams->fileRelInf))
	{
		fprintf(stdout, "Couldn't read fileRelInf\n");
		return 0;
	}
	if(!readStringFromCfg(argc, argv, szCfgFile, "fileRelSus", pParams->fileRelSus))
	{
		fprintf(stdout, "Couldn't read fileRelSus\n");
		return 0;
	}
	if(!readStringFromCfg(argc, argv, szCfgFile, "fileRelPri", pParams->fileRelPri))
	{
		fprintf(stdout, "Couldn't read fileRelPri\n");
		return 0;
	}
	if(!readStringFromCfg(argc, argv, szCfgFile, "outStub", pParams->outStub))
	{
		fprintf(stdout, "Couldn't read outStub\n");
		return 0;
	}
#ifdef _MSC_VER
	if(!mkdir(pParams->outStub))
#else
	if(!mkdir(pParams->outStub, 0777))
#endif
	{
		fprintf(stdout, "\tcreated directory %s for output\n",pParams->outStub);
	}
	else
	{
		fprintf(stdout, "\tdirectory %s already exists\n",pParams->outStub);
	}
	if(!readDoubleFromCfg(argc, argv, szCfgFile, "maxTime", &pParams->maxTime))
	{
		fprintf(stdout, "Couldn't read maxTime\n");
		return 0;
	}
	if(!readDoubleFromCfg(argc, argv, szCfgFile, "ratePriInf", &pParams->ratePriInf))
	{
		fprintf(stdout, "Couldn't read ratePriInf\n");
		return 0;
	}
	if(!readDoubleFromCfg(argc, argv, szCfgFile, "rateSecInf", &pParams->rateSecInf))
	{
		fprintf(stdout, "Couldn't read rateSecInf\n");
		return 0;
	}
	if(!readDoubleFromCfg(argc, argv, szCfgFile, "dispScale", &pParams->dispScale))
	{
		fprintf(stdout, "Couldn't read dispScale\n");
		return 0;
	}
	if(!readDoubleFromCfg(argc, argv, szCfgFile, "reportTime", &pParams->reportTime))
	{
		fprintf(stdout, "Couldn't read reportTime\n");
		return 0;
	}
	if(!readDoubleFromCfg(argc, argv, szCfgFile, "maxIncidence", &pParams->maxIncidence))
	{
		fprintf(stdout, "Couldn't read maxIncidence\n");
		return 0;
	}
	if(!readDoubleFromCfg(argc, argv, szCfgFile, "withinCellBulkUp", &pParams->withinCellBulkUp))
	{
		fprintf(stdout, "Couldn't read withinCellBulkUp\n");
		return 0;
	}
	if(!readDoubleFromCfg(argc, argv, szCfgFile, "withinCellMin", &pParams->withinCellMin))
	{
		fprintf(stdout, "Couldn't read withinCellMin\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "trueMinFlag", &pParams->trueMinFlag))
	{
		fprintf(stdout, "Couldn't read trueMinFlag\n");
		return 0;
	}
	/* Following code dumps out parameters as read in */
	{
		FILE	*paramsOut;
		char	outFile[_MAX_STATIC_BUFF_LEN];

		sprintf(outFile, "%s%cparamsOut.txt", pParams->outStub, C_DIR_DELIMITER);
		fprintf(stdout, "\twriting copy of parameters to %s\n", outFile);
		paramsOut = fopen(outFile, "wb");
		if(paramsOut)
		{
			fprintf(paramsOut, "pParams->cellThresh=%.6f\n", pParams->cellThresh);
			fprintf(paramsOut, "pParams->numIts=%d\n", pParams->numIts);
			fprintf(paramsOut, "pParams->filePropFull=%s\n", pParams->filePropFull);
			fprintf(paramsOut, "pParams->fileRelInf=%s\n", pParams->fileRelInf);
			fprintf(paramsOut, "pParams->fileRelSus=%s\n", pParams->fileRelSus);
			fprintf(paramsOut, "pParams->fileRelPri=%s\n", pParams->fileRelPri);
			fprintf(paramsOut, "pParams->outStub=%s\n", pParams->outStub);
			fprintf(paramsOut, "pParams->maxTime=%.6f\n", pParams->maxTime);
			fprintf(paramsOut, "pParams->ratePriInf=%.6f\n", pParams->ratePriInf);
			fprintf(paramsOut, "pParams->rateSecInf=%.6f\n", pParams->rateSecInf);
			fprintf(paramsOut, "pParams->dispScale=%.6f\n", pParams->dispScale);
			fprintf(paramsOut, "pParams->reportTime=%f\n", pParams->reportTime);
			fprintf(paramsOut, "pParams->maxIncidence=%.6f\n", pParams->maxIncidence);
			fprintf(paramsOut, "pParams->withinCellMin=%.6f\n", pParams->withinCellMin);
			fprintf(paramsOut, "pParams->withinCellBulkUp=%.6f\n", pParams->withinCellBulkUp);
			fprintf(paramsOut, "pParams->trueMinFlag=%d\n", pParams->trueMinFlag);
			fclose(paramsOut);
		}
		else
		{
			fprintf(stderr, "couldn't dump parameters...exiting\n");
			return 0;
		}
	}
	return 1;
}

/*
	Rather unwieldy parsing routine which reads in all data from GIS format
*/
int readLandscape(t_Landscape *pLandscape, char *filePropFull, char *fileRelInf, char *fileRelPri, char *filRelSus, double cellThresh, char *outStub)
{
	FILE	*fIn;
	char	*fileName,*inBuff,*pPtr;
	int		noData,retVal,thisX,thisY,i,headerCount;
	double	thisVal;
	char	outFile[_MAX_STATIC_BUFF_LEN];
	t_Cell* pTmp;

	fprintf(stdout, "readLandscape()\n");
	retVal = 0;
	memset(pLandscape,0,sizeof(t_Landscape));
	inBuff = malloc(_MAX_DYNAMIC_BUFF_LEN*sizeof(char));
	if(inBuff)
	{
		retVal = 1;
		for(i=0;i<4 && retVal;i++)
		{
			retVal = 0;
			fileName = NULL;
			switch(i)
			{
			case 0:
				fileName = filePropFull;
				break;
			case 1:
				fileName = fileRelInf;
				break;
			case 2:
				fileName = fileRelPri;
				break;
			case 3:
				fileName = filRelSus;
				break;
			default:
				fprintf(stderr, "readLandscape(): shouldn't get here...\n");
				break;
			}
			fprintf(stdout, "\t%s...", fileName);
			fIn = fopen(fileName,"rb");
			if(fIn)
			{
				/* strip GIS header */
				retVal = 1;
				headerCount = 0;
				while(retVal && headerCount < _GIS_HEADER_LENGTH && fgets(inBuff, _MAX_DYNAMIC_BUFF_LEN, fIn))
				{
					if(i == 0)
					{
						if(strncmp(inBuff, "ncols", strlen("ncols"))==0)
						{
							pPtr = strpbrk(inBuff,_GIS_WHITESPACE);
							if(pPtr)
							{
								pLandscape->numCols = atoi(pPtr);
							}
							else
							{
								retVal = 0;
							}
						}
						if(strncmp(inBuff, "nrows", strlen("nrows"))==0)
						{
							pPtr = strpbrk(inBuff,_GIS_WHITESPACE);
							if(pPtr)
							{
								pLandscape->numRows = atoi(pPtr);
							}
							else
							{
								retVal = 0;
							}
						}
					}
					headerCount++;
				}
				if(retVal && pLandscape->numCols > 0 && pLandscape->numRows > 0)
				{
					if(i == 0)
					{
						pLandscape->aCellLookup = malloc(sizeof(int)*pLandscape->numCols*pLandscape->numRows);
						if(!pLandscape->aCellLookup)
						{
							retVal = 0;
						}
					}
					thisY = 0;
					while(retVal && fgets(inBuff, _MAX_DYNAMIC_BUFF_LEN, fIn)) 	/* will fail on very long lines */
					{
						pPtr = strpbrk(inBuff, _GIS_NEWLINES);
						if(pPtr)
						{
							pPtr[0] = '\0';
						}
						else
						{
							fprintf(stderr, "readLandscape(): line too long\n");
							retVal=0;											/* no CR or NL means line not read in fully */
						}
						thisX = 0;
						pPtr = strtok(inBuff, _GIS_WHITESPACE);
						while(retVal && pPtr)
						{
							noData = 0;
							if(strcmp(pPtr,_GIS_NODATA)==0)
							{
								noData = 1;
							}
							else
							{
								thisVal = atof(pPtr);
							}
							if(i==0)	/* in the first pass through need to add this cell to the landscape */
							{
								if(noData || thisVal < cellThresh)
								{
									pLandscape->aCellLookup[gridToPos(thisX,thisY,pLandscape->numCols)] = _EMPTY_CELL;
								}
								else
								{
									if(pLandscape->numCells == pLandscape->numCellsSpace)
									{
										pLandscape->numCellsSpace += _LANDSCAPE_BLOCK_SIZE;

										/* have to rewrite perfectly idiomatic C to get around visual studio warnings re. memory leaks */
										pTmp = realloc(pLandscape->aCells, sizeof(t_Cell) * pLandscape->numCellsSpace);
										if (pTmp)
										{
											pLandscape->aCells = pTmp;
										}
										else
										{
											retVal = 0;
										}
									}
									pLandscape->aCellLookup[gridToPos(thisX,thisY,pLandscape->numCols)] = pLandscape->numCells;
									if (pLandscape->aCells)
									{
										pLandscape->aCells[pLandscape->numCells].xPos = thisX;
										pLandscape->aCells[pLandscape->numCells].yPos = thisY;
										pLandscape->aCells[pLandscape->numCells].propFull = thisVal;
										pLandscape->aCells[pLandscape->numCells].tInf = _UNDEF_TIME;
										pLandscape->aCells[pLandscape->numCells].tNext = _UNDEF_TIME;
										pLandscape->aCells[pLandscape->numCells].infBy = _EMPTY_CELL;
										pLandscape->aCells[pLandscape->numCells].infType = _EMPTY_CELL;
										pLandscape->totalFull += pLandscape->aCells[pLandscape->numCells].propFull;
										pLandscape->numCells++;
									}
								}
							}
							else
							{
								if(pLandscape->aCellLookup[gridToPos(thisX,thisY,pLandscape->numCols)] != _EMPTY_CELL)
								{
									if(!noData)
									{
										switch(i)
										{
										case 1:
											pLandscape->aCells[pLandscape->aCellLookup[gridToPos(thisX,thisY,pLandscape->numCols)]].relInf = thisVal;
											break;
										case 2:
											pLandscape->aCells[pLandscape->aCellLookup[gridToPos(thisX,thisY,pLandscape->numCols)]].relPri = thisVal;
											break;
										case 3:
											pLandscape->aCells[pLandscape->aCellLookup[gridToPos(thisX,thisY,pLandscape->numCols)]].relSus = thisVal;
											break;
										default:
											fprintf(stderr, "readLandscape(): shouldn't get here...\n");
											break;
										}
									}
									else
									{
										fprintf(stderr, "readLandscape(): NODATA when expecting value...\n");
										retVal = 0;
									}
								}
							}
							thisX++;
							pPtr = strtok(NULL, _GIS_WHITESPACE);
						}
						if(thisX != pLandscape->numCols)
						{
							fprintf(stderr, "%s, line %d: bad number of columns (%d)\n", fileName, thisY, thisX);
							retVal = 0;
						}
						thisY++;
					}
					if(thisY != pLandscape->numRows)
					{
						fprintf(stderr, "%s: bad number of rows (%d)\n", fileName, thisY);
						retVal = 0;
					}
				}
				else
				{
					fprintf(stderr, "readLandscape(): failed to parse gis header\n");
				}
				if(retVal)
				{
					fprintf(stdout, "successfully\n");
				}
				fclose(fIn);
			}
			else
			{
				fprintf(stderr, "couldn't open %s\n", fileName);
				retVal = 0;
				break;
			}
		}
		free(inBuff);
	}
	else
	{
		fprintf(stderr, "out of memory\n");
		retVal = 0;
	}
	/*
		Write out information on all cells that are active in the simulation (i.e. >= cellThresh)
	*/
	if(retVal)
	{
		FILE *fpOut;
		int	 i;

		sprintf(outFile, "%s%cactiveLandscape.txt", outStub, C_DIR_DELIMITER);
		fpOut = fopen(outFile, "wb");
		if(fpOut)
		{
			for(i=0;i<pLandscape->numCells;i++)
			{
				fprintf(fpOut, "%d %d %f %d\n", pLandscape->aCells[i].xPos, pLandscape->aCells[i].yPos, pLandscape->aCells[i].propFull, i);
			}
			fclose(fpOut);
		}
		fprintf(stdout, "\t%d valid cells (%d rows; %d cols)\n", pLandscape->numCells, pLandscape->numRows, pLandscape->numCols);
		fprintf(stdout, "\ttotalFull=%f\n", pLandscape->totalFull);
	}
	return retVal;
}

/*
	Primary rate of infection is stored in cumulative form to make it easier to find out which cell is infected next
*/
int setupPrimary(t_PriInf *pPriInf, t_Landscape *pLandscape, double ratePri)
{
	int		i,retVal;
	double	thisVal,cumVal;

	fprintf(stdout, "setupPrimary()\n");
	retVal = 0;
	memset(pPriInf,0,sizeof(t_PriInf));
	pPriInf->aCumPressure = malloc(sizeof(double)*pLandscape->numCells);
	if(pPriInf->aCumPressure)
	{
		cumVal = 0;
		for(i=0;i<pLandscape->numCells;i++)
		{
			/*
				Values in the GIS file are multipled by proportion of cell occupied, as well as relative susceptibility
			*/
			thisVal = pLandscape->aCells[i].propFull * pLandscape->aCells[i].relPri * pLandscape->aCells[i].relSus;
			cumVal += thisVal;
			pPriInf->aCumPressure[i] = cumVal;
		}
		pPriInf->totalPressure = cumVal;
		pPriInf->nextT = _UNDEF_TIME;
		pPriInf->ratePri = ratePri;
		retVal = 1;
		fprintf(stdout, "\ttotalPressure=%f\n", pPriInf->totalPressure);
	}
	return retVal;
}

/*
	Allocate memory for epidemic
*/
int setupEpidemic(t_Epidemic *pEpidemic, t_Landscape *pLandscape, double dispScale)
{
	int retVal;

	fprintf(stdout, "setupEpidemic()\n");
	retVal = 0;
	memset(pEpidemic,0,sizeof(t_Epidemic));
	pEpidemic->queueLen = 0;
	pEpidemic->aQueueCells = malloc(sizeof(int) * pLandscape->numCells);
	pEpidemic->totalInf = 0;
	pEpidemic->aInfCells = malloc(sizeof(int) * pLandscape->numCells);
	if(pEpidemic->aQueueCells && pEpidemic->aInfCells)
	{
		retVal = 1;
		fprintf(stdout, "\t%d cells\n", pLandscape->numCells);
	}
	return retVal;
}

/*
	Find the time of the next primary infection across the landscape
		given that the total rate of entry on an entirely susceptible landscape is pPriInf->ratePri
		per unit of time [that some will hit already infected cells is handled later]
*/
void setNextPossPriTime(t_PriInf *pPriInf, double thisTime)
{
	double randDbl;

	randDbl = uniformRandom();
	if(pPriInf->ratePri > 0.0)
	{
		pPriInf->nextT = thisTime - log(randDbl)/pPriInf->ratePri;
	}
	else
	{
		pPriInf->nextT=_VERY_LONG_TIME;
	}
}

double getNextPossPriTime(t_PriInf *pPriInf)
{
	return pPriInf->nextT;
}

/*
	Binary chop to figure out cell which is infected by primary infection
*/
int		whichCellPrimary(t_PriInf *pPriInf,int numCells)
{
	double  randDbl;//,runningSum;
	int		left,right,mid;//,cellID;

	randDbl = pPriInf->totalPressure * uniformRandom();
	left = 0;
	right = numCells-1;
	while((right-left) > 1)
	{
		mid=(left+right)/2;
		if(pPriInf->aCumPressure[mid] < randDbl)
		{
			left = mid;
		}
		else
		{
			right = mid;
		}
	}
	return right;
}

/*
	Figure out which cell is challenged by a potential secondary infection
*/
int whichCellSecondary(t_Dispersal *pDispersal, t_Landscape *pLandscape, int cellInfectFrom, t_Params *pParams)
{
	double	randDbl;
	int		posToChallenge,x,y,xOffset,yOffset,left,right,mid,cellToChallenge,cellQuad;

	cellToChallenge = _EMPTY_CELL;
	/*
		First need to find cell to challenge
	*/
	randDbl = 4.0 * uniformRandom();
	cellQuad = (int)randDbl;
	randDbl -= cellQuad;
	if(!(randDbl < pDispersal->inCell || randDbl > pDispersal->onLandscape))
	{
		left = 0;
		right = pDispersal->numProbs-1;
		while((right-left) > 1)
		{
			mid=(left+right)/2;
			if(pDispersal->aProbs[mid] < randDbl)
			{
				left = mid;
			}
			else
			{
				right = mid;
			}
		}
		/* right is the cell to challenge*/;
		posToGrid(right, pLandscape->numCols, &xOffset, &yOffset);
		/* but need to account for only storing one quarter of the kernel */
		switch(cellQuad)
		{
		case 0:
			/* do nothing */
			break;
		case 1:
			xOffset *= -1;
			break;
		case 2:
			xOffset *= -1;
			yOffset *= -1;
			break;
		case 3:
			yOffset *= -1;
			break;
		default:
			fprintf(stderr, "should never get here...\n");
			break;
		}
		x = pLandscape->aCells[cellInfectFrom].xPos + xOffset;
		if(x >= 0 && x < pLandscape->numCols)
		{
			y = pLandscape->aCells[cellInfectFrom].yPos + yOffset;
			if(y >=0 && y < pLandscape->numRows)
			{
				posToChallenge = gridToPos(x, y, pLandscape->numCols);
				cellToChallenge = pLandscape->aCellLookup[posToChallenge];
				if(cellToChallenge == _EMPTY_CELL)
				{
#ifdef _DEBUG_PRINT_MSG
					fprintf(stdout, "\t\t\t\tno hosts\n");
#endif
				}
			}
			else
			{
#ifdef _DEBUG_PRINT_MSG
				fprintf(stdout, "\t\t\t\toff landscape\n");
#endif
			}
		}
		else
		{
#ifdef _DEBUG_PRINT_MSG
			fprintf(stdout, "\t\t\t\toff landscape\n");
#endif
		}
	}
	else
	{
		if(randDbl > pDispersal->onLandscape)
		{
#ifdef _DEBUG_PRINT_MSG
			fprintf(stdout, "\t\t\t\toff landscape\n");
#endif
		}
		else
		{
#ifdef _DEBUG_PRINT_MSG
			fprintf(stdout, "\t\t\t\twithin cell\n");
#endif
		}
	}
	return cellToChallenge;
}

/*
	Peek at minimum element in the priority queue to find time of next secondary infection
*/
double getNextPossSecTime(t_Epidemic *pEpidemic, t_Landscape *pLandscape, double tPrimary)
{
	if(pEpidemic->queueLen > 0)
	{
		/* Peek at min element of queue (do not remove it) */
		return pLandscape->aCells[pEpidemic->aQueueCells[0]].tNext;
	}
	/*
		If nothing in the queue (i.e. no primary infection has happened) make sure
		the putative secondary comes after the next primary (and so won't be selected anyway)
	*/
	return tPrimary + 1.0;
}

/*
	Retrieve element from the priority queue
*/
int getCellInfectFrom(t_Epidemic *pEpidemic, t_Landscape *pLandscape)
{
	int minCell;

	if(pEpidemic->queueLen > 0)
	{
		minCell = pEpidemic->aQueueCells[0];
		removeMinElement(pEpidemic, pLandscape);
		return minCell;
	}
	return _EMPTY_CELL;
}

/*
	Find the time of the next secondary infection from a given cell
*/
void findNextSecondary(t_Landscape *pLandscape, int thisCell, double thisTime, double rateSecInf, t_Epidemic *pEpidemic, t_RunStats *pRunStats, double withinCellMin, double withinCellBulkUp, int trueMinFlag)
{
	double	randDbl,rateSec;
	double	deltaMin;			/* delay before next infection if cell were full of infection */
	double	tSinceInf;			/* time since this cell was infected */
	double	deltaReal;			/* delay before next infection accounting for logistic bulk up */
	double	logisticJ;			/* J= (1-withinCellMin)/withinCellMin */
	double  thisWCM;

	pRunStats->numFindNextSecondary++;
	randDbl = uniformRandom();
	/* find maximum rate of infection from this cell */
	rateSec = pLandscape->aCells[thisCell].propFull*pLandscape->aCells[thisCell].relInf*rateSecInf;
	if(rateSec > 0)
	{
		/* lengthen length of time until the infection to account for infectivity bulking up logistically */
		deltaMin = -log(randDbl)/rateSec;
		tSinceInf = thisTime - pLandscape->aCells[thisCell].tInf;
		if (trueMinFlag)
		{
			if (withinCellMin >= pLandscape->aCells[thisCell].propFull)
			{
				/* make it fully infected - i.e. incidence = pLandscape->aCells[thisCell].propFull - immediately */
				logisticJ = 0;
			}
			else
			{
				/* figure out the correct fraction initially infected to make the initial incidence = withinCellMin */
				thisWCM = withinCellMin / pLandscape->aCells[thisCell].propFull;
				logisticJ = (1.0 - thisWCM) / thisWCM;
			}
		}
		else
		{
			logisticJ = (1.0 - withinCellMin) / withinCellMin;
		}
		deltaReal = (1.0/withinCellBulkUp)*(log(exp(withinCellBulkUp*(tSinceInf+deltaMin)) + logisticJ*(exp(withinCellBulkUp*deltaMin)-1.0))) - tSinceInf;
		pLandscape->aCells[thisCell].tNext = thisTime + deltaReal;
		/* Add this entry to queue */
		insertElement(thisCell, pEpidemic, pLandscape);
	}
}

/*
	Book-keeping to handle a cell newly becoming infected
*/
void infectCell(t_Landscape *pLandscape, int thisCell, double thisTime, double rateSecInf, t_Epidemic *pEpidemic, int infType, int infBy, t_RunStats *pRunStats, double withinCellMin, double withinCellBulkUp, int trueMinFlag)
{
	/* set the time of infection */
	pLandscape->aCells[thisCell].tInf = thisTime;
	pLandscape->aCells[thisCell].infType = infType;
	pLandscape->aCells[thisCell].infBy = infBy;
	/* add to the list of all infections for later dumping */
	pEpidemic->aInfCells[pEpidemic->totalInf] = thisCell;
	pEpidemic->totalInf++;
	/* find time of next (potential) secondary infection from this cell */
	findNextSecondary(pLandscape, thisCell, thisTime, rateSecInf, pEpidemic, pRunStats, withinCellMin, withinCellBulkUp, trueMinFlag);
}

/*
	Find the incidence in a given cell based on the current time and its time of infection
*/
double getIncidence(t_Landscape *pLandscape, t_Epidemic *pEpidemic, double thisTime, int hostID, double withinCellMin, double withinCellBulkUp, int trueMinFlag)
{
	double tSinceInf,logisticJ,thisIncidence,thisWCM;

	tSinceInf = thisTime - pLandscape->aCells[hostID].tInf;
	if (tSinceInf >= 0.0)	/* getIncidence() can be called before a cell has become infected...if so ignore */
	{
		if (trueMinFlag)
		{
			if (withinCellMin >= pLandscape->aCells[hostID].propFull)
			{
				logisticJ = 0.0;
			}
			else
			{
				thisWCM = withinCellMin / pLandscape->aCells[hostID].propFull;
				logisticJ = (1.0 - thisWCM) / thisWCM;
			}
		}
		else
		{
			logisticJ = (1.0 - withinCellMin) / withinCellMin;
		}
		thisIncidence = pLandscape->aCells[hostID].propFull / (1 + logisticJ * exp(-withinCellBulkUp*tSinceInf));
	}
	else
	{
		thisIncidence = 0.0;
	}
	return thisIncidence;
}

/*
	Main routine to run an ensemble of epidemics and dump the results
*/
int runEpidemics(t_Params *pParams, t_Landscape *pLandscape, t_PriInf *pPriInf, t_Dispersal *pDispersal, t_Epidemic *pEpidemic)
{
	int			doneInf,firstInf,continueRunning,j,i,k,retVal,cellToChallenge,cellInfectFrom,thisReason;
	double		thisTime,nextPri,nextSec,randDbl,infectProb,nextReport,trueIncidence,maxFullIncidence,thisFinalIncidence;
	char		outFile[_MAX_STATIC_BUFF_LEN],dpcFile[_MAX_STATIC_BUFF_LEN];
	FILE		*fOut,*fEnd,*fSingleEnd,*fDPC;
	t_RunStats	runStats;

	trueIncidence = 0.0;
	fprintf(stdout, "runEpidemics()\n");
	retVal = 1;
	sprintf(outFile, "%s%cendTimes.txt", pParams->outStub, C_DIR_DELIMITER);
	fEnd = fopen(outFile, "wb");
	if(!fEnd)
	{
		fprintf(stderr, "couldn't open endTimes file for writing\n");
		return 0;
	}
	for(i=0;i<pParams->numIts;i++)
	{
		thisReason = 0;		/* will be set to 1 if simulation stops because hit threshold incidence */
		sprintf(dpcFile, "%s%c%s_dpc_%d.txt", pParams->outStub, C_DIR_DELIMITER, pParams->outStub, i);
		fDPC = fopen(dpcFile, "wb");
		if(fDPC)
		{
			fprintf(stdout, "\titeration %d\n", i);
			memset(&runStats, 0, sizeof(t_RunStats));
			thisTime = 0.0;
			setNextPossPriTime(pPriInf, thisTime);
			nextReport = 0.0;
			if (pParams->ratePriInf == 0.0)
			{
				firstInf = (int)((double)pLandscape->numCells*uniformRandom());
				infectCell(pLandscape, firstInf, 0.0, pParams->rateSecInf, pEpidemic, _PRI_INF_TYPE, _EMPTY_CELL, &runStats, pParams->withinCellMin, pParams->withinCellBulkUp, pParams->trueMinFlag);
				fprintf(stdout, "infecting %d at t=0.0\n", firstInf);
			}
			continueRunning = 1;
			maxFullIncidence = pParams->maxIncidence * pLandscape->totalFull;
			while (retVal && continueRunning)
			{
				doneInf = 0;
				while (nextReport <= thisTime)
				{
					trueIncidence = 0.0;
					for (j = 0; j < pEpidemic->totalInf; j++)
					{
						trueIncidence += getIncidence(pLandscape, pEpidemic, nextReport, pEpidemic->aInfCells[j], pParams->withinCellMin, pParams->withinCellBulkUp, pParams->trueMinFlag);
					}
					fprintf(stdout, "\t\tt=%.4f (infNum=%d propInf=%.4f trueInc=%.4f)\n", nextReport, pEpidemic->totalInf, pEpidemic->totalInf / (double)pLandscape->numCells, trueIncidence / pLandscape->totalFull);
					fprintf(fDPC, "%.4f %d %.4f %.4f\n", nextReport, pEpidemic->totalInf, pEpidemic->totalInf / (double)pLandscape->numCells, trueIncidence / pLandscape->totalFull);
					nextReport += pParams->reportTime;
				}
				nextPri = getNextPossPriTime(pPriInf);
				nextSec = getNextPossSecTime(pEpidemic, pLandscape, nextPri);
				if (nextPri < nextSec)
				{
					/* attempt a primary infection */
					if (nextPri < pParams->maxTime)
					{
						/* update time */
						thisTime = nextPri;
						/* find the cell to challenge */
						cellToChallenge = whichCellPrimary(pPriInf, pLandscape->numCells);
#ifdef _DEBUG_PRINT_MSG
						fprintf(stdout, "\t\t(primary) challenging %d at %.4f\n", cellToChallenge, thisTime);
#endif
						/*
							Note that have built relSus and area into the rate of primary infection
							for each cell, so just need to check whether it is already infected or not
						*/
						if (pLandscape->aCells[cellToChallenge].tInf >= 0.0)					/* already infected */
						{
#ifdef _DEBUG_PRINT_MSG
							fprintf(stdout, "\t\t\talready infected\n");
#endif
						}
						else
						{
							/* infect */
#ifdef _DEBUG_PRINT_MSG
							fprintf(stdout, "\t\t\tinfecting\n");
#endif
							infectCell(pLandscape, cellToChallenge, thisTime, pParams->rateSecInf, pEpidemic, _PRI_INF_TYPE, _EMPTY_CELL, &runStats, pParams->withinCellMin, pParams->withinCellBulkUp, pParams->trueMinFlag);
							doneInf = 1;
						}
						setNextPossPriTime(pPriInf, thisTime);
					}
					else
					{
						thisTime = pParams->maxTime;
					}
				}
				else
				{
					/* attempt a secondary infection */
					if (nextSec < pParams->maxTime)
					{
						thisTime = nextSec;
						/* find a cell to challenge, and challenge it if makes sense to */
						cellInfectFrom = getCellInfectFrom(pEpidemic, pLandscape);
#ifdef _DEBUG_PRINT_MSG
						fprintf(stdout, "\t\t(secondary) challenging from %d at %.4f\n", cellInfectFrom, thisTime);
#endif
						if (cellInfectFrom != _EMPTY_CELL)
						{
							runStats.numSecondaryAttempts++;
							cellToChallenge = whichCellSecondary(pDispersal, pLandscape, cellInfectFrom, pParams);
#ifdef _DEBUG_PRINT_MSG
							fprintf(stdout, "\t\t\tchallenging %d\n", cellToChallenge);
#endif
							if (cellToChallenge != _EMPTY_CELL)
							{
								runStats.numNonEmpty++;
								if (pLandscape->aCells[cellToChallenge].tInf >= 0.0)					/* already infected */
								{
#ifdef _DEBUG_PRINT_MSG
									fprintf(stdout, "\t\t\talready infected\n");
#endif
								}
								else
								{
									runStats.numNonInfected++;
									/* possibly infect, depending on relative susceptibility */
									infectProb = pLandscape->aCells[cellToChallenge].relSus*pLandscape->aCells[cellToChallenge].propFull;
#ifdef _DEBUG_PRINT_MSG
									fprintf(stdout, "\t\t\tp(infect)=%f\n", infectProb);
#endif
									randDbl = uniformRandom();
									if (randDbl < infectProb)
									{
										runStats.numSuccessful++;
#ifdef _DEBUG_PRINT_MSG
										fprintf(stdout, "\t\t\t\tinfecting\n");
#endif
										infectCell(pLandscape, cellToChallenge, thisTime, pParams->rateSecInf, pEpidemic, _SEC_INF_TYPE, cellInfectFrom, &runStats, pParams->withinCellMin, pParams->withinCellBulkUp, pParams->trueMinFlag);
										doneInf = 1;
									}
									else
									{
#ifdef _DEBUG_PRINT_MSG
										fprintf(stdout, "\t\t\t\tfailed to infect\n");
#endif
									}
								}
							}
							/* need to update the source cell's time of next secondary infection too */
							findNextSecondary(pLandscape, cellInfectFrom, thisTime, pParams->rateSecInf, pEpidemic, &runStats, pParams->withinCellMin, pParams->withinCellBulkUp, pParams->trueMinFlag);
						}
						else
						{
							fprintf(stderr, "\tsecondary infection from invalid cell...\n");
							retVal = 0;
						}
					}
					else
					{
						thisTime = pParams->maxTime;
					}
				}
				if (doneInf)
				{
					trueIncidence = 0.0;
					for (j = 0; j < pEpidemic->totalInf; j++)
					{
						trueIncidence += getIncidence(pLandscape, pEpidemic, thisTime, pEpidemic->aInfCells[j], pParams->withinCellMin, pParams->withinCellBulkUp, pParams->trueMinFlag);
					}
				}
#if 0
				if (maxInfected > 0)
				{
					if (pEpidemic->totalInf >= maxInfected)
					{
						continueRunning = 0;
					}
				}
#endif
				if (pParams->maxIncidence > 0.0)
				{
					if (trueIncidence >= maxFullIncidence)
					{
						continueRunning = 0;
						thisReason = 1;
					}
				}
				if (thisTime >= pParams->maxTime)
				{
					continueRunning = 0;
				}
			}
#if _CHECK_HEAP
			{
				FILE *fpTmp = fopen("checkHeap.txt", "wb");

				if (fpTmp)
				{
					checkHeap(fpTmp, pEpidemic, pLandscape);
					fclose(fpTmp);
				}
			}
#endif
			/* Do a final round of printing to the screen */
			{
				trueIncidence = 0.0;
				for (j = 0; j < pEpidemic->totalInf; j++)
				{
					trueIncidence += getIncidence(pLandscape, pEpidemic, thisTime, pEpidemic->aInfCells[j], pParams->withinCellMin, pParams->withinCellBulkUp, pParams->trueMinFlag);
				}
				fprintf(stdout, "\t\tt=%.4f (infNum=%d propInf=%.4f trueInc=%.4f)\n", thisTime, pEpidemic->totalInf, pEpidemic->totalInf / (double)pLandscape->numCells, trueIncidence / pLandscape->totalFull);
				fprintf(fDPC, "%.4f %d %.4f %.4f\n", thisTime, pEpidemic->totalInf, pEpidemic->totalInf / (double)pLandscape->numCells, trueIncidence / pLandscape->totalFull);
			}
			fprintf(fEnd, "%f\n", thisTime);
			sprintf(outFile, "%s%cendTime_%d.txt", pParams->outStub, C_DIR_DELIMITER, i);
			fSingleEnd = fopen(outFile, "wb");
			if (!fSingleEnd)
			{
				fprintf(stderr, "couldn't open endTimes file for writing\n");
				return 0;
			}
			fprintf(fSingleEnd, "%f\n", thisTime);
			fclose(fSingleEnd);

			/*
				Dump reason simulation stopped to a file
			*/
			sprintf(outFile, "%s%cendReason_%d.txt", pParams->outStub, C_DIR_DELIMITER, i);
			fSingleEnd = fopen(outFile, "wb");
			if (!fSingleEnd)
			{
				fprintf(stderr, "couldn't open endReason file for writing\n");
				return 0;
			}
			fprintf(fSingleEnd, "%d\n", thisReason);
			fclose(fSingleEnd);


			/*
				Dump all the information
			*/
			sprintf(outFile, "%s%c%s_%d.txt", pParams->outStub, C_DIR_DELIMITER, pParams->outStub, i);
			fOut = fopen(outFile, "wb");
			if (fOut)
			{
				for (j = 0; j < pEpidemic->totalInf; j++)
				{
					thisFinalIncidence = getIncidence(pLandscape, pEpidemic, thisTime, pEpidemic->aInfCells[j], pParams->withinCellMin, pParams->withinCellBulkUp, pParams->trueMinFlag)/pLandscape->aCells[pEpidemic->aInfCells[j]].propFull;
					trueIncidence = 0.0;
					for (k = 0; k <= j; k++)
					{
						trueIncidence += getIncidence(pLandscape, pEpidemic, pLandscape->aCells[pEpidemic->aInfCells[j]].tInf, pEpidemic->aInfCells[k], pParams->withinCellMin, pParams->withinCellBulkUp, pParams->trueMinFlag);
					}
					fprintf(fOut, "%d %d %.4f %d %d %d %.4f %.4f %.4f %.4f %d %.4f %d %.4f %.4f\n",
						pLandscape->aCells[pEpidemic->aInfCells[j]].xPos,
						pLandscape->aCells[pEpidemic->aInfCells[j]].yPos,
						pLandscape->aCells[pEpidemic->aInfCells[j]].tInf,
						pLandscape->aCells[pEpidemic->aInfCells[j]].infType,
						(pLandscape->aCells[pEpidemic->aInfCells[j]].infBy == _EMPTY_CELL) ? _EMPTY_CELL : pLandscape->aCells[pLandscape->aCells[pEpidemic->aInfCells[j]].infBy].xPos,
						(pLandscape->aCells[pEpidemic->aInfCells[j]].infBy == _EMPTY_CELL) ? _EMPTY_CELL : pLandscape->aCells[pLandscape->aCells[pEpidemic->aInfCells[j]].infBy].yPos,
						pLandscape->aCells[pEpidemic->aInfCells[j]].propFull,
						pLandscape->aCells[pEpidemic->aInfCells[j]].relInf,
						pLandscape->aCells[pEpidemic->aInfCells[j]].relSus,
						pLandscape->aCells[pEpidemic->aInfCells[j]].relPri,
						(j + 1),
						(j + 1.0) / (double)pLandscape->numCells,
						pEpidemic->aInfCells[j],
						trueIncidence / pLandscape->totalFull,
						thisFinalIncidence);
				}
				fclose(fOut);
			}
			/*
				Blank all the information so start next simulation totally afresh
			*/
			pEpidemic->queueLen = 0;
			pEpidemic->totalInf = 0;
			for (j = 0; j < pLandscape->numCells; j++)
			{
				pLandscape->aCells[j].tInf = _UNDEF_TIME;
				pLandscape->aCells[j].tNext = _UNDEF_TIME;
				pLandscape->aCells[j].infBy = _EMPTY_CELL;
				pLandscape->aCells[j].infType = _EMPTY_CELL;
			}
			pPriInf->nextT = _UNDEF_TIME;
			/*
				Print out information on runstats
			*/
			fprintf(stdout, "\trunStats:\n\t\tnumSecondaryAttempts=%ld\n\t\tnumFindNextSecondary=%ld\n\t\tnumNonEmpty=%ld\n\t\tnumNonInfected=%ld\n\t\tnumSuccessful=%ld\n", runStats.numSecondaryAttempts, runStats.numFindNextSecondary, runStats.numNonEmpty, runStats.numNonInfected, runStats.numSuccessful);
			fclose(fDPC);
		}
	}
	if (i == pParams->numIts)
	{
		/*
			Dump out a file of the last simulation that was run this time
			Can be useful to keep track of which results are which if repeatedly do simulations with same file name
		*/
		sprintf(outFile, "%s%clastRunNumber.txt", pParams->outStub, C_DIR_DELIMITER);
		fOut = fopen(outFile, "wb");
		if (fOut)
		{
			fprintf(fOut, "%d", pParams->numIts);
			fclose(fOut);
		}
	}
	fclose(fEnd);
	return retVal;
}

/*
	Initialise the dispersal kernel
*/
int setupDispersal(t_Dispersal *pDispersal, t_Landscape *pLandscape, double dispScale, double *pRateSecInf)
{
	int		i,retVal,x,y;
	double	checkDisp, thisDistSq;

	fprintf(stdout, "setupDispersal()\n\t");
	checkDisp = 0.0;
	retVal = 0;
	pDispersal->numProbs = pLandscape->numCols * pLandscape->numRows;
	pDispersal->aProbs = malloc(sizeof(double) * pDispersal->numProbs);
	if(pDispersal->aProbs)
	{
		fprintf(stdout, "Functional form type kernel\n\t\t");
		for(i=0;i<pDispersal->numProbs;i++)
		{
			posToGrid(i, pLandscape->numCols, &x, &y);
			thisDistSq = (double) x * (double) x + (double) y * (double) y;
			pDispersal->aProbs[i] = exp(-sqrt(thisDistSq)/dispScale);
			pDispersal->aProbs[i] /= (dispScale*dispScale*2.0*_PI);
			if(x==0)
			{
				pDispersal->aProbs[i] /= 2.0;
			}
			if(y==0)
			{
				pDispersal->aProbs[i] /= 2.0;
			}
			pDispersal->aProbs[i] *= 4.0;
			checkDisp += pDispersal->aProbs[i];
			if(i % _SETUP_DISPERSAL_PRINT_DOT == 0)
			{
				fprintf(stdout, ".");
			}
		}
		pDispersal->inCell = pDispersal->aProbs[0];
		pDispersal->onLandscape = checkDisp;
		retVal = 1;
		/*
			Store the dispersal kernel as a list of cumulative sums
			(to make it easy to find what infects what)
		*/
		for(i=1;i<pDispersal->numProbs;i++)
		{
			pDispersal->aProbs[i] += pDispersal->aProbs[i-1];
		}
		fprintf(stdout, "\n\t\tinCell=%f onLandscape=%f\n", pDispersal->inCell, pDispersal->onLandscape);
		/*
			Note that it is possible to have onLandscape > 1.0 for certain dispersal kernels (particularly long ranged ones)
			This can happen even thought the kernel is normalised because of discretisation error in
				adding up the kernel but approximating the value of each cell by the value at its centre
			If this happens then renormalise so that all cells can be hit
				(note onLandscape < 1.0 is ok; corresponds to appreciable fraction of dispersals going so far that do not stay within landscape)
		*/
		if (pDispersal->onLandscape > 1.0)
		{
			double dOldRateSec;

			dOldRateSec = *pRateSecInf;
			*pRateSecInf = dOldRateSec * pDispersal->onLandscape;
			for (i = 0; i<pDispersal->numProbs; i++)
			{
				pDispersal->aProbs[i] /= pDispersal->onLandscape;
			}
			pDispersal->onLandscape = 1.0;
			pDispersal->inCell = pDispersal->aProbs[0];
			fprintf(stdout, "\tonLandscape > 1.0 -> Renormalising kernel (and rateSecInf)...\n");
			fprintf(stdout, "\t\tinCell=%f onLandscape=%f\n", pDispersal->inCell, pDispersal->onLandscape);
			fprintf(stdout, "\t\trateSecInf=%f (was %f)\n", *pRateSecInf, dOldRateSec);
		}
	}
	return retVal;
}

int main(int argc, char **argv)
{
	t_Params		sParams;
	t_Landscape		sLandscape;
	t_PriInf		sPriInf;
	t_Epidemic		sEpidemic;
	t_Dispersal		sDispersal;
	clock_t			beforeClock;
	clock_t			afterClock;
	double			totalSeconds;

	seedRandom();
	if(readParams(&sParams, argc, argv))
	{
		if(readLandscape(&sLandscape,sParams.filePropFull,sParams.fileRelInf,sParams.fileRelPri,sParams.fileRelSus,sParams.cellThresh,sParams.outStub))
		{
			if(setupPrimary(&sPriInf,&sLandscape,sParams.ratePriInf))
			{
				if(setupDispersal(&sDispersal, &sLandscape, sParams.dispScale, &sParams.rateSecInf))
				{
					if(setupEpidemic(&sEpidemic,&sLandscape,sParams.dispScale))
					{
						beforeClock = clock();
						runEpidemics(&sParams,&sLandscape,&sPriInf,&sDispersal,&sEpidemic);
						afterClock = clock();
						totalSeconds = ((double)afterClock-(double)beforeClock)/((double)CLOCKS_PER_SEC);
						fprintf(stdout, "%d iterations in %.3f seconds\n", sParams.numIts, totalSeconds);
					}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}