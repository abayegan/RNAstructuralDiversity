#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include "misc.h"
using namespace std;


double ** BasePairProbabilities(char * , int);
double energyOfStructure(char * , char *);
double structureProbability(double);
double expectedPositionalEntropy(double ** , int);
double morganHiggsStructuralDiversity (double ** , int);
double viennaStructuralDiversity (double ** , int);
double expectedBpDistance(double ** , int , char *);
double ensembleDefect(double ** , int , char *);
double expectedNumborOfBasepairs(double ** , int);
double expectedPropOfTargetContacts(double ** , int , char *);
double ensembleHammingDistance(double **, int, char *);
double ensembleBasePairDistance(double ** , int , char *);
double positionalEntropyDistance(double ** , int, char *);
void printPositionalEntropy(double ** , int , char *, int );
int * getBasePairs(char * , int);
void usage(char *);

extern double mfe;
extern char * mfeStr;
int dangleFlag = 2, t99Flag = 0,cutPoint=-1;
double T= 37.;
double kT = (T + 273.15)*1.98717/1000.; /* kT in kcal/mol */


int main(int argc , char *argv[])
{
	char * targetStr, * targetSeq, * seq;
	int targetStrFlag=0 , targetSeqFlag=0 , posEntFlag = 0, n ;
	double strEng;
	int verboseFlag = 0 , normFlag = 0;
	if (argc < 2)
	{
		usage(argv[0]);
	}
	/*print help*/
	if (!strcmp(argv[1],"-h"))
	{
		printf("\ncompute various structural diversity measures:\n\n"
		"\t-s sequence: The input RNA sequence\n\n"
		"\t-d 0|2: set dangling end contributions to 0 or 2 (default is 2)\n\n"
		"\t-c targetStructure: target structure for computation of expected base pair distance,ensemble defect,"
		"expected number of native contacts and expected ensemble distance. if it is not defined the"
		" minimum free energy structure is used to calculate these measures\n\n"
		"\t-q targetSequence: target sequence for computation of ensemble diversity between two sequences\n\n"
		"\t-t temperature: fold sequences at the given temperature\n\n"
		"\t-u 99|04 : fold sequnces using Turner99 or Turner2004 energy parameters" 
		"\t-n : normalize measures by length\n\n"
		"\t-v : verbose mode\n\n"
		"\t-h : print help\n\n");
		usage(argv[0]);
	}
	
	/*parse input command */
	//targetStr = (char* ) malloc(sizeof(char)*(strlen(seq)+1));
	if(argv[1][0]!= '-')
		usage(argv[0]);
	if (argc>1)
	{
		for (int i=1;i<argc;i++)
		{
			if (argv[i][0] == '-')
			{
				switch (argv[i][1])
				{
					case 's':
						if(argc>i+1 && argv[i+1][0]!='-')
						{
							seq = argv[i+1];
							n = strlen(seq);
						}
						else
						{
							printf("\nerror in the input sequence");
							usage(argv[0]);
						}
						break;
					case 'd':
						if(argc>i+1 && !strcmp(argv[i+1],"0"))
							dangleFlag=0;
						else if (argc>i+1 && !strcmp(argv[i+1],"2"))
							dangleFlag=2;
						else
						{
							printf("\nerror in the input dangle flag");
							usage(argv[0]);
						}
						break;
					case 'q' : 
						if(argc>i+1 && argv[i+1][0]!='-' && strlen(argv[i+1]) == n)
						{
							targetSeq = argv[i+1];
							targetSeqFlag = 1;
						}
						else
						{
							printf("\nerror in the input target sequence");
							usage(argv[0]);
						}
						break;
					case 'c' : 
						if(argc>i+1 && argv[i+1][0]!='-')
						{
							targetStr = argv[i+1];
							targetStrFlag = 1;
						}
						else
						{
							printf("\nerror in the input target structure");
							usage(argv[0]);
						}
						break;
					case 'v' :
						verboseFlag = 1;
						break;
					case 'n' : 
						normFlag = 1;
						break;
					case 'p' :
						posEntFlag = 1;
						break;
					case 't' :
						if(argc>i+1 && argv[i+1][0]!='-')
						{
							T = atof(argv[i+1]);
							kT=(T + 273.15)*1.98717/1000.;
						}
						else
						{
							printf("\nerror in the input temperature");
							usage(argv[0]);
						}
						break;
					case 'u':
						if(argc>i+1 && !strcmp(argv[i+1],"99"))
							t99Flag = 1;
						else if (argc>i+1 && !strcmp(argv[i+1],"04"))
							t99Flag = 0;
						else
						{
							printf("\nerror in the input Turner flag");
							usage(argv[0]);
						}
						break;
					default : 
						usage(argv[0]);
				}
			}
		}
	}
	size_t found = string(seq).find('&');
	if(found!=string::npos){
		cout<<"point0000"<<endl;
		cutPoint=found;
		for(int i=found+1;i<strlen(seq);i++)
			seq[i-1] = seq[i];
	}
	CheckSequence(seq);
	printf("seq:%s\n",seq);
	double ** bpr = BasePairProbabilities(seq,n);
	/*target structure is given */
	if(targetStrFlag)
	{
		strEng = energyOfStructure(seq,targetStr); 
	}
	/*otherwise use minimum free energy structure as the target*/
	else
	{
		targetStr = mfeStr;
		strEng = mfe;
	}
	string tstr = targetStr;
	int tLength = 2 * count(tstr.begin(), tstr.end(), '(');
	/*concise print*/
	if(verboseFlag == 0)
	{
		if(normFlag == 0)
		{
			printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",n,strEng,structureProbability(strEng), expectedPositionalEntropy(bpr,n),
			morganHiggsStructuralDiversity(bpr,n),viennaStructuralDiversity(bpr,n),expectedBpDistance(bpr,n,targetStr), 
			ensembleDefect(bpr,n,targetStr), expectedNumborOfBasepairs(bpr,n),expectedPropOfTargetContacts(bpr,n,targetStr)/2);
			if(targetSeqFlag)
				printf("\t%lf\t%lf\t%lf\n", ensembleHammingDistance(bpr,n,targetSeq),ensembleBasePairDistance(bpr,n,targetSeq),positionalEntropyDistance(bpr,n,targetSeq));
			else
				printf("\n");
		}
		else
		{
			printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",n,strEng/n,structureProbability(strEng), expectedPositionalEntropy(bpr,n)/n,
			morganHiggsStructuralDiversity(bpr,n)/n,viennaStructuralDiversity(bpr,n)/n,expectedBpDistance(bpr,n,targetStr)/n, 
			ensembleDefect(bpr,n,targetStr)/n, expectedNumborOfBasepairs(bpr,n)/n,expectedPropOfTargetContacts(bpr,n,targetStr)/tLength);
			if(targetSeqFlag)
				printf("\t%lf\t%lf\t%lf\n", ensembleHammingDistance(bpr,n,targetSeq)/n,ensembleBasePairDistance(bpr,n,targetSeq)/n,positionalEntropyDistance(bpr,n,targetSeq)/n);
			else
				printf("\n");
		}
	}
	else if(verboseFlag == 1)
	{
		if(normFlag == 0)
		{
			printf("sequence: %s\tlength: %d\ntarget structure: %s\nenergy of target structure: %lf\nprobability of target structure: %lf\n"
			"expected Positional entropy: %lf\nMorgan-Higgs structureal diversity: %lf\nVienna structural diversity: %lf\n"
			"expected base pair distance: %lf\nensemble defect: %lf\nexpected number of base pairs: %lf\nexpected proportion of native contacts: %lf\n",
			seq, n, targetStr,strEng,structureProbability(strEng), expectedPositionalEntropy(bpr,n),
			morganHiggsStructuralDiversity(bpr,n),viennaStructuralDiversity(bpr,n),expectedBpDistance(bpr,n,targetStr), 
			ensembleDefect(bpr,n,targetStr), expectedNumborOfBasepairs(bpr,n),expectedPropOfTargetContacts(bpr,n,targetStr)/2);
			if(targetSeqFlag)
				printf("ensemble hamming distance: %lf\nensemble base pair distance: %lf\nexpected positional entropy distance: %lf\n", ensembleHammingDistance(bpr,n,targetSeq),
				ensembleBasePairDistance(bpr,n,targetSeq),positionalEntropyDistance(bpr,n,targetSeq));
		}
		else
		{
			printf("sequence: %s\tlength: %d\ntarget structure: %s\nnormalized energy of target structure: %lf\nprobability of target structure: %lf\n"
			"normalized expected Positional entropy: %lf\nnormalized Morgan-Higgs structureal diversity: %lf\nnormalized Vienna structural diversity: %lf\n"
			"normalized expected base pair distance: %lf\nnormalized ensemble defect: %lf\nnormalized expected number of base pairs: %lf\nnormalized expected proportion of native contacts: %lf\n",
			seq, n, targetStr,strEng/n,structureProbability(strEng), expectedPositionalEntropy(bpr,n)/n,
			morganHiggsStructuralDiversity(bpr,n)/n,viennaStructuralDiversity(bpr,n)/n,expectedBpDistance(bpr,n,targetStr)/n, 
			ensembleDefect(bpr,n,targetStr)/n, expectedNumborOfBasepairs(bpr,n)/n,expectedPropOfTargetContacts(bpr,n,targetStr)/tLength);
			if(targetSeqFlag)
				printf("ensemble hamming distance: %lf\nensemble base pair distance: %lf\nexpected positional entropy distance: %lf\n", ensembleHammingDistance(bpr,n,targetSeq)/n,
				ensembleBasePairDistance(bpr,n,targetSeq)/n,positionalEntropyDistance(bpr,n,targetSeq)/n);
		}
		
	}
	if (posEntFlag==1)
		printPositionalEntropy(bpr,n,targetSeq,targetSeqFlag);
	cout <<"last point"<<endl;
	return 0;
}

void usage(char * name)
{
	printf("\nUSAGE: %s [-s sequence] [-d 0|2] [-c tragetStructure] [-t temperature] [-q tragetSequence] [-u 99|04] [-n] [-v]\ntry %s -h for help\n\n",name, name);
	exit(1);
}

