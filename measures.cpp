#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <string.h>
#include "misc.h"
extern "C"{
#include "utils.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"
#include "RNAstruct.h"
#include "read_epars.h"
#include "cofold.h"
#include "part_func_co.h"
//#include "data_structures.h"
}

using namespace std;

extern double kT;
double ensEng; /* the Gibbs free energy of ensemble */
double mfe;
char * mfeStr;
extern int dangleFlag,cutPoint;
extern double T;
	
int * getBasePairs(char * str, int n)
{
	vector <int> stack;
	int * bp = new int[n];
	for(int i=0; i<n;i++)
	{
		if (str[i] == '(' )
			stack.push_back(i);
		else if (str[i] == ')')
		{
			int m = stack.back();
			bp[i] = m;
			bp[m] = i;
			stack.pop_back();
		}
		else
			bp[i] = -1;
	}
		
	return bp;
}

double ** BasePairProbabilities(char * seq , int n)
{
	char *struct1; /* the base pair brobability information in a psuedo dot bracket notation */
	FLT_OR_DBL *bppm;
	mfeStr = (char* )space(sizeof(char)*(n+1));
	dangles = dangleFlag;
	temperature = T;
	//cout<<"point1"<<endl;
	if(cutPoint==-1){
		//cout<<"point2"<<endl;
		//cout<<"temp "<<temperature<<endl;
		//cout<<"d "<<dangles<<endl;
		mfe = fold(seq,mfeStr);
		//printf("seq:%s,temp:%f,d:%d,mfe:%f\n",seq,T,dangleFlag,mfe);
		//cout<<"point5"<<mfe<<endl;
	 	pf_scale = exp((-mfe/kT)/n);
		update_pf_params(n);
		ensEng = pf_fold(seq, struct1);
		bppm = export_bppm();
		//cout<<"point6"<<endl;
		//~ for(int i=0;i<n;i++){
			//~ for(int j=i;j<n;j++){
				//~ printf	("%.5f\t",bppm[iindx[i+1] - (j+1)  ]);
				//~ cout<<endl;
			//~ }
		//~ }
	}
	else{
		//cout<<"cutPoint"<<cutPoint<<endl;
		cut_point = cutPoint;
	    //initialize_cofold(string(seq).size());
	    mfe = cofold(seq,mfeStr);
	   // printf("seq:%s,temp:%f,d:%d,mfe:%f",seq,T,dangleFlag,mfe);
	    pf_scale = exp((-mfe/kT)/strlen(seq));
	    update_co_pf_params(strlen(seq));
	    cofoldF cof = co_pf_fold(seq,struct1);
	    printf("cof.F0AB:%d,cofFAB%d cofFcAB%d cofFA:%d\n",cof.F0AB,cof.FAB,cof.FcAB,cof.FA);
	    ensEng = cof.FAB;
	    bppm = export_co_bppm();
	    //~ for(int i=0;i<strlen(seq);i++){
			//~ for(int j=i;j<strlen(seq);j++){
				//~ printf	("%.5f\t",bppm[iindx[i+1] - (j+1)  ]);
				//~ cout<<endl;
			//~ }
		//~ }
	}
	//cout<<"point 10 "<<endl;	
	double **bppr = new double *[n];
	for(int i = 0; i < n; i++)
		bppr[i] = new double[n+1];
	//cout<<"point 11 "<<endl;	
	for(int i=0;i<n;i++)
		for(int j=0;j<n+1;j++)
			bppr[i][j]=0;		
	//cout<<"point 12 "<<endl;	
	double pr = 0;
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			bppr[i][j]=bppm[iindx[i+1] - (j+1)  ];
			bppr[j][i]=bppm[iindx[i+1] - (j+1)  ];
		}
	}
	//cout<<"point 13		 "<<endl;	
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			pr +=bppr[i][j];
		}
		bppr[i][n]= 1-pr;
		pr=0;
	}
	//print base pair probabilities
	/*cout<< "\nbase pair probabilities:" << endl;
	for (int i=0; i<n;i++)
	{
		for (int j=0;j<n+1;j++)
			printf	("%.5f\t",bppr[i][j]);
		cout<<endl;
	}*/
	free_arrays();
    free_pf_arrays();
    free_co_arrays();
    free_co_pf_arrays();
	return bppr;
}


double energyOfStructure(char * seq , char * str)
{
	return energy_of_structure(seq , str , 0);
}

double structureProbability(double eng)
{
	double boltPr = exp(-eng/kT) / exp(-ensEng/kT);
	return boltPr;
}

double expectedPositionalEntropy(double ** bppr, int n)
{
	double h=0;
		for (int i=0;i<n;i++)
		for (int j=0;j<n+1;j++)
		{
			if(bppr[i][j]>0)
				h += bppr[i][j]*log(bppr[i][j]);
		}
	return -h;
}

double morganHiggsStructuralDiversity (double **bppr, int n)
{
	double sum=0;
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
			sum += pow(bppr[i][j],2);
		sum += pow(bppr[i][n],2);
	}
	return n-sum; 	
}

double viennaStructuralDiversity (double **bppr, int n)
{
	double sum = 0;
	for (int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			sum += bppr[i][j]*(1-bppr[i][j]);
	return sum;
}

double expectedBpDistance(double **bppr, int n , char * s)
{
	double sum = 0;
	int * bpList = new int [n];
	bpList = getBasePairs(s,n);
	for (int i=0;i<n;i++)
	{
		for(int j=i+1;j<n;j++)
		{
			if(bpList[i]!=j )
			  sum += bppr[i][j];
			else 
			  sum += 1 - bppr[i][j];
		}
	}

	return sum;
}

double ensembleDefect(double **bppr, int n , char * s)
{
	double sum = 0;
	int * bpList = new int [n];
	bpList = getBasePairs(s,n);
	
	for(int i=0;i<n;i++)
	{
		if (bpList[i] == -1)
			sum += bppr[i][n];
		else
			sum += bppr[i][bpList[i]];
	}
	return n-sum;
}

double expectedNumborOfBasepairs(double **bppr, int n)
{
	double sum = 0;
	for (int i=0; i<n-3; i++)
		for(int j=i+4;j<n;j++)
			sum += bppr[i][j];
	return sum;
}

double expectedPropOfTargetContacts(double ** bppr , int n , char * s)
{
	double sum = 0;
	int * bpList = new int [n];
	bpList = getBasePairs(s,n);
	
	for(int i=0;i<n;i++)
	{
		if(bpList[i]!=-1)
			sum += bppr[i][bpList[i]];
	}
	return sum;
}

double ensembleHammingDistance(double ** bppr , int n , char * targetSeq)
{
	double sum = 0;
	double ** bppr_b = BasePairProbabilities(targetSeq,strlen(targetSeq));
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
			sum += pow(bppr[i][j]-bppr_b[i][j],2);
		sum += pow(bppr[i][n]-bppr_b[i][n],2);
	}
	return sqrt(0.5*sum);
}

double ensembleBasePairDistance(double ** bppr , int n , char * targetSeq)
{
	double sum = 0;
	double ** bppr_b = BasePairProbabilities(targetSeq,strlen(targetSeq));
	for(int i=0;i<n;i++)
		for (int j=i+1;j<n;j++)
			sum += pow(bppr[i][j]-bppr_b[i][j],2);
	return sqrt(sum);
}

double positionalEntropyDistance(double ** bppr , int n , char * targetSeq)
{
	double ** bppr_b = BasePairProbabilities(targetSeq,strlen(targetSeq));
	double h_ab = 0;
	double * h_a = new double [n];
	double * h_b = new double [n];
	for(int i=0;i<n;i++)
	{
		h_a[i] = 0;
		h_b[i] = 0;
	}
	for (int i=0;i<n;i++)
		for (int j=0;j<n+1;j++)
		{
			if(bppr[i][j]>0)
				h_a[i] += -bppr[i][j]*log(bppr[i][j]);
			if(bppr_b[i][j]>0)
				h_b[i] += -bppr_b[i][j]*log(bppr_b[i][j]);
		}
	for (int i=0;i<n;i++)
		h_ab += pow(h_a[i] - h_b[i],2);
	return sqrt(h_ab);
}

void printPositionalEntropy(double ** bppr , int n , char * targetSeq, int tsFlag)
{
	double * h_a = new double [n];
	double * h_b = new double [n];
	for(int i=0;i<n;i++)
		h_a[i] = 0;
	for (int i=0;i<n;i++)
		for (int j=0;j<n+1;j++)
		{
			if(bppr[i][j]>0)
				h_a[i] += -bppr[i][j]*log(bppr[i][j]);
		}
	if (tsFlag==1)
	{
		double ** bppr_b = BasePairProbabilities(targetSeq,strlen(targetSeq));
		for(int i=0;i<n;i++)
			h_b[i] = 0;
		for (int i=0;i<n;i++)
			for (int j=0;j<n+1;j++)
			{
				if(bppr_b[i][j]>0)
					h_b[i] += -bppr_b[i][j]*log(bppr_b[i][j]);
			}
	}
	printf("position\tpositionalEnt1\tpositionalEnt2\tabsolute difference\n");
	for(int i=0;i<n;i++)
	{
		if(tsFlag==0)
			printf("%d\t%lf\n",i,h_a[i]);
		else if(tsFlag==1)
			printf("%d\t%lf\t%lf\t%lf\n",i,h_a[i],h_b[i],fabs(h_a[i]-h_b[i]));
	}
}
	
