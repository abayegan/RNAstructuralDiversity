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
#include "cofold.h"
#include "part_func_co.h"
}

using namespace std;

extern double kT;
extern int dangleFlag;
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

double ** BasePairProbabilities(char * seq , int n, double * mfe, char* mfeStr, double *ensEng , int cutPoint)
{
	FLT_OR_DBL *bppm;
	dangles = dangleFlag;
	temperature = T;
	if(cutPoint==-1){
		*mfe = fold(seq,mfeStr);
		//printf("seq:%s,temp:%f,d:%d,mfe:%f\n",seq,T,dangleFlag,*mfe);
	 	pf_scale = exp((-*mfe/kT)/n);
		//update_pf_params(n);
		*ensEng = pf_fold(seq, NULL);
		bppm = export_bppm();

	}
	else{
		cut_point = cutPoint;
	    *mfe = cofold(seq,mfeStr);
	   // printf("seq:%s,temp:%f,d:%d,mfe:%f",seq,T,dangleFlag,mfe);
	    pf_scale = exp((-*mfe/kT)/n);
	    cofoldF cof = co_pf_fold(seq,NULL);
	   // printf("cof.F0AB:%d,cofFAB:%d cofFcAB:%d cofFA:%d\n",cof.F0AB,cof.FAB,cof.FcAB,cof.FA);
	    *ensEng = cof.FAB;
	    bppm = export_co_bppm();

	}	
	double **bppr = new double *[n];
	for(int i = 0; i < n; i++)
		bppr[i] = new double[n+1];
	for(int i=0;i<n;i++)
		for(int j=0;j<n+1;j++)
			bppr[i][j]=0;		
	double pr = 0;
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			bppr[i][j]=bppm[iindx[i+1] - (j+1)  ];
			bppr[j][i]=bppm[iindx[i+1] - (j+1)  ];
		}
	}
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			pr +=bppr[i][j];
		}
		bppr[i][n]= 1-pr;
		pr=0;
	}
	//print base pair probabilities
	//cout<< "\nbase pair probabilities:" << endl;
	//for (int i=0; i<n;i++)
	//{
		//for (int j=0;j<n+1;j++)
			//printf	("%.5f\t",bppr[i][j]);
		//cout<<endl;
	//}
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

double structureProbability(double eng, double ensEng)
{
	double boltPr = exp(-eng/kT) / exp(-ensEng/kT);
	return boltPr;
}

double expectedPositionalEntropy(double ** bppr, int n, int type)
{	
	//-------type=1|2|3: regular|binary|ternary positional entropy
	
	double h=0;
	int i,j;
	if(type==1){
		for (i=0;i<n;i++)
			for (j=0;j<n+1;j++)
			{
				if(bppr[i][j]>0)
					h += bppr[i][j]*log(bppr[i][j]);
			}
	}
	else if(type==2){
		for (i=0;i<n;i++){
			if(bppr[i][n]>0)
				h+= bppr[i][n]*log(bppr[i][n]);
			if(1-bppr[i][n]>0)
				h+= (1-bppr[i][n])* log(1-bppr[i][n]);
		}
	}
	else{
		for (i=0;i<n;i++){
			double rsum=0,lsum=0;
			for (j=i+1;j<n;j++)
					rsum+= bppr[i][j];
			for (j=i-1;j>0;j--)
					lsum+= bppr[i][j];
			if(rsum!=0)
				h+= rsum*log(rsum);
			if(lsum!=0)
				h+= lsum*log(lsum);
			if(bppr[i][n]!=0)
				h+= bppr[i][n]*log(bppr[i][n]);
		}
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
	int bpcnt = 0;
	for(int i=0;i<n;i++)
	{
		if(bpList[i]!=-1)
			sum += bppr[i][bpList[i]];
	}
	
	return sum/2;
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
	//-------alternative approach-------------
	//double sum1 = 0;	
	//int i,j;
	//for(int i=0;i<n;i++)
	//{
		//if (bpList[i] == -1)
			//sum1 += 1 - bppr[i][n]; //probability of i being paired
		//else
			//for (j=0;j<n;j++)
				//if (bpList[i] == j)
					//sum1 += 1 - bppr[i][j]; // probability that (i,j) do not exist
	//}
	//cout << "ensDefAlt:" << sum1 << endl;
	////return sum;

}

double epsilonNearExpectedBpDistance0(double **bppr, int n , char * s, int epsilon)
{
	double sum=0;
	int * bpList = new int [n];
	bpList = getBasePairs(s,n);
	int found;
	for (int i=0;i<n;i++)
	{
		for(int j=i+1;j<n;j++)
		{
			if (bpList[i]==j){
				sum+= 1+bppr[i][j];
				for(int y=j-epsilon;y<=j+epsilon;y++)
					if(i<=y  && 0<=y && y<n)
						sum -= bppr[i][y];
				for(int x=i-epsilon;x<=i+epsilon;x++)
					if(x<=j  && 0<=x && x<n)
						sum -= bppr[x][j];
			}
			else if(!((bpList[i]!=-1 && bpList[i]<=j+epsilon && bpList[i]>=j-epsilon) || (bpList[j]!=-1 && bpList[j]<=i+epsilon && bpList[j]>=i-epsilon))){
				sum += bppr[i][j];
			}
		}
	}
	return sum;
}

double epsilonNearExpectedBpDistance(double **bppr, int n , char * s, int epsilon)
{
	double sum=0;
	int * bpList = new int [n];
	bpList = getBasePairs(s,n);
	int found=1;
	for (int i=0;i<n-3;i++)
	{
		for(int j=i+4;j<n;j++)
		{
			double sum1=0,sum2=0;
			if (bpList[i]==j){
				sum += 2 - bppr[i][n] - bppr[j][n] + bppr[i][j] + bppr[i][n]*bppr[j][n];
				for(int y=j-epsilon;y<=j+epsilon;y++)
					if(i<=y && y>=0 && y<n)
						sum -= bppr[i][y];
				for(int x=i-epsilon;x<=i+epsilon;x++)
					if(x<=j && x>=0 && x<n)
						sum -= bppr[x][j];
			}
			else if(!((bpList[i]!=-1 && bpList[i]<=j+epsilon && bpList[i]>=j-epsilon) || (bpList[j]!=-1 && bpList[j]<=i+epsilon && bpList[j]>=i-epsilon))){
				//printf("+bppr[%d][%d]\n",i,j);
				sum += bppr[i][j];
			}
		}
	}
	
	return sum;
}

double epsilonNearEnsembleDefect0(double **bppr, int n , char * s,int epsilon)
{
	double sum=0;
	int k;
	int * bpList = new int [n];
	
	//---------------Peter's Definition-------------------------
	bpList = getBasePairs(s,n);
	for(int i=0;i<n;i++)
		if (bpList[i]==-1){
			bpList[i] = i;
		}
	for(int i=0;i<n;i++)
		bppr[i][i] =bppr[i][n];
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			if(fabs(j-bpList[i])>epsilon){
				sum+=bppr[i][j];
			}
		}
	}
	
	//---------------mine-----------------------------------------
	bpList = getBasePairs(s,n);
	double sum1 = 0,sum2 = 0, sum3=0;	
	int i,j;
	for(int i=0;i<n;i++)
	{
		if (bpList[i] == -1){
			sum1 += 1 - bppr[i][n];//probability of i being paired
			//printf("%i,%f + ",i,1-bppr[i][n]);
		}
		else
		{
			sum2 = 0;
			for (j=bpList[i]-epsilon;j<=bpList[i]+epsilon;j++){
				if(0<=j && j<n){
					sum2 += bppr[i][j]; 
				}
			}
			sum3 += 1-sum2; //probability of i being paired in [j-epsilon,j+epsilon]
			//printf("%d,%d,%f + ",i,j-1,1-sum2);
		}			
	}
	//cout <<"Peter:"<<sum<<endl;
	//cout << "mine:"<<sum1+sum3<<endl;
	return sum;
}

double epsilonNearEnsembleDefect(double **bppr, int n , char * s,int epsilon)
{
	int * bpList = new int [n];
	bpList = getBasePairs(s,n);
	double sum1 = 0,sum2 = 0, sum3=0;	
	int i,j;
	for(int i=0;i<n;i++)
	{
		if (bpList[i] == -1)
			sum1 += 1 - bppr[i][n];//probability of i being paired
		else
		{
			sum2 = 0;
			for (j=bpList[i]-epsilon;j<=bpList[i]+epsilon;j++){
				if(0<=j && j<n){
					sum2 += bppr[i][j]; 
				}
			}
			sum3 += 1-sum2; //probability of i being paired in [j-epsilon,j+epsilon]
		}			
	}
	return sum1+sum3;
}

double ensembleHammingDistance(double ** bppr , double ** bppr_b ,int n)
{
	double sum = 0;
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
			sum += pow(bppr[i][j]-bppr_b[i][j],2);
		sum += pow(bppr[i][n]-bppr_b[i][n],2);
	}
	return sqrt(0.5*sum);
}

double ensembleBasePairDistance(double ** bppr , double ** bppr_b  , int n)
{
	double sum = 0;
	for(int i=0;i<n;i++)
		for (int j=i+1;j<n;j++){
			//cout <<i<<','<<j<<'\t' << bppr[i][j] << '\t' << bppr_b[i][j] <<endl; 
			sum += pow(bppr[i][j]-bppr_b[i][j],2);
		}
	return sqrt(sum);
}

double positionalEntropyDistance(double ** bppr , double ** bppr_b  , int n)
{
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
		h_ab += h_a[i]/n - h_b[i]/n;
	return h_ab;
}

void printPositionalEntropy(double ** bppr ,double ** bppr_b  , int n, int tsFlag)
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
	if(tsFlag){
		
	for(int i=0;i<n;i++)
		h_b[i] = 0;
	for (int i=0;i<n;i++)
		for (int j=0;j<n+1;j++)
		{
			if(bppr_b[i][j]>0)
				h_b[i] += -bppr_b[i][j]*log(bppr_b[i][j]);
		}
	}
	if(tsFlag==0)
		printf("position\tpositionalEnt\n");	
	else
		printf("position\tpositionalEnt1\tpositionalEnt2\tdifference(Ent1-Ent2)\n");	
	for(int i=0;i<n;i++)
	{
		if(tsFlag==0)
			printf("%d\t%lf\n",i,h_a[i]);
		else if(tsFlag==1)
			printf("%d\t%lf\t%lf\t%lf\n",i,h_a[i]/n,h_b[i]/n,h_a[i]/n-h_b[i]/n);
	}
}

double ensembleVariationDistance(double ** bppr , double ** bppr_b  , int n)
{
	double sum=0;
	for (int i=0;i<n;i++)
		for(int j=i+1;j<n;j++)
			sum += fabs(bppr[i][j] - bppr_b[i][j]);
	return sum/2;
}

double RodrigoRobustness(double ** bppr , int n , char * s)
{
	string ss = string(s);
	float sum=0.0;
	char * mm = new char[n+1];
	//----compute the average expected base pair distance from all 1-mutatns to the given seq
	for (int i=0; i<n;i++)
	{
		string nuc = "ACGU";
		for (int j=0;j<4;j++)
		{
			if (nuc[j]==s[i]){
				nuc.erase(nuc.begin()+j);
				continue;
			}				
		}
		for (int k=0;k<3;k++){
			//-----generate the mutant---
			string mut=ss;
			mut.replace(i,1,1,nuc[k]);
			strcpy(mm,mut.c_str());
			//-----fold the mutant------
			char * tMFEStr;
			int cutp=-1;
			double ensEngTarget,mfeTarget;
			tMFEStr= (char* )space(sizeof(char)*(n+1));
			size_t f = mut.find('&');
			if(f!=string::npos){
				cutp=f+1;
				for(int i=f+1;i<n;i++)
					mm[i-1] = mm[i];
				mm[n-1]='\0';
			}	
			double ** bppr2 = BasePairProbabilities(mm,n,&mfeTarget,tMFEStr,&ensEngTarget,cutp);
			sum += ensembleBasePairDistance(bppr,bppr2,n);
		}
	}
	double avgbpdist = sum/(3*n);
	//-----compute robustness
	double rob = 1 - avgbpdist/(n/2);
	//printf("robustness:%f\n",rob);
	delete [] mm;
	return rob;
}

double expected2DcontactOrder(double ** bppr , int n , char * s)
{
	double sum=0;
	for (int i=0;i<n;i++)
		for (int j=i+1;j<n;j++)
			sum += bppr[i][j]*(j-i);
	return sum;
			
}

void expectedHeight(double **bppr , int n)
{
	double * h = new double [n];
	for(int i=0;i<n;i++)
		h[i] = 0;
	//--------base case--------
	for (int i=1;i<n;i++)
		h[0] += bppr[0][i];
		
	for(int k=0;k<n-1;k++)
	{
		//printf("K=%d\n",k);
		h[k+1] += h[k];
		for (int y=k+2;y<n;y++){
			h[k+1] += bppr[k+1][y];
			//printf("+: %d,%d,%lf\n",k+1,y,bppr[k+1][y]);
		}
		for (int x=0; x<=k ; x++){
			h[k+1] -= bppr[x][k+1];
			//printf("-: %d,%d,%lf\n",x,k+1,bppr[x][k+1]);
		}
	}
	printf("\n****************expected height at each position***********************\n");
	for(int i=0;i<n;i++)
			printf("%d\t%lf\n",i,h[i]);
}

void expectedModifiedHeight(double ** bppr, int n)
{
	double * h = new double [n];
	for(int i=0;i<n;i++)
		h[i] = 0;	
		
	//--------base case--------
	for (int i=1;i<n;i++)
		h[0] += bppr[0][i]*(1/float(i));
		
	for(int k=0;k<n-1;k++)
	{
		//printf("K=%d\n",k);
		h[k+1] += h[k];
		for (int y=k+2;y<n;y++){
			h[k+1] += bppr[k+1][y]*(1/float((y-k-1)));
			//printf("+: %d,%d,%lf\n",k+1,y,bppr[k+1][y]);
		}
		for (int x=0; x<=k ; x++){
			h[k+1] -= bppr[x][k+1] * (1/float((k+1-x)));
			//printf("-: %d,%d,%lf\n",x,k+1,bppr[x][k+1]);
		}
	}
	printf("\n****************expected modified height at each position***********************\n");
	for(int i=0;i<n;i++)
			printf("%d\t%lf\n",i,h[i]);
	
}

void printTriProbabilities(double **bppr , int n){
	int i,j;
	double * lp = new double [n];
	double * rp = new double [n];
	double * up = new double [n];
	for(i=0;i<n;i++)
		lp[i]=rp[i]=up[i]=0.0;
		
	for(i=0;i<n;i++){
		for(j=i+1;j<n;j++){
			rp[i]+=bppr[i][j];
		}
		for(j=0;j<i;j++){
			lp[i]+=bppr[i][j];
		}
		up[i] = 1.0-lp[i]-rp[i];
	}
	printf("\n******Probability of being paired to the right,left,and being unpaired i.e \"().\" respectively******\n");
	for(int i=0;i<n;i++)
			printf("%d\t%lf\t%lf\t%lf\n",i,rp[i],lp[i],up[i]);
}
