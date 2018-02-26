#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <vector>
#include <sstream>
extern "C"{
#include "read_epars.h"
#include "utils.h"   //for space()
}
#include "misc.h"
#include "measures.h"

using namespace std;

void usage(char *);


int dangleFlag = 2, engFlag = 2004,cutPoint=-1,epsilon=2, cutPoint2=-1;
double T= 37.;
double kT = (T + 273.15)*1.98717/1000.; /* kT in kcal/mol */
int MAXSEQNUM=10000; //maximum number of seqs in a fasta file
const int FIELDSNUM=17;


int main(int argc , char *argv[])
{	
	char * targetStr, * targetSeq, * seq0, * seq, * mfeStr, * targetMFEStr2;
	int targetStrFlag=0 , targetSeqFlag=0 , posEntFlag = 0,seqFlag=0, outFlag=0, verboseFlag=0 , normFlag=0, Hflag=0 , MHflag=0, PRflag=0 , n , n2 ;
	double mfe,strEng,targetMFE2, ensEng, ensEng2;
	string fastaName="";
	string out;	
	vector <string> seqs;
	vector<string> strs;
	double results[FIELDSNUM];
	for (int i=0;i<FIELDSNUM;i++)
		results[i]=-1;
	double ** bpr2;
	
	if (argc < 2)
	{
		usage(argv[0]);
	}
	/*print help*/
	if (!strcmp(argv[1],"-h"))
	{
		
		printf("\nUSAGE: %s -s sequence -f fasta  [-d 0|2] [-c tragetStructure] [-t temperature] [-q tragetSequence] [-u turner99|turner04|andronescu07] [-x epsilon] [-o output] [-n] [-v]\ntry %s -h for help\n\n",argv[0], argv[0]);
		printf("\ncompute various structural diversity measures:\n\n"
		"\t-s sequence: The input RNA sequence\n\n"
		"\t-f fasta: The input fasta file. Target structure for each rna can be define after its sequence. Otherwise MFE structure will be used. \n\n"
		"\t-d 0|2: set dangling end contributions to 0 or 2 (default is 2)\n\n"
		"\t-c targetStructure: target structure for computation of expected base pair distance,ensemble defect,expected number of native contacts\n"
		"\tand expected ensemble distance. if it is not defined the"
		" minimum free energy structure is used to calculate these measures\n\n"
		"\t-q targetSequence: target sequence for computation of ensemble diversity between two sequences\n\n"
		"\t-t temperature: fold sequences at the given temperature\n\n"
		"\t-u turner99|turner04|andronescu07 : energy parameter used to fold the sequences\n\n" 
		"\t-n : normalize measures by length\n\n"
		"\t-x : epsilon for calculation of epsilon-near expected bp distance and ensemble defect (default is 2)\n\n"
		"\t-v : verbose mode\n\n"
		"\t-h : print help\n\n"
		"\t-o : measures to be computed separated by comma. Each measure is defined by an integer value:\n"
		"\t1: expected positional entropy(if -n is not used computes the sum over all poistions)\n\t2:Morgan Higgs diversity\n"
		"\t3: Vienna diversity\n\t4: expected base pair distance\n"
		"\t5: ensemble defect\n\t6: expected number of base pairs\n"
		"\t7: expected proportion of target contacts\n\t8: expected 2D contact order(if -n is not used computes the sum over all poistions)\n"
		"\t9: epsilon-near BP distance\n\t10: epsilon-near ensemble defect\n"
		"\t11: robustness\n\t12: ensemble Hamming distance\n"
		"\t13: ensemble base pair distance\n\t14: positional entropy distance\n"
		"\t15: expected height\n\t16: expected modified height\n"
		"\t17: For each position the probability of bein basepaired to left,right and being unpaired is reported respectively."
		"\tExample: using -o 1,2,5,7 will only compute the corresppnding measures. \n\n");
		exit(1);
	}
	
	/*------------------parse input command---------------------------*/
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
							seq0 = argv[i+1];
							n = strlen(seq0);
							seqFlag=1;
						}
						else
						{
							printf("\nerror in the input sequence");
							usage(argv[0]);
						}
						break;
					case 'f':
						if(argc>i+1 && argv[i+1][0]!='-')
						{
							fastaName = argv[i+1];
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
						if(argc>i+1 && argv[i+1][0]!='-')
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
						if(argc>i+1 && !strcmp(argv[i+1],"turner99"))
							engFlag =  1999;
						else if (argc>i+1 && !strcmp(argv[i+1],"turner04"))
							engFlag = 2004;
						else if(argc>i+1 && !strcmp(argv[i+1],"andronescu07"))
							engFlag = 2007;
						else
						{
							printf("\nerror in the input Turner flag");
							usage(argv[0]);
						}
						break;
					case 'x':
						if(argc>i+1 && argv[i+1][0]!='-')
							epsilon = atoi(argv[i+1]);
						else
						{
							printf("\nerror in x flag");
							usage(argv[0]);
						}
						break;
					case 'o':
						if(argc>i+1 && argv[i+1][0]!='-'){
							out = argv[i+1];
							outFlag=1;
						}
						else
						{
							printf("\nerror in output flag");
							usage(argv[0]);
						}
						break;
					default:
						usage(argv[0]);
				}
			}
		}
	}
	
	if (seqFlag==0 && fastaName==""){
		cerr << "error: either a sequence or a fasta file should be defined in the input!\n";
		exit(1);
	}
	
	//-------load the correct energy parameter file
	char * execPath = getExecPath(argv[0]);
	//printf("energy file path: %s\n",execPath);	
	char turner99[100] = "/energy_params/rna_turner1999.par";
	char turner04[100] = "/energy_params/rna_turner2004.par";
	char andronescu07[100] = "/energy_params/rna_andronescu2007.par";
	
	if(engFlag==1999)
		strcat(execPath,turner99);
	else if (engFlag==2004)
		strcat(execPath , turner04);
	else if (engFlag==2007)
		strcat(execPath , andronescu07);
	else
		usage(argv[0]);
	read_parameter_file(execPath);
	
	//-------read the sequence and target structures
	if(fastaName!=""){
		parseFasta(fastaName,seqs,strs);
		//for (int i=0;i<seqs.size();i++)
			//printf("%s\n\%s\n",seqs.at(i).c_str(),strs.at(i).c_str());
		if (seqs.size() >= MAXSEQNUM){
			char err[100];
			sprintf(err,"error: too much sequences in the fasta file. Maximum is %d\n", MAXSEQNUM);
			cerr << err;
			exit(1);
		}
	}
	else{
		seqs.push_back(seq0);
		if(targetStrFlag)
			strs.push_back(targetStr);
		else
			strs.push_back("MFE");
		}
	
	//------parse the output flag
	vector<int> outVec;
	if(outFlag==1){
		stringstream ss(out);
	    string item;
	    while (getline(ss, item,',')) {
			int o = atoi(item.c_str());
	        outVec.push_back(o);
	        if(o>11 && o<15 && targetSeqFlag==0){
				cerr << "error: target sequence must be defined!\n";
				exit(1);
			}
		}
	}
	else{ //default output is all measures except the ones that need a target seq
		for (int i=1;i<=9;i++)
			outVec.push_back(i);
	}
	
	
	//for (int i=0;i<outVec.size();i++)
		//cout<< outVec.at(i)<<endl;
	char s1[30],s2[30];
	sprintf(s1,"%d.near.BP.dist",epsilon);
	sprintf(s2,"%d.near.ens.defect",epsilon); 
	string outTitle[FIELDSNUM+1] = {" ","exp.pos.ent","Morgan.Higgs.div","Vienna.div","exp.BP.dist","ens.defect","exp.num.of.BPs","exp.native.contacts","exp.2D.contact.order"
									,s1, s2,"robustness", "ens.Hamming.dist", "ens.BP.dist","pos.ent.dist"};
	if(verboseFlag == 0){
		printf("sequence\tlength\ttarget.str\teng.of.target\tprob.of.target\t");
		for (int i=0;i<outVec.size();i++)
			printf("%s\t",outTitle[outVec.at(i)].c_str());
		printf("\n");
	}
	
	if(targetSeqFlag)
	{
		n2 = strlen(targetSeq);
		targetMFEStr2= (char* )space(sizeof(char)*(n+1));
		size_t found2 = string(targetSeq).find('&');
			if(found2!=string::npos){
				cutPoint2=found2+1;
				for(int i=found2+1;i<n2;i++)
					targetSeq[i-1] = targetSeq[i];
				targetSeq[n2-1]='\0';
			}	
		bpr2 = BasePairProbabilities(targetSeq,n2,&targetMFE2,targetMFEStr2,&ensEng2,cutPoint2);
		//printf("seq:%s,temp:%f,d:%d,mfe:%f,MFESTR:%s\n",targetSeq,T,dangleFlag,targetMFE2,targetMFEStr2);
	}
		
	for (int i=0;i<seqs.size();i++){
		//seq = new char[(seqs.at(i).length())+1];
		seq = (char* )space(sizeof(char)*(seqs.at(i).length()+1));
		//targetStr = new char[(seqs.at(i).length())+1];
		targetStr = (char*)space(sizeof(char)*(seqs.at(i).length()+1));
		strcpy(seq, seqs.at(i).c_str());
		strcpy(targetStr, strs.at(i).c_str());
		//printf("%s\n",seq);
		n = strlen(seq);
		if(targetSeqFlag && n2!=n){
			cerr << "error: input sequences and the target sequence must have the same length!\n";
			exit(1);
		}
		//-------check for hybridization----------
		size_t found = string(seq).find('&');
		if(found!=string::npos){
			cutPoint=found+1;
			for(int i=found+1;i<strlen(seq);i++)
				seq[i-1] = seq[i];
			seq[n-1]='\0';
			n = n-1;
			if (targetStr!="MFE")
				if(targetStr[found]=='&'){
					for(int i=found+1;i<strlen(targetStr);i++)
						targetStr[i-1] = targetStr[i];
					targetStr[n] = '\0';
				}
				else
				{
					printf("given sequence and structure are not compatible!\n");
					usage(argv[0]);
				}
		}	
		
		CheckSequence(seq);
		mfeStr = (char* )space(sizeof(char)*(n+1));
		double ** bpr = BasePairProbabilities(seq,n,&mfe,mfeStr,&ensEng,cutPoint);
		//printf("seq:%s,temp:%f,d:%d,mfe:%f,MFESTR:%s\n",seq,T,dangleFlag,mfe,mfeStr);
		
		if(string(targetStr).find("MFE") == string::npos)
		{
			strEng = energyOfStructure(seq,targetStr); 
		}
		/*otherwise use minimum free energy structure as the target*/
		else
		{
			targetStr = mfeStr;
			strEng = mfe;
		}
		//cout << i <<'\t'<< seq << '\t' << targetStr << endl;
		string tstr = targetStr;
		int tLength = count(tstr.begin(), tstr.end(), '(');
		int num = outVec.size();
		for (int j=0;j<num;j++){
			int idx = outVec.at(j);
			switch (idx){
				case 1:
					results[1] = expectedPositionalEntropy(bpr,n,1);
					break;
				case 2:
					results[2] = morganHiggsStructuralDiversity(bpr,n);
					break;
				case 3:
					results[3] = viennaStructuralDiversity(bpr,n);
					break;
				case 4:
					results[4] = expectedBpDistance(bpr,n,targetStr);
					break;	
				case 5:
					results[5] = ensembleDefect(bpr,n,targetStr);
					break;		
				case 6:
					results[6] = expectedNumborOfBasepairs(bpr,n);
					break;
				case 7:
					results[7] = expectedPropOfTargetContacts(bpr,n,targetStr);
					break;
				case 8:
					results[8] = expected2DcontactOrder(bpr,n,targetStr);
					break;
				case 9:
					results[9] = epsilonNearExpectedBpDistance(bpr,n,targetStr,epsilon);
					break;
				case 10:
					results[10] = epsilonNearEnsembleDefect(bpr,n,targetStr,epsilon);
					break;
				case 11:
					results[11] = RodrigoRobustness(bpr,n,seq);
					break;
				case 12:
					results[12] = ensembleHammingDistance(bpr,bpr2,n);
					break;
				case 13:
					results[13] = ensembleBasePairDistance(bpr,bpr2,n);
					break;
				case 14:
					results[14] = positionalEntropyDistance(bpr,bpr2,n);
					break;
				case 15:
					Hflag=1;
					break;
				case 16:
					MHflag=1;
					break;
				case 17:
					PRflag=1;
					break;
			}
		}
		if(normFlag!=0)
		{			
			for (int j=0;j<num;j++)
				results[outVec.at(j)]/=n;
		}
		if(verboseFlag == 0)
		{
			printf("%s\t%d\t%s\t%lf\t%lf\t",seq,n,targetStr,strEng,structureProbability(strEng, ensEng));
			for (int k=0;k<num;k++){
				if(outVec.at(k)<15)
					printf("%lf\t",results[outVec.at(k)]);
			}
			printf("\n");
		}
		else if(verboseFlag == 1)
		{
			printf("sequence: %s\nlength: %d\ntarget structure: %s\nenergy of target structure: %lf\nprobability of target structure: %lf\n",
			seq, n, targetStr,strEng,structureProbability(strEng,ensEng));
			for (int k=0;k<num;k++)
				if(outVec.at(k)<15)
					printf("%s: %lf\n",outTitle[outVec.at(k)].c_str(),results[outVec.at(k)]);
			printf("%s\n", "------------------------------------------------------------");
		}
		//if (posEntFlag==1)
			//printPositionalEntropy(bpr,bpr2,n,targetSeqFlag);
		if (Hflag==1)
			expectedHeight(bpr,n);
		if (MHflag==1)
			expectedModifiedHeight(bpr,n);
		if (PRflag==1)
			printTriProbabilities(bpr,n);
		}
	return 0;
}

void usage(char * name)
{
	printf("\nUSAGE: %s -s sequence -f fasta  [-d 0|2] [-c tragetStructure] [-t temperature] [-q tragetSequence] [-u turner99|turner04|andronescu07] [-x epsilon] [-o output] [-n] [-v]\ntry %s -h for help\n\n",name, name);
	exit(1);
}

