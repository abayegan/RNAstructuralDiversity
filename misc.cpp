#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include <ctype.h>
#include <vector>
#include <string.h>
#include <fstream>
#include"misc.h"

#ifdef _WIN32
	#include <windows>
#elif defined(__APPLE__) || defined(__linux)  || defined(__unix)  || defined(__posix) 
	#include <unistd.h>
	#include <limits.h>
#endif

using namespace std;

void parseFasta(string fname, vector<string> &seqs, vector<string> &strs){
	ifstream input(fname.c_str());
    if(!input.good()){
        cerr << "Error opening '"<<fname<< endl;
    }
    string line, name, content,str;
    while(getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){
            if( !name.empty() ){
                seqs.push_back(content);
                if (!str.empty())
					strs.push_back(str);
				else
					strs.push_back("MFE");
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
            str.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } 
            else if (line[0]=='(' || line[0]==')' || line[0]=='.' )
				str += line;
			else
                content += line;
            }
        }
    if( !name.empty() ){ // Print out what we read from the last entry
       seqs.push_back(content);
       if (!str.empty())
			strs.push_back(str);
		else
			strs.push_back("MFE");
    }
	return;
}

/*The following function check if the input string is a valid RNA sequence*/
int CheckSequence(char * sequence){     //safe
  int i;
  for (i=1;i<strlen(sequence);++i)
    {
     if (toupper(sequence[i])!='A' && toupper(sequence[i])!='U' && toupper(sequence[i])!='C' && toupper(sequence[i])!='G'&& toupper(sequence[i])!='T') //check if there are invalid characters
        {
         printf("The input string should only contain A,U,T,C,G!\n");
         exit(1);
        }
     else if (toupper(sequence[i])=='T') //change T to U
        sequence[i]='U';
     else
        sequence[i]=toupper(sequence[i]); // change lower case to upper case
    }
  return 0;
}
char* getExecPath(char* argv0){
	char* exec_path=(char*) malloc (sizeof(char)*PATH_MAX);
	exec_path[0]='.';
	exec_path[1]='\0';
	int path_length=1;
	FILE* file;
	#if defined(_WIN32)
	/* GetModuleFileName */
	/* _pgmptr */
	#elif defined(__APPLE__) || defined(__linux)  || defined(__unix)  || defined(__posix) 
	char buff[PATH_MAX];
	int bufsize = PATH_MAX-1;
	if(file = fopen("/proc/self/exe", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/self/exe", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}		
	}
	else if(file = fopen("/proc/curproc/file", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/curproc/file", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}		
	}
	else if(file = fopen("/proc/self/path/a.out", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/self/path/a.out", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}		
	}
	else{
		exec_path= argv0;
	}
	char* slash_pointer= strrchr(exec_path,'/');
	int slash_position = (int)(slash_pointer - exec_path);
	if(slash_position != path_length){			
		exec_path[slash_position]='\0';
	}
	else{
		exec_path[0]='.';
		exec_path[1]='\0';
	}

	#endif
	

	return exec_path;
}


