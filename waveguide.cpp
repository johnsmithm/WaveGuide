#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>

#include "READ.h"
#include "CG.h"
#include "Timer.h"


using namespace std;


SparceMatrix A;
SparceMatrix M;
VecD u,f;

void run_test(){
	f.assign(3,0.); f[0]=3.; f[1]=9.; f[2] = -6.;
	u.assign(3,0.);
	A.assign(3,map<size_t,double>());
	A[0][0] = 4.; A[0][1] = -1.; A[0][2]=-1.;
	A[1][0] = -2.;A[1][1] = 6.; A[1][2] = 9.;
	A[2][0] = -1.; A[2][1] = 1.; A[2][2] = 7.;
	//SolveEquationGS(A,u,f);
	multiplyVectorMatrix(u,A,f);
	printTest(u);
	multiplayVectorScalar(f,2.);printTest(f);
	cerr<<eurclidianNorm(f)<<"\n";
	cerr<<multiplyVectorVector(f,u)<<"\n";

}

int main(int argc, char *argv[]){
	(void) argc; //to suppress Warnings about unused argc
	assert(argc>0);
	double ro = stod(argv[1]);
	double eps = stod(argv[2]);
	cerr<<ro<<'\n';  
	cerr<<eps<<'\n'; 
	//run_test();return 0;
	readGrid("inputs/unit_circle.txt");
	//getLocalMatrix(0);
	getMassMatrix(M);
	return 0;
    
    //getTestMatrix("reference-outputs/A-ref.txt",A);    
    //getTestMatrix("reference-outputs/M-ref.txt",M);


    u.assign(A.size(),1.);
    f.assign(A.size(),0.);    
    inversePowerIteration( A, M, f, u, eps);
	


}
