#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>

double ro;

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
	ro = stod(argv[1]);
	double eps = stod(argv[2]);
	readGrid(/*"unit_circle_fine.txt"*/);
	writeK(vertexes,"ksq.txt");
	getMassMatrix(M);
	writeSparceMatrix(M,"M.txt");
	getStiffnessMatrix(A);
	writeSparceMatrix(A,"A.txt");
    u.assign(A.size(),1.);
    f.assign(A.size(),0.);    
    inversePowerIteration( A, M, f, u, eps);
    writeVector(u,vertexes,"eigenmode.txt");

}
