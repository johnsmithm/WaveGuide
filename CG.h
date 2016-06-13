#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <vector>
#include <functional>

#define x first
#define y second
#define ABS(x) (x < 0 ? -x : x) 

void printTest(VecD u);

double residuum(SparceMatrix &A, vector<double> &u, vector<double> &f){
	double res = 0.,tmp = 0.;
	for(size_t i=0;i<u.size();++i){
		tmp = f[i];
		for(auto cell : A[i]){
			assert(cell.x>=0 && cell.x<u.size());
			tmp -= cell.y*u[cell.x];
		}
		res += tmp*tmp;
	}
	return sqrt(res);
}

void SolveEquationGS(SparceMatrix &A, vector<double> &u, vector<double> &f, double eps){
	for(;residuum(A,u,f)>eps;){
		for(size_t i=0;i<u.size();++i){
			double val = f[i], mid = 0.0001;
			//cerr<<f[i]<<" \n";
			for(auto cell : A[i]){
				//cerr<<"i:"<<cell.x<<" v:"<<cell.y<<" ";
				assert(cell.x>=0 && cell.x<u.size());
				if(cell.x != i)
					val -= cell.y*u[cell.x];
				else mid = cell.y;
				//cerr<<" val:"<<val<<" mid:"<<mid<<" ";
			}
			//cerr<<"\n";
			//cerr<<val<<" "<<mid<<" \n";
			u[i] = val/mid;
		}
		//printTest(u);
	}
}

vector<double> tmp;

void multiplyVectorMatrix(VecD &f, SparceMatrix A, VecD u){
	assert(f.size() == A.size());
	for(size_t i=0;i<A.size();++i){
		f[i] = 0.;
		for(auto cell : A[i]){
			assert(cell.x < u.size() && cell.x>=0);
			f[i] += cell.y*u[cell.x];
		}
	}
}

void multiplayVectorScalar(VecD &a, double s){
	for(size_t i=0;i<a.size();++i)
		a[i] *= s;
}

double eurclidianNorm(VecD a){
	double val = 0.;
	for(size_t i=0;i<a.size();++i)
		val += a[i]*a[i];
	return sqrt(val);
}

double multiplyVectorVector(VecD a, VecD b){
	assert(b.size() == a.size());
	double val = 0.;
	for(size_t i=0;i<a.size();++i)
		val += a[i]*b[i];
	return val;
}

void printTest(VecD a){
	for(size_t  i=0;i<a.size()&&i<10;++i)
		cerr<<a[i]<<" ";
	cerr<<"\n";
}

void inversePowerIteration(SparceMatrix A, SparceMatrix M, VecD f, VecD &u, double eps){
	double eigenvalue = 100., eigenvalueOld = 1.;
	tmp.assign(f.size(),0.);

	int nr = 0;
	while( ABS((eigenvalue-eigenvalueOld)/eigenvalueOld) > 0.0000000001){
		eigenvalueOld = eigenvalue;
		multiplyVectorMatrix(f,M,u);
		SolveEquationGS(A,u,f,eps);
		//printTest(u);        
		multiplayVectorScalar(u,1./eurclidianNorm(u));
		multiplyVectorMatrix(tmp,A,u);
		eigenvalue = multiplyVectorVector(u,tmp);
		multiplyVectorMatrix(tmp,M,u);
		eigenvalue /= multiplyVectorVector(u,tmp);
		//cout.setprecision(9);
		cerr<<"Step:"<<++nr<<" EigenValue:"<<eigenvalue<<"\n";
	}
	
	//cerr<<"EigenValue:"<<eigenvalue<<"\n";
}