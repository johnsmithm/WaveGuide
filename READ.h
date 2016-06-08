#include <math.h>
#include <map>
#include <vector>
#include<fstream>
#include <string>


using namespace std;


typedef pair<double,double> pr;
typedef pair<double,pr > tr;
typedef vector<map<size_t,double>> SparceMatrix;
typedef vector<double> VecD;
#define x first
#define y second
#define pb push_back
#define mp make_pair


void getTestMatrix(string path, SparceMatrix &m){
	ifstream in(path);
	size_t q,p;
	double val;
	while(in>>q){
		in>>p>>val;
		assert(q<0);
		assert(p<0);
		while(m.size()<=q)
			m.push_back(map<size_t,double>());

		m[q][p] = val;
	}
	in.close();
}

vector<pr > vertexes;
vector<vector<int> > triangles;
int nrV,nrT;

void readGrid(string path){
	ifstream in(path);
	string s;
	in>>nrV;
	getline(in,s);
	getline(in,s);
	for(int i=0;i<nrV;++i){
		int q; double w,p;
		in>>q>>w>>p;
		vertexes.pb(mp(w,p));
	}
	in>>nrT;
	for(int i=0;i<18;++i)in>>s;

	for(int i=0;i<nrT;++i){
		vector<int> r(3,0);
		in>>r[0]>>r[1]>>r[2];
		triangles.pb(r);
	}
	in.close();
}

 #include "Source/Colsamm.h"
 using namespace _COLSAMM_;



 vector< vector< double > > getLocalMatrix(int i){

 	vector< vector< double > > my_local_matrix;
	vector<double> corners(6, 0.0);
	ELEMENTS::Triangle my_element;
 	// array corners contains the x- and y-coordinates of the
	// triangle corners in the order x0, y0, x1, y1, x2, y2
	corners[0] = vertexes[triangles[i][0]].x; corners[1] = vertexes[triangles[i][0]].y;
	corners[2] = vertexes[triangles[i][1]].x; corners[3] = vertexes[triangles[i][1]].y;
	corners[4] = vertexes[triangles[i][2]].x; corners[5] = vertexes[triangles[i][2]].y;
	// pass the corners to the finite element
 	my_element(corners);
 	my_local_matrix = my_element.integrate(v_() * w_());
 	/*for(int i=0;i<3;++i){
 		for(int j=0;j<3;++j)cerr<<my_local_matrix[i][j]<<" ";
 		cerr<<"\n";
 	}*/
 	return my_local_matrix;

 }


 void getMassMatrix(SparceMatrix &m){
 	vector< vector< double > > my_local_matrix;
 	m.assign(nrV,map<size_t,double>());


 	for(int i1=0;i1<nrT;++i1){
 		my_local_matrix = getLocalMatrix(i1);
 		for(int i=0;i<3;++i){
	 		for(int j=0;j<3;++j){
	 			assert(triangles[i1][i] < nrV && triangles[i1][j] < nrV);
	 		
	 			if(m[triangles[i1][i]].find(triangles[i1][j]) != m[triangles[i1][i]].end()){
	 				m[triangles[i1][i]][triangles[i1][j]] = m[triangles[i1][i]][triangles[i1][j]] + my_local_matrix[i][j];
	 			}
	 			else{
	 				m[triangles[i1][i]][triangles[i1][j]] = my_local_matrix[i][j];
	 			}
	 		}
	 		
	 	}
 	}
 	
 	for(auto cell : m[0])
 		cerr<<cell.x<<" "<<cell.y<<"\n";
 }