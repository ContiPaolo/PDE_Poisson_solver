# include <iostream>
using namespace std;
#include <cassert>
#include <vector>
#include <utility>
#include <cmath>
#include "Validation.hpp"

//To construct the tableau use it to populate the CCS vectors
 
typedef vector<int> V;

class Stencils{
public:
  int N; int M;
  int m_dim;
  int num_elem;
  V* tab;

  Stencils(int N1, int M1 ){
    N=N1; M= M1;
    m_dim = (N-1)*(M-1);
    num_elem= 5*(M-1)*(N-1)-2*(N+M)+4;
    tab = new V[m_dim]; assert(tab);
    for(int i=0; i <m_dim ; ++i){
      int sU = i-N+1;
      int sL = i-1;
      int sM = i;
      int sR = i+1;
      int sD = i+N-1;

      pair<int,int> v = bij(i, N, M);
      double n = v.first;
      double m = v.second;
      //Case where the ith element has 2 stencils: i.e. at the 4 corners.
      if (n==1 && m == M-1){
        tab[i].push_back(sM); tab[i].push_back(sR); tab[i].push_back(sD);
	assert(tab[i].size() == 3);
      }
      else if (n == N-1 && m == M-1){
	tab[i].push_back(sL); tab[i].push_back(sM); tab[i].push_back(sD);
	assert(tab[i].size() == 3);
      }
      else if (n==1 && m==1){
	tab[i].push_back(sU); tab[i].push_back(sM); tab[i].push_back(sR);
	assert(tab[i].size() == 3);
      }
      else if (n== N-1 && m ==1){
	tab[i].push_back(sU); tab[i].push_back(sL); tab[i].push_back(sM);
	assert(tab[i].size() == 3);
      }
      //Case where the ith element has 3 stencils:
      else if (m == M-1 && n!=1 && n!=N-1){ 
	tab[i].push_back(sL); tab[i].push_back(sM); tab[i].push_back(sR); tab[i].push_back(sD);	
	assert(tab[i].size() == 4);
      }
      else if (n==1 && m!=1 && m!= M-1){ //note: don't need to put the && since that case is before...
	tab[i].push_back(sU);tab[i].push_back(sM); tab[i].push_back(sR); tab[i].push_back(sD);
        assert(tab[i].size() == 4);
      }
      else if (m==1 && n!=1 && n!=N-1){
	tab[i].push_back(sU); tab[i].push_back(sL); tab[i].push_back(sM); tab[i].push_back(sR);
        assert(tab[i].size() == 4);
      }
      else if (n==N-1 && m!=1 && m!=M-1){
	tab[i].push_back(sU); tab[i].push_back(sL); tab[i].push_back(sM); tab[i].push_back(sD);
        assert(tab[i].size() == 4);
      } 
      else{ //Case where the ith element has 4 stencils:
	tab[i].push_back(sU); tab[i].push_back(sL); tab[i].push_back(sM); tab[i].push_back(sR); tab[i].push_back(sD);
        assert(tab[i].size() == 5);
      }
    }
  }
  //Prints the row/column indicies corresponding to the non-zero entries of the ith row/column
  void show_ith(int i){
    V::iterator it;
    for (it = tab[i].begin(); it !=tab[i].end(); ++it){
      cout <<*it << "  ";
    }
  }

  int get_elem(int i , int j){
    return tab[i][j];
  }
  int size_stencil(int i){
    return tab[i].size();
  }

//   //For testing when N and M are small: Prints all of the stencils in V:
//   void show_all(){
//     for (int i = 0; i<m_dim; ++i){
//       cout << "i="<< i << ": " ;
//       show_ith(i); cout << endl;
//     }
//   }
};

void CCS(Stencils & A, double * & Ax, int * & Ap, int * & Ai){
  int cur_num = 0;
  double h1 = 1./A.N;
  double h2 = 1./A.M;
  
  for(int i=0; i<A.m_dim; ++i){
    int ni = A.size_stencil(i);
    Ap[i+1] = Ap[i] +ni;

    for(int j=0; j<ni; ++j){
      int row_index = A.get_elem(i,j); // Extract from the ith vector in tab, the jth element
      Ai[cur_num] = row_index;
      //populating Ax:
      if(i == row_index){ //elem is on the diag
	Ax[cur_num] = 2.*(1./(h1*h1) + 1./(h2*h2)); }
      else if(row_index == i+1 || row_index == i-1){ //row_index is a left or right neighbour of i
	Ax[cur_num] = -1./pow(h1,2.0);}
      else {Ax[cur_num]= -1./pow(h2,2.0);} //j is a upper or lower neighbour of i
      ++cur_num;
  }
 }
  assert(Ap[A.m_dim] == A.num_elem);
}
// //For testing if N and M are small:
// void printCCS(double * & Ax, int * & Ap, int * & Ai, int nz, int dim){
//   cout << "Ap : ";
//   for(int i =0; i<dim+1; ++i){
//     cout << Ap[i] << "  " ;
//   }
//   cout << '\n' << "Ai: ";
//   for(int i = 0; i <nz; ++i){
//     cout <<Ai[i] << "  ";
//   }
//   cout << '\n'<< "Ax : ";
//   for(int i = 0; i<nz; ++i){
//     cout <<Ax[i] << "  ";
//   }
//   cout << endl;
// }
