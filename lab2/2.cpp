#include <iostream>
#include <mpi.h>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;
unsigned long long fact[40];


int Ifsafe(int rank,unsigned long long i , int n){

  int d,perm[n];

  // generating ith permutation of n natural numbers
  for(int j = 0; j < n; ++j) {
    perm[j] = i / fact[n - 1 - j];
    i = i % fact[n - 1 - j];
  }

  for(int j = n - 1; j > 0; --j) {
    for(int k = j - 1; k >= 0; --k) {
      if(perm[k] <= perm[j]) perm[j]++;
    }
  }

  

  // checking diagonals is valid or not if diagonal is not valid this for loop will return 0 
  for(int j = 0; j < n; j++) {
    int val = perm[j];

    for(int k = j+1,  d = 1; k < n; k++, d++) {
      if(val - d == perm[k] || val + d  == perm[k]) return 0;
    
    }

    for(int k = j-1, d = 1; k >= 0; k--, d++) {
      if(val - d == perm[k] || val + d  == perm[k]) return 0;
    }

  }

 // if this permutation is valid then return 1
 return 1;
}



int main(int argc, char **argv)
{
	int n = atoi(argv[1]),result=0;
	int rank,n0_process,subtot=0,k=10;
	fact[0]=1;
	unsigned long long i;
	
	auto start = high_resolution_clock::now();

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&n0_process);

	// calculating factorial and storing in global array fact[] 
	#pragma omp parallel 
        for(i=1;i<=n;i++)
        {
            fact[i]=fact[i-1]*i;
        }
    
	    
	MPI_Barrier(MPI_COMM_WORLD);
    //  assigning equal number of permutating to every process
    for(i=rank;i< fact[n];i+=n0_process){
    	subtot += Ifsafe(rank,i,n); 
    }
	
    MPI_Reduce(&subtot, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); // adding subtotal of each process 
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)cout<<"\nTotal number of solution = "<<result<<"\n";

	MPI_Finalize();

	auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start);
  std::cout <<"time =" <<duration.count() << std::endl;
	
	return 0;
}