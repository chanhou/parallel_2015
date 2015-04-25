#include <omp.h> 
#include <stdio.h> 
#include <iostream>
using namespace std;
int main() { 

while(true)	{
	#pragma omp parallel for
	for(int i=0; i<10; ++i){
		cout<<i<<endl;
	}
	// printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads()); 
	cout<<endl<<omp_get_max_threads();
}

return 0;

}
