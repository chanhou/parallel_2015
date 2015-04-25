#include <iostream>
#include <fstream>
#include <time.h> 
using namespace std;

int main(int argc,char **argv){
	srand(time(NULL));

	int dim ;
	int nPoints ;
	int nValues ;
	int nDivisions ;
	
	ofstream file ( argv[5] );
	// read matrix size
	dim = stoi(argv[1]);
	nPoints = stoi(argv[2]);
	nValues = stoi(argv[3]);
	nDivisions = stoi(argv[4]);
	file<<dim<<" "<<nValues<<" "<<nPoints<<" "<<nDivisions<<endl;

	for(int p=0;p<nPoints;p++){	
	
		for(int i=0;i<dim;i++){
			file<<(double)rand()/RAND_MAX *10<<" ";
		}

		for(int q=0;q<nValues;q++){
			file<<(double)rand()/RAND_MAX *1000;
			if (q != (nValues-1)) file<<" ";	
		}
		
		file<<endl;
	}


	file.close();

	return 0;
}

