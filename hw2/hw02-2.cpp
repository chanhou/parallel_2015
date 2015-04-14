#include <iostream>
#include <fstream>
#include <time.h> 
using namespace std;

int main(int argc,char **argv){
	srand(time(NULL));

	int nr1 ;
	int nr2 ;
	int nc1 ;
	int nc2 ;
	
	ofstream file ( argv[5] );
	// read matrix size
	nr1 = stoi(argv[1]);
	nr2 = stoi(argv[3]);
	nc1 = stoi(argv[2]);
	nc2 = stoi(argv[4]);
	file<<nr1<<" "<<nc1<<" "<<nr2<<" "<<nc2<<endl;

	for(int i=0;i<nr1;i++){
		for(int q=0;q<nc1;q++){
			file<<(double)rand()/RAND_MAX *400000000;
			if (q != (nc1-1)) file<<" ";
		}
		file<<endl;
		
	}

	for(int i=0;i<nr2;i++){
		for(int q=0;q<nc2;q++){
			file<<(double)rand()/RAND_MAX *400000000;
			if (q != (nc2-1)) file<<" ";
		}
		if (i!=(nr2-1)) file<<endl;
	}

	file.close();

	return 0;
}

