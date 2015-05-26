#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>

#include "dnn.h"

using namespace std;

int main(int argc, char** argv){

	vector< vector <double> > my_array;
	vector< vector <double> > my_array_x;
	vector <double>  my_array_y;
	vector< vector <double> > val_x;
	vector <double>  val_y;
	// vector< vector< vector <double> > > eee;
	// vector<int> label_array;

	vector< vector <double> > testing;
	vector< vector <double> > test_x;
	vector <double>  test_y;
	vector< vector <double> > test_val_x;
	vector <double>  test_val_y;


	srand (time(NULL));	

	read_file( argv[1], my_array, 123);
	cv_split( my_array, my_array_x, my_array_y , val_x , val_y, 0.0);

	read_file( argv[2], testing, 123);
	cv_split( testing, test_x, test_y , test_val_x , test_val_y, 0.0);
	// cout<<my_array_x.size()<<" "<<my_array_x[0].size()<<endl;
	// dnn * net = new dnn( my_array_x.size() , my_array_x[0].size() , 2, 100,100);
	dnn * net = new dnn( my_array_x.size() , 123 , 2, 100,100);
	// dnn( int num, int dimension, int klass, int hidd1, int hidd2 );
	net->training(my_array_x, my_array_y, test_x, test_y, 1, 0.1, 0.001);




	return 0;
}