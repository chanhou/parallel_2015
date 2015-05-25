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

	srand (time(NULL));

	read_file( argv[1], my_array);

	cv_split( my_array, my_array_x, my_array_y , val_x , val_y, 0.1);

	dnn * net = new dnn( my_array_x.size() , my_array_x[0].size() , 2, 10,10);
	// dnn( int num, int dimension, int klass, int hidd1, int hidd2 );
	net->training(my_array_x, my_array_y, 1, 0.1);


	// // double * abc ;

	// // sort_index(my_array, 2);
	// // for(auto x:my_array){
	// // 	cout<<x[1]<<" "<<x[2]<<endl;
	// // }

	// // cout<<"read data done"<<endl;

	// // abc = decision(my_array);

	// // cout<<abc[0]<<" "<<abc[1]<<endl;

	// tree* ddd = new tree();	

	// ddd -> root -> dim = abc[0];
	// ddd -> root -> theta = abc[1];

	// eee = split(my_array, abc[0], abc[1]);

	// // cout<<eee[0].size()<<endl;
	// // cout<<eee[1].size()<<endl;
	// // cout<<my_array.size()<<endl;
	// // cout<<atof(argv[2])<<endl;

	// // tree_node * sss;

	// // for(int rr=0;rr<2;rr++){
	// // 	sss = build_tree(eee[rr], abc[0], abc[1], atof(argv[2]));
	// // 	if(rr==0)
	// // 		ddd->root->left = sss;
	// // 	else
	// // 		ddd->root->right = sss;
	// // }
	// ddd->root->left = build_tree(eee[0], abc[0], abc[1], atof(argv[2]));
	// ddd->root->right = build_tree(eee[1], abc[0], abc[1], atof(argv[2]));

	// // do one time for branching
	// // then use
	// // left = build_tree()
	// // right = build_tree()
	// // to build the tree

	// FILE *writefl;
	// writefl = fopen ( "tree_pred_func.cpp", "w");

	// ddd->print_file(writefl);


	return 0;
}