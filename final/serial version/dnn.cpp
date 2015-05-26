#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

#include <fstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

#include <stack>

#include <algorithm>

#include <random> // random generator

#include "stopWatch.h"
#include "dnn.h"

#include <cassert>
#include <chrono>

// int forest_predict(double *attr);

using namespace std;

using std::string;

#define MAX_FEATURE (1024+1)

dnn::dnn( int num, int dimension, int klass, int hidd1, int hidd2  ) {
	N = num;
	dim = dimension;
	k = klass;
	h1 = hidd1;
	h2 = hidd2;
	// model['h'] = h # size of hidden layer 1
	// model['h2']= h2# size of hidden layer 2
	// model['W1']= 0.1 * np.random.randn(D,h)
	// model['b1'] = np.zeros((1,h))
	// model['W2'] = 0.1 * np.random.randn(h,h2)
	// model['b2']= np.zeros((1,h2))
	// model['W3'] = 0.1 * np.random.randn(h2,K)
	// model['b3'] = np.zeros((1,K))

	// std::random_device rd;
	// std::mt19937 gen(rd());
	// std::normal_distribution<> d(0,1.);
	// d(gen)

	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);

	std::normal_distribution<double> distribution (0.0,1.0);
	

	w1.resize(dim);
	for(auto &x: w1){
		for(int i=0;i <h1 ; ++i){
			x.push_back( 0.01*distribution(generator) );
			// cout<<0.1*distribution(generator)<<" ";
		}
	}

	// cout<<0.1*distribution(generator)<<endl;

	// b1.resize(1);
	// for(auto &x: b1){
	// 	x.resize( h1 ); // zero
	// }
	b1.resize(h1);
	
	w2.resize(h1);
	for(auto &x: w2){
		for(int i=0;i <h2 ; ++i){
			x.push_back( 0.01*distribution(generator) );
		}
	}
	// b2.resize(1);
	// for(auto &x: b2){
	// 	x.resize( h2 ); 
	// }
	b2.resize(h2);

	w3.resize(h2);
	for(auto &x: w3){
		for(int i=0;i < k ; ++i){
			x.push_back( 0.01*distribution(generator) );
		}
	}
	// b3.resize(1);
	// for(auto &x: b3){
	// 	x.resize( k ); 
	// }
	b3.resize(k);

}


void read_file(char *file, vector< vector <double> >& my_array , 
	int max_fea){

	double* features = new double[MAX_FEATURE];
	
	ifstream fin;
	string istring;
	fin.open(file);

	int max_index = -INFINITY;
	int num_data=0;
	while (std::getline(fin, istring)) {
		char *cstring, *tmp;
		int label;
		num_data++;
		// memset ( void * ptr, int value, size_t num );
		memset(features, 0, sizeof(double) * MAX_FEATURE);  // replace the value of the memory size

		cstring = new char[istring.size() + 1];
		// char * strncpy ( char * destination, const char * source, size_t num );
		strncpy(cstring, istring.c_str(), istring.size()+1);

		tmp =  strtok(cstring, ": "); // split by token
		label = atoi(tmp); // convert string to int
		tmp = strtok(NULL, ": "); // split by token

		while(tmp != NULL) {
			int id = atoi(tmp);
			if (id > max_index) max_index = id;
			tmp = strtok(NULL, ": ");
			features[id] = atof(tmp); // convert string to float
			tmp = strtok(NULL, ": ");
		}

		// cout<<features[14]<<endl;

		delete[] cstring;
	}
	delete[] features;

	// cout<< max_index<<" "<<num_data<<endl;

	fin.close();

	// double* my_array = new double[max_index*num_data];
	// int* label_array = new int[num_data];

	// vector< vector <double> > my_array(num_data, vector<double>(max_index,0));
	// vector< vector <double> > my_array;
	// vector<int> label_array(num_data,0);
	// vector<int> label_array;

	// my_array(num_data, vector<double>(max_index,0));
	my_array.resize(num_data);
	for(auto &x: my_array){
		// x.resize(max_index + 1); // insert label
		x.resize(max_fea +1);
	}

	// label_array.resize(num_data);

	fin.open(file);

	num_data = 0;

	while (std::getline(fin, istring)) {
		char *cstring, *tmp;
		int label;
		// memset ( void * ptr, int value, size_t num );
		// memset(my_array, 0, sizeof(double) * MAX_FEATURE);  // replace the value of the memory size

		cstring = new char[istring.size() + 1];
		// char * strncpy ( char * destination, const char * source, size_t num );
		strncpy(cstring, istring.c_str(), istring.size()+1);

		tmp =  strtok(cstring, ": "); // split by token
		label = atoi(tmp); // convert string to int
		tmp = strtok(NULL, ": "); // split by token

		// label_array[num_data] = label;
		if(label==-1)
			my_array[num_data][ 0 ] = 0;
		else
			my_array[num_data][ 0 ] = label;

		while(tmp != NULL) {
			int id = atoi(tmp);
			tmp = strtok(NULL, ": ");
			// my_array[ max_index * num_data + ( id - 1) ] = atof(tmp); // convert string to float
			my_array[ num_data ][ id ] = atof(tmp); // convert string to float
			tmp = strtok(NULL, ": ");
		}

		num_data ++;

		// cout<<features[14]<<endl;

		delete[] cstring;
	}

	fin.close();

	// cout<< my_array[0]<<endl;

	// for(auto &x : my_array[0]){
	// 	cout << x<<",";
	// }
	// cout<<endl;

	// for(auto &x : label_array){
	// 	cout << x<<",";
	// }
}

void cv_split( 
	std::vector< std::vector <double> > &my_array ,
	std::vector< std::vector <double> > &my_array_x,
	std::vector <double>  &my_array_y, 
	std::vector< std::vector <double> > &val_x ,
	std::vector <double>  &val_y, double percent){

	int num_train = my_array.size() * (1-percent);
	int num_val = my_array.size() - num_train;
	int max_index = my_array[0].size() -1;

	// randomize array
	random_shuffle ( my_array.begin(), my_array.end() );

	my_array_x.resize(num_train);
	for(auto &x: my_array_x){
		x.resize(max_index); 
	}

	val_x.resize(num_val);
	for(auto &x: val_x){
		x.resize(max_index); 
	}

	my_array_y.resize(num_train);
	val_y.resize(num_val);

	for(int i=0; i< my_array.size(); ++i){
		if (i<num_train)
			my_array_y[i] = my_array[i][0];
		else
			val_y[i - num_train] = my_array[i][0];

		for (int q = 0; q< max_index; ++q){
			if(i< num_train)
				my_array_x[i][q] = my_array[i][q+1];
			else
				val_x[i - num_train ][q] = my_array[i][q+1];
		}
	}

}

vector< vector< vector <double> > > split( vector< vector <double> >& my_array, int dim, double theta){
	
	vector< vector< vector <double> > > sppp(2);

	for(int i=0 ; i<(int) my_array.size(); ++i){
		if(my_array[i][dim] < theta ){
			sppp[0].push_back( my_array[i] );
		}
		else{
			sppp[1].push_back( my_array[i] );
		}
	}

	return sppp;
}

void dnn::training( vector< vector <double> >& train_x,
	vector <double> & train_y , 
	vector< vector <double> >& val_x,
	vector <double> & val_y , 
	int epochs, double yida, double reg){

	vector < vector <double>>  hidden_layer_1;
	vector < vector <double>>  hidden_layer_2;
	vector < vector <double>>  scores;

	vector < vector <double>>  temp;

	// vector <double>  correct_logprobs(train_x.size(), 0);
	vector < vector <double>>  dw3;
	vector < vector <double>>  dw2;	
	vector < vector <double>>  dw1;
	vector <double>  db3;
	vector <double>  db2;
	vector <double>  db1;

	vector < vector <double>>  dhidden2;
	vector < vector <double>>  dhidden1;

	double num_examples = (double)train_x.size();

	stopWatch timer;

	for(int iter=0; iter< epochs; ++iter){
		double t1;
		
		timer.start();
		/*
		feed forward propagation
		*/		
		hidden_layer_1 = ( matrix_matrix(train_x, w1, false, false) );
		m_v_add(hidden_layer_1, b1);
		sigmoid(hidden_layer_1);

		hidden_layer_2 = matrix_matrix(hidden_layer_1, w2, false, false);

		m_v_add(hidden_layer_2, b2);
		sigmoid(hidden_layer_2);

		scores = matrix_matrix(hidden_layer_2, w3, false, false);
		m_v_add(scores, b3);

		// exp_scores = np.exp(scores)
		// probs = exp_scores / np.sum(exp_scores, axis=1, keepdims=True) # [N x K]
		/*
		compute probability

		softmax
		#     e = T.exp(X)
		#     return e / T.sum(e, axis=1).dimshuffle(0, 'x')

        */
		for (int i=0 ; i<scores.size();++i){
			double temp = 0;
			for(int j=0; j< scores[0].size(); ++j){ // class num
				temp += exp(scores[i][j]); // normalize term
				// cout<<exp(scores[i][j])<<" ";
				// cout<<scores[i][j]<<" ";
			}
			// cout<<endl;
			// cout<<endl<<temp<<endl;
			for (int j=0; j< scores[0].size(); ++j){
				scores[i][j] = exp(scores[i][j]) / temp; // transform to probalitity
				// cout<<scores[i][j]<<" ";
			}
			// cout<<endl;
			// break;
		}
		// cout<<endl<<"======="<<endl;

		// # compute the loss: average cross-entropy loss and regularization
		// corect_logprobs = -np.log(probs[range(num_examples),y])
		// data_loss = np.sum(corect_logprobs)/num_examples
		// reg_loss = 0.5*reg*np.sum(W1*W1) + 0.5*reg*np.sum(W2*W2)+ 0.5*reg*np.sum(W3*W3)
		// loss = data_loss + reg_loss

		/*
		compute cost function
		*/
		// z      = [ 0.34, 0.66 ]
		// target = [ 1 , 0 ]
		// def multinominal_cross_entropy(z, target):
		//     loss = - T.mean( target * T.log(z) + (1 - target) * T.log(1 - z))
		// multinominal_cross_entropy(p_y_given_x, target) 

		double data_loss = 0;
		for(int i=0; i<scores.size();++i){
			data_loss += log(scores[i][train_y[i]]) ;
		}
		// for(int i=0; i<scores.size();++i){
		// 	if(train_y[i] == 0 ) // y = 0
		// 		data_loss += train_y[i] * log(scores[i][0]) + (1-train_y[i])*log(1-scores[i][0]);
		// 	else // y= +1
		// 		data_loss += train_y[i] * log(scores[i][1]) + (1-train_y[i])*log(1-scores[i][1]);
		// }
		data_loss = (-1)*data_loss / num_examples;
		
		double reg_loss = 0;

		temp = matrix_matrix(w1,w1,true,false);
		reg_loss += 0.5*reg*sum(temp);

		temp = matrix_matrix(w2,w2,true,false);
		reg_loss += 0.5*reg*sum(temp);

		temp = matrix_matrix(w3,w3,true,false);
		reg_loss += 0.5*reg*sum(temp);

		data_loss = data_loss + reg_loss;

		// for(auto &x: scores){
		// 	for (auto &y: x){
		// 		cout<< y<<" ";
		// 	}
		// 	cout<<endl;
		// }


		// # compute the gradient on scores
		for(int i=0;i<scores.size(); ++i){
			scores[i][ train_y[i] ] -= 1;
			for(int q=0;q<scores[0].size(); ++q)
				scores[i][q] /= num_examples;
		}

		// for(auto &x: scores){
		// 	for (auto &y: x){
		// 		cout<< y<<" ";
		// 	}
		// 	cout<<endl;
		// }



		// # BACKPROP HERE
		// dw3 = (hidden_layer_2.T).dot(dscores)
		// db3 = np.sum(dscores, axis=0, keepdims=True)

		// cout<<hidden_layer_2.size()<<" "<<hidden_layer_2[0].size()<<endl;
		// cout<<scores.size()<<" "<<scores[0].size()<<endl;

		dw3 = matrix_matrix(hidden_layer_2, scores, true,false);
		// for(auto &x: hidden_layer_2){
		// 	for (auto &y: x){
		// 		cout<< y<<" ";
		// 	}
		// 	cout<<endl;
		// }

		// cout<<dw3.size()<<" "<<dw3[0].size()<<endl;

		for(int i=0; i<scores[0].size(); ++i){
			db3.push_back(0);
			for(int q=0; q< scores.size(); ++q){
				db3[i] += scores[q][i];
			}
		}

		// cout<<db3.size()<<" "<<endl;
		// for(auto &x: db3){
		// 		cout<< x <<" ";
			
		// 	cout<<endl;
		// }

		// #backprop sigmoid nonlinearity here
		// dhidden2 = dscores.dot(w3.T)*sigmoid_grad(hidden_layer_2)

		// cout<<scores.size()<<" "<<scores[0].size()<<endl;
		// cout<<w3.size()<<" "<<w3[0].size()<<endl;
		dhidden2 = matrix_matrix( scores, w3, false, true );

		// cout<<dhidden2.size()<<" "<<dhidden2[0].size()<<endl;

		sigmoid_grad(hidden_layer_2);

		// cout<<hidden_layer_2.size()<<" "<<hidden_layer_2[0].size()<<endl;

		for(int i=0; i< dhidden2.size(); ++i){
			for(int q =0; q< dhidden2[0].size(); ++q){
				dhidden2[i][q] *= hidden_layer_2[i][q];
			}
		}

		// dw2 = (hidden_layer.T).dot(dhidden2)
		// db2 = np.sum(dhidden2, axis=0)
		dw2 = matrix_matrix(hidden_layer_1 , dhidden2, true , false);
		for(int i=0; i<dhidden2[0].size(); ++i){
			db2.push_back(0);
			for(int q=0; q< dhidden2.size(); ++q){
				db2[i] += dhidden2[q][i];
			}
		}

		// dhidden = dhidden2.dot(w2.T)*sigmoid_grad(hidden_layer_1)
		dhidden1 = matrix_matrix(dhidden2, w2, false, true);
		sigmoid_grad(hidden_layer_1);
		for(int i=0; i< dhidden1.size(); ++i){
			for(int q =0; q< dhidden1[0].size(); ++q){
				dhidden1[i][q] *= hidden_layer_1[i][q];
			}
		}

		// dw1 =  np.dot(X.T, dhidden)
		// db1 = np.sum(dhidden, axis=0)
		dw1 = matrix_matrix(train_x, dhidden1, true, false);

		for(int i=0; i<dhidden1[0].size(); ++i){
			db1.push_back(0);
			for(int q=0; q< dhidden1.size(); ++q){
				db1[i] += dhidden1[q][i];
			}
		}

		// # add regularization
		// dW3+= reg * W3
		// dW2 += reg * W2
		// dW1 += reg * W1
		for(int i=0; i<dw3.size(); ++i)
			for (int q=0; q<dw3[0].size(); ++q)
				dw3[i][q] += reg * w3[i][q];
		for(int i=0; i<dw2.size(); ++i)
			for (int q=0; q<dw2[0].size(); ++q)
				dw2[i][q] += reg * w2[i][q];
		for(int i=0; i<dw1.size(); ++i)
			for (int q=0; q<dw1[0].size(); ++q)
				dw1[i][q] += reg * w1[i][q];

		// cout<<"done"<<endl;

		// # update
		// W1 += -step_size * dW1
		// b1 += -step_size * db1
		// W2 += -step_size * dW2
		// b2 += -step_size * db2
		// W3 += -step_size * dW3
		// b3 += -step_size * db3

		// cout<<w1.size()<<" "<<w1[0].size()<<endl;
		// cout<<dw1.size()<<" "<<dw1[0].size()<<endl;

		// for(auto &x: w1){
		// 	for (auto &y: x){
		// 		cout<< y<<" ";
		// 	}
		// 	cout<<endl;
		// }
		// cout<<endl<<endl<<endl;
		// cout<<endl<<endl<<endl;
		// cout<<endl<<endl<<endl;
		// cout<<"************************"<<endl;
		// cout<<endl<<endl<<endl;
		update_2(w1, dw1, yida);
		// for(auto &x: w1){
		// 	for (auto &y: x){
		// 		cout<< y<<" ";
		// 	}
		// 	cout<<endl;
		// }
		update_1(b1, db1, yida);
		update_2(w2, dw2, yida);
		update_1(b2, db2, yida);
		update_2(w3, dw3, yida);
		update_1(b3, db3, yida);

		// cout<<yida<<endl;

		timer.stop();
		t1 = timer.elapsedTime();
		// cout<<t1<<endl;
		if (iter%1==0){
			printf("epochs %d: loss %f, val: %f, t: %f \n", 
				iter, data_loss , predict(val_x, val_y), t1 );
		}	

	}




}

double dnn::predict( std::vector< std::vector <double> >  & train_x,
		std::vector <double>   & train_y ){
	// # evaluate training set accuracy
	// if NONLINEARITY == 'RELU':
	//     hidden_layer = relu(np.dot(X, W1) + b1)
	//     hidden_layer2 = relu(np.dot(hidden_layer, W2) + b2)
	// elif NONLINEARITY == 'SIGM':
	//     hidden_layer = sigmoid(np.dot(X, W1) + b1)
	//     hidden_layer2 = sigmoid(np.dot(hidden_layer, W2) + b2)
	// scores = np.dot(hidden_layer2, W3) + b3
	// predicted_class = np.argmax(scores, axis=1)
	// print 'training accuracy: %.2f' % (np.mean(predicted_class == y))  

	vector < vector <double>>  hidden_layer_1;
	vector < vector <double>>  hidden_layer_2;
	vector < vector <double>>  scores;

	vector <int>  predict_class;
	
	
	// cout<<train_x.size()<<" "<<train_x[0].size()<<endl;
	// cout<<w1.size()<<" "<<w1[0].size()<<endl;

	hidden_layer_1 = ( matrix_matrix(train_x, w1, false, false) );
	// cout<<"done"<<endl;
	// cout<<hidden_layer_1.size()<<" "<<hidden_layer_1[0].size()<<endl;
	// cout<<b1.size()<<" "<<endl;

	m_v_add(hidden_layer_1, b1);
	// cout<<"done"<<endl;
	sigmoid(hidden_layer_1);

	// cout<<"done"<<endl;

	hidden_layer_2 = matrix_matrix(hidden_layer_1, w2, false, false);
	m_v_add(hidden_layer_2, b2);
	sigmoid(hidden_layer_2);

	// for(auto&x:hidden_layer_2){
	// 	for(auto&z:x){
	// 		cout<<z<<" ";
	// 	}
	// }

	// cout<<hidden_layer_2.size()<<" "<<hidden_layer_2[0].size()<<endl;
	// cout<<w3.size()<<" "<<w3[0].size()<<endl;
	scores = matrix_matrix(hidden_layer_2, w3, false, false);
	m_v_add(scores, b3);

	// for(auto&x:scores){
	// 	for(auto&z:x){
	// 		cout<<z<<" ";
	// 	}
	// }

	
	for(int i=0; i< scores.size(); ++i){
		predict_class.push_back(0);
		double temp = -INFINITY;
		for (int q=0; q<scores[0].size(); ++q){
			// cout<<scores[i][q]<<endl;
			if (scores[i][q]>temp){
				temp = scores[i][q];
				predict_class[i] = q; // y label
			}
		}
	}

	double correct =0;
	for (int i=0; i<predict_class.size(); ++i){
		// cout<<train_y[i]<<" ";
		if (predict_class[i] == train_y[i])
			correct++;
	}

	return correct/(double)predict_class.size();

}

void sigmoid (std::vector< std::vector <double> >& array){
	// x = 1/(1+np.exp(-x))
	for (int x=0; x < array.size(); ++x){
		for(int y = 0; y<array[0].size(); ++y){
			array[x][y] = 1 / ( 1 + exp(-1*array[x][y]) );
		}
	}
}

void sigmoid_grad(std::vector< std::vector <double> >& array){
	// (x)*(1-x)
	for (int x=0; x < array.size(); ++x){
		for(int y = 0; y<array[0].size(); ++y){
			array[x][y] = array[x][y] * ( 1- array[x][y] );
		}
	}
}

double sum(std::vector< std::vector <double> >& array1){
	double sum=0;
	for (int i=0;i<array1.size(); ++i){
		for(int q=0; q< array1[0].size(); ++q){
			sum += array1[i][q];
		}
	}
	return sum;
}

void update_1 ( std::vector <double> & array1,
	std::vector <double> & array2,
	double b ){

	for (int i=0; i< array1.size(); ++i){
			array1[i] += -b* array2[i] ;
	}
}

void update_2 (std::vector< std::vector <double> >& array1,
	std::vector< std::vector <double> >& array2,
	double b ){
	// cout<<"done"<<endl;
	// cout<<array1.size()<<" "<<array1[0].size();
	// cout<<array2.size()<<" "<<array2[0].size();
	for (int i=0; i< array1.size(); ++i){
		for(int q =0 ;q < array1[0].size(); ++q){
			array1[i][q] += -b * array2[i][q] ;
		}
	}
}



void m_v_add (std::vector< std::vector <double> >& array,
	std::vector <double> & b ){

	assert(array[0].size()==b.size());
	// cout<<"==="<<endl;
	for (int i=0; i< array.size(); ++i){ // 20
		for(int q =0 ;q < array[0].size(); ++q){ //100
			array[i][q] += b[q];
		}
	}
}

std::vector< std::vector <double> > matrix_matrix(std::vector< std::vector <double> >& array1, 
	std::vector< std::vector <double> >& array2, bool transpose_1, bool transpose_2){

	int x = array1.size();
	int y = array2[0].size();

	if (transpose_1){
		
		int z = array1[0].size();
		vector< vector <double> > array1_T(z, vector<double> (x,0));
		for(int i=0; i<x; ++i ){
			for(int j=0; j<z; ++j){
				array1_T[j][i] = array1[i][j];
			}
		}

		// array1 transpose
		assert(array1.size()==array2.size());

		// cout<<x<<" "<<array1[0].size()<<" "<<array2.size()<<" "<<y<<endl;
		vector< vector <double> > result(z, vector<double> (y,0));
		for (int i=0; i < array1_T.size(); ++i){
			for(int k = 0; k<array2[0].size(); ++k){
				for (int j=0; j< array1_T[0].size(); ++j){
					result[i][k] += array1_T[i][j] * array2[j][k];
				}
			}
		}

		array1_T.clear();
		return result;
	}
	else if (transpose_2){
		int z = array2.size();
		vector< vector <double> > array2_T(y, vector<double> (z,0));
		for(int i=0; i<z; ++i ){
			for(int j=0; j<y; ++j){
				array2_T[j][i] = array2[i][j];
			}
		}

		// array 2 transpose
		assert(array1[0].size()==array2[0].size());
		// cout<<x<<" "<<array1[0].size()<<" "<<array2.size()<<" "<<y<<endl;
		vector< vector <double> > result(x, vector<double> (z,0));
		for (int i=0; i < array1.size(); ++i){
			for(int k = 0; k<array2_T[0].size(); ++k){
				for (int j=0; j< array1[0].size(); ++j){
					result[i][k] += array1[i][j] * array2_T[j][k];
				}
			}
		}

		array2_T.clear();
		return result;
	}
	else{

		assert(array1[0].size()==array2.size());

		vector< vector <double> > result(x, vector<double> (y,0));
		for (int i=0; i < array1.size(); ++i){
			for(int k = 0; k<array2[0].size(); ++k){
				for (int j=0; j< array1[0].size(); ++j){
					result[i][k] += array1[i][j] * array2[j][k];
				}
			}
		}
		return result;
	}

	
	// cout<<result.size()<<" "<<result[0].size()<<endl;
	
}


// void sort_index( vector< vector <double> >& my_array , int colindex ){

// 	// bool compareTwoRows(vector< double > rowA, vector< double > rowB){
// 	// 	return (  (rowA[ colindex ]<rowB[ colindex ]) );
// 	// // return ( (rowA[0]<rowB[0]) || ((rowA[0]==rowB[0])&&(rowA[1]<rowB[1])) );
// 	// }

// 	auto glambda = [ colindex ]( vector< double > a, vector< double > b) { return a[colindex] < b[colindex] ; };

// 	sort(my_array.begin(), my_array.end(), glambda );

// }

// tree_node* build_tree(std::vector< std::vector <double> > &my_array, 
// 	int colindex , double theta, double epsilon){

// 	// if terminate 
// 	// return ans

// 	// else 
// 	// find decision of theta and dim
// 	// split data
// 	// run build dnn
// 	// create a node 
// 	// assign left node and right node
// 	// return the node


// 	int *ggg = check_terminate(my_array,   colindex , theta, epsilon);
// 	// cout<<ggg[0]<<" "<<ggg[1]<<endl;
// 	if(ggg[0]==1){
// 		// cout<<"hi"<<endl;
// 		if( ggg[1]==0 ){
// 			// cout<<"hi hi"<<endl;
// 			delete[] ggg;
// 			return NULL;
// 		}
// 		else{
// 			// cout<<"hi~"<<endl;
// 			tree_node* ddd  = new tree_node;
// 			ddd->answer = ggg[1];
// 			ddd->left = NULL;
// 			ddd->right = NULL;

// 			delete[] ggg;
// 			return ddd;
// 		}
// 	}
// 	else{
// 		// cout<<"no hi"<<endl;
// 		vector< vector< vector <double> > > eee;

// 		double * abc = new double [2];

// 		// cout<<"total "<<my_array.size()<<endl;

// 		abc = decision(my_array);

// 		// // split array

// 		tree_node* ddd  = new tree_node;

// 		// cout<<abc[0]<<endl;

// 		ddd -> dim = abc[0];
// 		ddd -> theta = abc[1];

// 		eee = split(my_array, abc[0], abc[1]);

// 		for(int oo=0;oo< my_array.size(); ++oo){
// 			my_array[oo].clear();
// 			vector<double> ().swap(my_array[oo]);
// 		}
// 		my_array.clear();
// 		vector< vector <double> > ().swap( my_array);

// 		// cout<<"eee "<<endl;

// 		// cout<<eee[0].size()<<endl;
// 		// cout<<eee[1].size()<<endl;

// 		// tree_node * fff ;

// 		// wrong!!!!
// 		// for (int qq =0;qq<2; ++qq){
// 		// 	// if( eee[qq].size()==0 )
// 		// 	fff = build_tree(eee[qq], abc[0] , abc[1], epsilon);
// 		// 	if(qq==0)
// 		// 		ddd->left = fff;
// 		// 	else
// 		// 		ddd->right = fff;
// 		// }

// 		ddd->left = build_tree(eee[0], abc[0] , abc[1], epsilon);
// 		ddd->right = build_tree(eee[1], abc[0] , abc[1], epsilon);

// 		for(int ff=0;ff<2;++ff){
// 			for( int i=0; i < (int)eee[ff].size(); i++){
// 				eee[ff][i].clear();
// 				vector<double> ( ).swap(eee[ff][i]);
// 			}
// 			eee[ff].clear();
// 			vector< vector< double> >( ).swap(eee[ff]);
// 		}
// 		eee.clear();
// 		vector< vector< vector< double> > >( ).swap(eee);


// 		delete[] ggg;
// 		delete[] abc;

// 		return ddd;		
// 	}


// }

// double confusion(  double a , double b ){
// 	return 1- pow(a/(a+b),2) - pow(b/(a+b),2);
// }



// double * decision(std::vector< std::vector <double> >& array){

// 	double * abc = new double(2);

// 	double best_theta = 0;
// 	double best_dim = -1;
// 	double best_error = INFINITY;
// 	double error;
// 	vector< vector< vector <double> > > sppp;

// 	// cout<<" array size "<<array.size()<<endl;

// 	// for(auto x:array){
// 	// 	cout<<x[1]<<" "<<x[2]<<endl;
// 	// }
// 	// cout<<endl;

// 	for(int i=1; i < (int)(array[0].size() ) ;++i){ // dimension 
// 		sort_index(array, i);
// 		// for(auto x:array){
// 		// 	cout<<x[1]<<" "<<x[2]<<endl;
// 		// }
		

// 		double theta_prepare;
// 		// cout<<"dim "<<i<<endl;

// 		for (int row=0; row < (int)(array.size()+1 ); ++row ){ // for each data
// 			if (row == 0){
// 				// cout<<"000 "<<endl;
// 				theta_prepare = ( (array[row][i]) + (-1.) ) /2.0;
// 			}
// 			else if ((row) == (int)array.size()  ){
// 				// cout<<"111 "<<endl;
// 				theta_prepare = ( (array[row-1][i]) + ( 1.) ) /2.0;
// 			}
// 			else{
// 				// cout<<"222 "<<array[row-1][i]<<" "<<array[row][i]<<endl;

// 				if( array[row-1][i] == array[row][i] ) continue;

// 				theta_prepare = ( (array[row-1][i]) + (array[row][i]) ) /2. ;
// 			}

// 			sppp = split( array, i, theta_prepare);

// 			// cout<<"theta_prepare "<<theta_prepare<<endl;
// 			// cout<<"size "<<sppp[0].size()<<endl;
// 			// cout<<sppp[1].size()<<endl;

// 			// if(sppp[0].size()==0 || sppp[1].size()==0){
// 			// 	cout<<"+++++++++++++++"<<endl;
// 			// 	cout<<array.size()<<" "<<array[0].size() <<" "<< i<<" "<<theta_prepare<<endl;
// 			// }

// 			error = 0;

// 			for (int q=0; q<2; q++){
// 				double c=0., d= 0.;

// 				for(auto &y: sppp[q]){
// 					if((int)y[0]==1) c++;
// 					else d++;
// 				}
// 				// if (c+d == sppp[q].size()) cout<<"00000000000000"<<endl;
// 				// cout<<sppp[q].size()<<endl;
// 				// cout<<"co "<<c<<" mi "<<d<<endl;

// 				if (sppp[q].size()==0) error += 0;
// 				else error += (c+d)/((double)array.size())*confusion(c,d);
// 				// error += (c+d)/((double)array.size())*confusion(c,d);
// 			}

// 			// cout<<"error "<<error<<endl;

// 			if( error < best_error) {
// 				// cout<<"=================="<<endl;
// 				// cout<<"dim "<<i<<" error "<<error<<endl;
// 				best_error = error;
// 				best_theta = theta_prepare;
// 				best_dim = i;
// 			}

// 			for(int ff=0;ff<2;++ff){
// 				for( int i=0; i < (int)sppp[ff].size(); i++){
// 				sppp[ff][i].clear();
// 				vector<double> ( ).swap(sppp[ff][i]);
// 				}
// 				sppp[ff].clear();
// 				vector< vector< double> >( ).swap(sppp[ff]);
// 			}
// 			sppp.clear();
// 			vector< vector< vector< double> > >( ).swap(sppp);


// 		}
// 	}

// 	abc[0] = best_dim;
// 	abc[1] = best_theta;

// 	// cout<<abc[0]<<endl;
// 	// cout<<"best error "<<best_error<<endl;
// 	// cout<<"best "<<abc[0]<<" "<<abc[1]<<endl;

// 	return abc;
// }

// int * check_terminate( std::vector< std::vector <double> >& my_array,  
// 	int dim , double theta ,double epsilon ){

// 	int * ggg = new int (2);

// 	int plus=0;
// 	int minus=0;
// 	for(auto &y: my_array){
// 		if((int)y[0]==1) plus ++;
// 		else minus ++;
// 	}
	
// 	// int count = 0;
// 	// bool waa = false;

// 	// for(int tt=0; tt< my_array.size()-1;++tt){
// 	// 	// cout<<y[dim]<<",";
// 	// 	if ( my_array[tt][dim] == my_array[tt+1][dim] ) count++;
// 	// }

// 	// if (count == ((int)my_array.size()-1)) waa = true;

// 	// cout<<"waa "<<waa<<endl;
	
// 	// cout<<plus<<" "<<minus<<endl;	
// 	// cout<<confusion(plus,minus)<<endl;

// 	// if(my_array.size() == 3){
// 	// 	cout<<"---> "<<plus<<" "<<minus<<endl;
// 	// 	cout<<"---> "<<confusion(plus,minus)<<endl;
		
// 	// 	for(auto &y: my_array){
// 	// 		cout<<y[dim]<<",";
// 	// 	}
// 	// 	cout<<endl;
// 	// }

// 	if (my_array.size() == 0){
// 		// cout<<"Null happen"<<endl;
// 		// return ;

// 		// cout<<"null00 "<<plus<<" "<<minus<<endl;

// 		ggg[0] = 1;
// 		ggg[1] = 0;
// 	}
// 	else if ( plus== 0 || minus == 0 ){
// 		// tree_node* ddd  = new tree_node;
// 		// cout<<"plus or minus happen"<<endl;
// 		ggg[0] = 1;
// 		if( plus == 0)
// 			ggg[1] = -1;
// 		else
// 			ggg[1] = 1;

// 		// cout<<"gggg "<<plus<<" "<<minus<<endl;

// 		// ddd->left = NULL;
// 		// ddd->right = NULL;

// 		// return ddd;
// 	}
// 	else if( confusion(plus,minus) <= epsilon ){ //  termination criteria
// 		// tree_node* ddd  = new tree_node;
// 		// cout<<"confusion happen"<<endl;
// 		ggg[0] = 1;
// 		if(plus>minus)
// 			ggg[1] = 1;
// 		else
// 			ggg[1] = -1;

// 		// cout<<"confusion "<<plus<<" "<<minus<<endl;
// 		// ddd->left = NULL;
// 		// ddd->right = NULL;

// 		// return ddd;
// 	}
// 	else if (plus == minus ){
// 		ggg[0] = 1;
// 		if(rand() %2 == 0) ggg[1] = 1;
// 		else ggg[1] = -1;
// 	}
// 	// else if (waa){
// 	// 	ggg[0] = 1;
// 	// 	if(plus>minus)
// 	// 		ggg[1] = 1;
// 	// 	else
// 	// 		ggg[1] = -1;
// 	// }
// 	else{
// 		ggg[0] = 0;
// 		ggg[1] = 0;
// 	}

// 	return ggg;


// }