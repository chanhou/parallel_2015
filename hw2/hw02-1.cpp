#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <cmath>

#include "stopWatch.h"

// #include <typeinfo>

using namespace std;

void split(const string& s, char c, vector<string>& v) {
	string::size_type i = 0;
	string::size_type j = s.find(c);

	if (j != string::npos){
		while (j != string::npos) {
			//if(j==(s.length())) break;
			//if(s.substr(i, j-i).find_first_not_of(' ') == string::npos) break;
			//cout<<"hh"<<endl;
			v.push_back(s.substr(i, j-i));
			i = ++j;
			j = s.find(c, j);
			//cout<<"s length"<<s.length()<<endl;
			//cout<<"j"<<j<<endl;
			
			if(j==(s.length()-1)) {
				v.push_back(s.substr(i, j-i));
				break;
			}
			

			if (j == string::npos)
				v.push_back(s.substr(i, s.length()));
			
		}	
	}
	else {
		// cout<<"hh"<<endl;
		v.push_back(s.substr(0, s.length()));
	}
}

void print(double** a, int row, int column){
	for (int i=0;i<row;i++){
		for (int j=0;j<column;j++){
			cout<<a[i][j]<<" ";
		}
		cout<<endl;
	}
}

int main(int argc,char **argv){

	string sa;
	string line;
	// vector<string> v;
	bool start = true;

	double **a;
	double **b;
	int nr1 ;
	int nr2 ;
	int nc1 ;
	int nc2 ;
	int count_a =0;
	int count_b =0;
	vector<string> v;

	
	ifstream file;
	file.open( argv[1] );
	
	while ( getline (file,line) ){
		v.clear();
		// get the number in the string by split ' '
		split(line,' ',v);
		// cout<<line<<endl;
		// cout<<(count_b+count_a)<<endl;
		if (start){ // initialize, allocate matrix memory

			// read matrix size
			nr1 = stoi(v[0]);
			nr2 = stoi(v[2]);
			nc1 = stoi(v[1]);
			nc2 = stoi(v[3]);
			count_a = 0;
			count_b = 0;

			// allocate array
			a = new double *[nr1];
			b = new double *[nr2];

			for (int i = 0;i<nr1;i++){
				a[i] = new double[nc1];
			}
			for (int i = 0;i<nr2;i++){
				b[i] = new double[nc2];
			}

			start = false;
			// int a[][];
		}
		else{
			// read number into matrix
			if (nr1>count_a){ // first matrix
				// if (v.size()==0){ // one number only
				// 	// cout<<line<<endl;
				// 	a[count_a][0] = stod(line);
				// 	// cout<<"a"<<a[count_a][0]<<endl;
				// }
				// else{
					for (int i=0;i<(int)v.size();++i){
						// cout << v[i] << '\n';	
						a[count_a][i] = stod(v[i]);
						//cout<<"a"<<a[count_a][i]<<endl;
					}		
				// }
				count_a ++;
			}
			else{ // second matrix
				// if (v.size()==0){
				// 	// cout<<line<<endl;
				// 	b[count_b][0] = stod(line);
				// 	// cout<<"b"<<b[count_b][0]<<endl;
				// }
				// else{
					for (int i=0;i<(int)v.size();++i){
						// cout << v[i] << '\n';	
						b[count_b][i] = stod(v[i]);
						//cout<<"b"<<b[count_b][i]<<endl;
					}		
				// }
				count_b ++;
			}
			// cout<<(count_b+count_a)<<endl;
			if ((count_b+count_a) == (nr1+nr2)) break;
		}
	}
	file.close();

	// cout<<endl<<"array A"<<endl;
	
	// print(a,nr1,nc1);

	// cout<<endl<<"array B"<<endl;

	// print(b,nr2,nc2);

	// initialize result
	double result[nr1][nc2];
	double temp;
	double flops = 0;
	double mem_r = 0;
	double mem_w = 0;
	double mem_size =(nr1*nc1+nr2*nc2+nr1*nc2)*8/1024.;
	double norm = 0;
	double t = 0;
	int count = 0;
	double secs;
	stopWatch timer;

	// doing matrix product

	do{
		timer.start();
		mem_r = 0;
		mem_w = 0;
		flops = 0;
		for(int i =0;i<nr1;i++){
			for(int j=0;j<nc2;j++){
				temp = 0;
				for (int k=0;k<nc1;k++){
					temp += a[i][k]*b[k][j];
					mem_r += 2;
					flops += 2;
				}
				result[i][j] = temp;
				mem_w += 1;
			}
		}
		timer.stop();
		t += timer.elapsedTime();	
		count++;
	} while (count<=50);

	// cout<<endl<<"result"<<endl;
	for (int i=0;i<nr1;i++){
		for (int j=0;j<nc2;j++){
			//cout<<result[i][j]<<" ";
			norm += pow(result[i][j],2);
		}
		// cout<<endl;
	}
	norm = pow(norm,0.5);


	secs = (double) (t/count);
	// cout<<count<<" "<<t<<endl;

	cout<<nr1<<", "<<nc1<<", "<<nr2<<", "<<nc2<<", ";
	cout<<mem_size<<", "<<norm<<", "<<secs<<", "<<flops<<", ";
	cout<<mem_r<<", "<<mem_w<<", "<< flops/secs <<", ";
	cout<<mem_r*8/1024./1024./secs<<", "<<mem_w*8/1024./1024./secs;
	cout<<endl;


	// cout<<endl<<"result"<<endl;
	// cout<<"FPOperation:"<<flops<<endl;
	// cout<<"mem_r:"<<mem_r<<endl;
	// cout<<"mem_w:"<<mem_w<<endl;

	return 0;
}
