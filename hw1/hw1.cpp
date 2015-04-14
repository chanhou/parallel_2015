#include <iostream>
#include <fstream>

#include <string>
#include <vector>

// #include <typeinfo>

using namespace std;

void split(const string& s, char c, vector<string>& v) {
	string::size_type i = 0;
	string::size_type j = s.find(c);

	while (j != string::npos) {
		v.push_back(s.substr(i, j-i));
		i = ++j;
		j = s.find(c, j);

		if (j == string::npos)
			v.push_back(s.substr(i, s.length()));
	}
}

void print(int** a, int row, int column){
	for (int i=0;i<row;i++){
		for (int j=0;j<column;j++){
			cout<<a[i][j]<<" ";
		}
		cout<<endl;
	}
}

int main(){

	string sa;
	string line;
	// vector<string> v;
	bool start = true;

	int **a;
	int **b;
	int row_a ;
	int row_b ;
	int column_a ;
	int column_b ;
	int count_a ;
	int count_b ;

	cout << "Pls enter file name"<<endl;

	// input file name
	cin >> sa ;
	
	ifstream file;
	file.open( sa+".txt" );
	
	if (file.is_open()){
		while ( getline (file,line) ){
			vector<string> v;
			// get the number in the string by split ' '
			split(line,' ',v);
			cout<<line<<endl;
			if (start){ // initialize, allocate matrix memory

				// read matrix size
				row_a = stoi(v[0]);
				row_b = stoi(v[2]);
				column_a = stoi(v[1]);
				column_b = stoi(v[3]);
				count_a = 0;
				count_b = 0;

				// allocate array
				a = new int *[row_a];
				b = new int *[row_b];

				for (int i = 0;i<row_a;i++){
					a[i] = new int[column_a];
				}
				for (int i = 0;i<row_b;i++){
					b[i] = new int[column_b];
				}

				start = false;
				// int a[][];
			}
			else{
				// read number into matrix
				if (row_a>count_a){ // first matrix
					if (v.size()==0){ // one number only
						// cout<<line<<endl;
						a[count_a][0] = stoi(line);
						// cout<<"a"<<a[count_a][0]<<endl;
					}
					else{
						for (int i=0;i<(int)v.size();++i){
							// cout << v[i] << '\n';	
							a[count_a][i] = stoi(v[i]);
							// cout<<"a"<<a[count_a][i]<<endl;
						}		
					}
					count_a ++;
				}
				else{ // second matrix
					if (v.size()==0){
						// cout<<line<<endl;
						b[count_b][0] = stoi(line);
						// cout<<"b"<<b[count_b][0]<<endl;
					}
					else{
						for (int i=0;i<(int)v.size();++i){
							// cout << v[i] << '\n';	
							b[count_b][i] = stoi(v[i]);
							// cout<<"b"<<b[count_b][i]<<endl;
						}		
					}
					count_b ++;
				}
			}
		}
		file.close();
	}

	cout<<endl<<"array A"<<endl;
	
	print(a,row_a,column_a);

	cout<<endl<<"array B"<<endl;

	print(b,row_b,column_b);

	// initialize result
	int result[row_a][column_b];
	int temp;
	int fpo = 0;
	int reads = 0;
	int writes = 0;
	
	//Creates an instance of ofstream, and opens example.txt
	ofstream a_file ( sa+"_result.txt" );

	// doing matrix product
	for(int i =0;i<row_a;i++){
		for(int j=0;j<column_b;j++){
			temp = 0;
			for (int k=0;k<column_a;k++){
				temp += a[i][k]*b[k][j];
				reads += 2;
				fpo += 2;
			}
			result[i][j] = temp;
			writes += 1;
		}
	}

	cout<<endl<<"result"<<endl;
	a_file<<"Array:"<<endl;
	for (int i=0;i<row_a;i++){
		for (int j=0;j<column_b;j++){
			cout<<result[i][j]<<" ";
			a_file<<result[i][j]<<" ";
		}
		cout<<endl;
		a_file<<endl;
	}

	cout<<endl<<"result"<<endl;
	cout<<"FPOperation:"<<fpo<<endl;
	cout<<"Reads:"<<reads<<endl;
	cout<<"Writes:"<<writes<<endl;

	// Outputs to example.txt through a_file
	// a_file<<"This text will now be inside of example.txt";
	a_file<<endl<<"FPOperation:"<<fpo<<endl;
	a_file<<"Reads:"<<reads<<endl;
	a_file<<"Writes:"<<writes<<endl;

	return 0;
}