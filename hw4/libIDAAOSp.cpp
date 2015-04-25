#include <iostream>
#include "IDALib.h"
#include "math.h"

#include <omp.h>

using namespace std;


// read/write to known coordinates
double read(const int DIM, const int nPoints, const int whichPt, 
	const int whichDim, const double *coords){

	return coords[ whichPt * DIM + whichDim ];
};


void  write(const int DIM, const int nPoints, const int whichPt, 
	const int whichDim, double *coords, const double val){
	// record raw data by row
	// x1,y1,x2,y2
	coords[ whichPt * DIM + whichDim ] = val; 
};

// read/write to grid coordinates
double readGrid(const int DIM, const int nGridPoints, 
	const int whichGridPt, const int whichDim, const double *gridCoords){

	// return gridCoords[ whichGridPt * DIM + whichDim ];
	return gridCoords[ whichGridPt +  nGridPoints * whichDim ];

};

void  writeGrid(const int DIM, const int nGridPoints, 
	const int whichGridPt, const int whichDim, double *gridCoords, 
	const double val){

	// gridCoords[ whichGridPt * DIM + whichDim ] = val; // x1,y1,x2,y2
	gridCoords[ whichGridPt +  nGridPoints * whichDim ] = val; // x1,x2,...,y1,y2,...

};

// read/write to known values
double readAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, const double *values){

	return values[whichPt * noAttr + whichAttr];
};

void writeAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, double *values, 
	const double val){
	// v1_1,v1_2,v1_3,v2_1,v2_2,v2_3
	values[ whichPt * noAttr + whichAttr ] = val;

};

// read/write to grid values
double readGridAttribute(const int noAttr, const int nPoints, 
	const int whichPt, 
	const int whichAttr, const double *values){

	// g1_1,g1_2,g1_3,g2_1,g2_2,g2_3
	return values[ whichPt* noAttr + whichAttr ];
	// return values[ whichGridPt * DIM + whichDim ];

};


void writeGridAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, double *values, 
	const double val){

	// g1_1,g1_2,g1_3,g2_1,g2_2,g2_3
	values[ whichPt* noAttr + whichAttr  ] = val;

};

void computeBounds(const int DIM, const int nPoints, 
	const double *knownCoords, double(*bounds)[2]){

	/*
	para
	*/
	// operation too small, no need do parallel
	for(int d=0; d<DIM;d++){
		bounds[d][0] = INFINITY; // min
		bounds[d][1] = -INFINITY; // max
		// cout<<bounds[d][1]<<endl;
	}

	/*
	para
	*/
	#pragma omp parallel for collapse(2)
	for(int i=0; i< nPoints ; ++i){ 

		for(int d=0; d<DIM;d++){

			if ( knownCoords[ i*DIM + d ] < bounds[d][0]){
				bounds[d][0] = knownCoords[ i*DIM + d ]; // min
			}
			else if ( knownCoords[ i*DIM + d ] > bounds[d][1] ){
				bounds[d][1] = knownCoords[ i*DIM + d ]; // max	
			}
			
		}
	}



};


void computeDistances(const int DIM, const int nPoints, 
	const double *knownCoords, const int nGrids, 
	const double *gridCoords, double *distances){ // distance (data1,grid1,...)
	
	/*
	para
	*/
	// (g1,d1), (g1,d2),...,(g2,d1),(g2,d2),...
	// cant use collapse since is not perfect retangular
	double temp;

	#pragma omp parallel for
	for(int i=0;i<nGrids;++i){
		// cout<<"Hello from thread, nthreads \n"<<omp_get_thread_num()<<", "<<omp_get_num_threads(); 
		for(int j=0;j< nPoints; ++j){

			distances[ i*nPoints + j ] = 0;
			// temp = 0;
			// #pragma omp parallel for reduction(+:temp)
			for(int d=0; d< DIM; d++){

				distances[ i*nPoints + j ] += pow(( knownCoords[ j*DIM + d ] 
					- gridCoords[ i + nGrids * d ]) , 2);
				// temp += pow(( knownCoords[ j*DIM + d ] 
				// 	- gridCoords[ i + nGrids * d ]) , 2);

			}

			distances[ i*nPoints + j ] = pow(distances[ i*nPoints + j ], 0.5);
			// distances[ i*nPoints + j ] = pow( temp , 0.5);
		}
	}
};


void computeWeights(const int nGrids, const int nPoints, 
	double *distances, double *weightSum, const double p ){

	/*
	para
	*/
	// g1 = d(g1,d1)+d(g1,d2)+....
	double temp;
	for( int i=0; i< nGrids; i++){
		// weightSum[i] = 0;
		temp = 0;
		#pragma omp parallel for reduction(+:temp)
		for(int j=0;j< nPoints; j++){
			distances[ i*nPoints + j ] = pow(1/distances[ i*nPoints + j ], p);
			// weightSum[i] += distances[ i*nPoints + j ];
			temp += distances[ i*nPoints + j ];
		}
		weightSum[i] = temp;
	}
};


void computeInterpolation(const int nValues, const int nGrids, 
	const int nPoints, const double *distances, 
	const double *weightSum, const double *knownValues, 
	double *gridValues){

	// double *knownCoords = new double[DIM * nPoints];
	// double *knownValues = new double[nPoints * nValues];
	// double (*bounds)[2] = new double[DIM][2];					// size of DIMS * 2
	// const int nGrids = (int)pow(nDivisions, DIM);				// # of grid points
	// double *gridCoords = new double[(size_t) pow(nDivisions, DIM) * DIM]; // multiply DIM since the dimension coordinate
	// double *distances = new double[nGrids * nPoints];
	// double *weightSum = new double[nGrids]; // why is grid? not nPoints? 
	// double *gridValues = new double[nGrids * nValues];

	double sum_w_attr = 0;

	#pragma omp parallel for
	for(int i=0; i<nGrids; ++i){
		bool distan_zero = false;
		for(int v=0; v<nValues; ++v){
			// g1_1,g1_2,g1_3,...,g2_1,g2_2,g2_3
			gridValues[ i*nValues + v ] = 0;

			for(int p=0; p<nPoints; ++p){
				gridValues[ i*nValues + v ] +=  distances[ i*nPoints + p ] * knownValues[ p*nValues + v];
				
				if (isinf(distances[ i*nPoints + p ])){ //isinf
					gridValues[ i*nValues + v ] = knownValues[ p*nValues + v];
					distan_zero = true;
					break;
				}
			}

			if(!distan_zero) {
				gridValues[ i*nValues + v ] /= weightSum[i];
			}
		}
	}


};

// void computeInterpolation(const int nValues, const int nGrids, 
// 	const int nPoints, const double *distances, 
// 	const double *weightSum, const double *knownValues, 
// 	double *gridValues){

// 	// double *knownCoords = new double[DIM * nPoints];
// 	// double *knownValues = new double[nPoints * nValues];
// 	// double (*bounds)[2] = new double[DIM][2];					// size of DIMS * 2
// 	// const int nGrids = (int)pow(nDivisions, DIM);				// # of grid points
// 	// double *gridCoords = new double[(size_t) pow(nDivisions, DIM) * DIM]; // multiply DIM since the dimension coordinate
// 	// double *distances = new double[nGrids * nPoints];
// 	// double *weightSum = new double[nGrids]; // why is grid? not nPoints? 
// 	// double *gridValues = new double[nGrids * nValues];

// 	double temp = 0;
// 	/*
// 	para
// 	*/

// 	#pragma omp parallel for 
// 	for(int i=0; i<nGrids; ++i){

// 		for(int j=0; j<nValues; ++j){
// 			// g1_1,g1_2,g1_3,...,g2_1,g2_2,g2_3

// 			// gridValues[ i*nValues + j ] = 0;
// 			temp = 0;
// 			bool distan_zero = false;
// 			int record_p = 0;

// 			#pragma omp parallel for reduction(+: temp)
// 			for(int p=0; p<nPoints; ++p){
// 				// gridValues[ i*nValues + j ] +=  distances[ i*nPoints + p ] * knownValues[ p*nValues + j];
// 				temp +=  distances[ i*nPoints + p ] * knownValues[ p*nValues + j];
				
// 				if (isinf(distances[ i*nPoints + p ])){ //isinf
// 					// temp = knownValues[ p*nValues + j];
// 					distan_zero = true;
// 					record_p = p;
// 					// break;
// 				}
// 			}

// 			if(!distan_zero) {
// 				gridValues[ i*nValues + j ] = temp/weightSum[i];
// 				// gridValues[ i*nValues + j ] /= weightSum[i];
// 			}
// 			else{ // distan_zero is true
// 				gridValues[ i*nValues + j ] = knownValues[ record_p * nValues + j];
// 			}
// 		}
// 	}


// };