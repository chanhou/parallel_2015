#include <iostream>
#include "IDALib.h"
#include "math.h"

using namespace std;


// read/write to known coordinates
double read(const int DIM, const int nPoints, const int whichPt, 
	const int whichDim, const double *coords){

	return coords[ whichPt + whichDim * nPoints ];

};
void  write(const int DIM, const int nPoints, const int whichPt, 
	const int whichDim, double *coords, const double val){

	// record by column
	// x1,x2,y1,y2
	coords[ whichPt + whichDim * nPoints ] = val;

};

// read/write to grid coordinates
double readGrid(const int DIM, const int nGridPoints, 
	const int whichGridPt, const int whichDim, const double *gridCoords){

	return gridCoords[ whichGridPt +  nGridPoints * whichDim ];

};
void  writeGrid(const int DIM, const int nGridPoints, 
	const int whichGridPt, const int whichDim, double *gridCoords, 
	const double val){


	// gridCoords[ whichGridPt * DIM + whichDim ] = val; 
	gridCoords[ whichGridPt +  nGridPoints * whichDim ] = val; // x1,x2,...,y1,y2,...

};

// read/write to known values
double readAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, const double *values){


	return values[ whichPt + whichAttr * nPoints];

};
void   writeAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, double *values, 
	const double val){

	// v1_1,v2_1,v3_1,....v1_2,v2_2,v3_2,.....
	values[ whichPt + whichAttr * nPoints] = val;

};

// read/write to grid values
double readGridAttribute(const int noAttr, const int nPoints, 
	const int whichPt, 
	const int whichAttr, const double *values){

	return values[ whichPt + whichAttr * nPoints ];

};
void   writeGridAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, double *values, 
	const double val){

	// g1_1,g2_1,g3_1,...,g1_2,g2_2,g3_2
	values[ whichPt + whichAttr * nPoints ] = val;

};

void computeBounds(const int DIM, const int nPoints, 
	const double *knownCoords, double(*bounds)[2]){

	for(int d=0; d<DIM;d++){
		bounds[d][0] = INFINITY; // min
		bounds[d][1] = -INFINITY; // max
		// cout<<bounds[d][1]<<endl;
	}

	for(int i=0; i< nPoints ; ++i){

		for(int d=0; d<DIM;d++){

			if ( knownCoords[ i + d * nPoints ] < bounds[d][0]){
				bounds[d][0] = knownCoords[ i + d * nPoints ]; // min
			}
			else if ( knownCoords[ i + d * nPoints ] > bounds[d][1] ){
				bounds[d][1] = knownCoords[ i + d * nPoints]; // max	
			}
			
		}
	}

};


void computeDistances(const int DIM, const int nPoints, 
	const double *knownCoords, const int nGrids, 
	const double *gridCoords, double *distances){

	// d(g1,d1), d(g2,d1),...,d(g1,d2),d(g2,d2),...

	for(int i=0;i<nGrids;++i){
		for(int j=0;j< nPoints; ++j){

			distances[ i + j * nGrids ] = 0;



			for(int d=0; d< DIM; d++){
				// change known coords and grid coords
				distances[ i + j * nGrids ] += pow(( knownCoords[ j + d * nPoints ] 
					- gridCoords[ i + nGrids * d ]) , 2);	
			}

			distances[ i + j * nGrids ] = pow(distances[ i + j * nGrids ], 0.5);
		}
	}

};



void computeWeights(const int nGrids, const int nPoints, 
	double *distances, double *weightSum, const double p ){

	// g1 = d(g1,d1)+d(g1,d2)+....

	for( int i=0; i< nGrids; i++){
		weightSum[i] = 0;
		for(int j=0;j< nPoints; j++){

			distances[ i + j * nGrids ] = pow(1/distances[ i + j * nGrids ], p);
			weightSum[i] += distances[ i + j * nGrids ];
			// weightSum[i] += pow(1/distances[ i*nPoints + j ], p);
		}
	}


};


void computeInterpolation(const int nValues, const int nGrids, 
	const int nPoints, const double *distances, 
	const double *weightSum, const double *knownValues, 
	double *gridValues){

	double sum_w_attr = 0;
	for(int i=0; i<nGrids; ++i){
		bool distan_zero = false;
		for(int v=0; v<nValues; ++v){
			// g1_1,g2_1,g3_1,..., g1_2,g2_2,g3_2,...
			
			gridValues[ i + v * nGrids  ] = 0;

			for(int p=0; p<nPoints; ++p){

				gridValues[ i + v * nGrids ] +=  distances[ i + p * nGrids ] * knownValues[ p + v * nPoints];

				if (isinf(distances[ i*nPoints + p ]) ){
					gridValues[ i + v * nGrids ] = knownValues[ p + v * nPoints ];
					distan_zero = true;
					break;
				}

			}
			
			if(!distan_zero) {
				gridValues[ i + v * nGrids ] /= weightSum[i];
			}

		}
	}

};
