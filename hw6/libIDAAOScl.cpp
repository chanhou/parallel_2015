#include <iostream>
#include "math.h"
#include "IDALib.h"

using namespace std;

typedef float MY_DATA_TYPE;


// read/write to known coordinates
MY_DATA_TYPE read(const int DIM, const int nPoints, const int whichPt, 
	const int whichDim, const MY_DATA_TYPE *coords){

	return coords[ whichPt * DIM + whichDim ];
};


void  write(const int DIM, const int nPoints, const int whichPt, 
	const int whichDim, MY_DATA_TYPE *coords, const MY_DATA_TYPE val){
	// record raw data by row
	// x1,y1,x2,y2
	coords[ whichPt * DIM + whichDim ] = val; 
};

// read/write to grid coordinates
MY_DATA_TYPE readGrid(const int DIM, const int nGridPoints, 
	const int whichGridPt, const int whichDim, const MY_DATA_TYPE *gridCoords){

	// return gridCoords[ whichGridPt * DIM + whichDim ];
	return gridCoords[ whichGridPt +  nGridPoints * whichDim ];

};

void  writeGrid(const int DIM, const int nGridPoints, 
	const int whichGridPt, const int whichDim, MY_DATA_TYPE *gridCoords, 
	const MY_DATA_TYPE val){

	// gridCoords[ whichGridPt * DIM + whichDim ] = val; // x1,y1,x2,y2
	gridCoords[ whichGridPt +  nGridPoints * whichDim ] = val; // x1,x2,...,y1,y2,...

};

// read/write to known values
MY_DATA_TYPE readAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, const MY_DATA_TYPE *values){

	return values[whichPt * noAttr + whichAttr];
};

void writeAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, MY_DATA_TYPE *values, 
	const MY_DATA_TYPE val){
	// v1_1,v1_2,v1_3,v2_1,v2_2,v2_3
	values[ whichPt * noAttr + whichAttr ] = val;

};

// read/write to grid values
MY_DATA_TYPE readGridAttribute(const int noAttr, const int nPoints, 
	const int whichPt, 
	const int whichAttr, const MY_DATA_TYPE *values){

	// g1_1,g1_2,g1_3,g2_1,g2_2,g2_3
	return values[ whichPt* noAttr + whichAttr ];
	// return values[ whichGridPt * DIM + whichDim ];

};


void   writeGridAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, MY_DATA_TYPE *values, 
	const MY_DATA_TYPE val){

	// g1_1,g1_2,g1_3,g2_1,g2_2,g2_3
	values[ whichPt* noAttr + whichAttr  ] = val;

};

void computeBounds(const int DIM, const int nPoints, 
	const MY_DATA_TYPE *knownCoords, MY_DATA_TYPE(*bounds)[2]){

	for(int d=0; d<DIM;d++){
		bounds[d][0] = INFINITY; // min
		bounds[d][1] = -INFINITY; // max
		// cout<<bounds[d][1]<<endl;
	}

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
	const MY_DATA_TYPE *knownCoords, const int nGrids, 
	const MY_DATA_TYPE *gridCoords, MY_DATA_TYPE *distances){ // distance (data1,grid1,...)
	// (g1,d1), (g1,d2),...,(g2,d1),(g2,d2),...
	for(int i=0;i<nGrids;++i){
		for(int j=0;j< nPoints; ++j){
			distances[ i*nPoints + j ] = 0;
			for(int d=0; d< DIM; d++){
				distances[ i*nPoints + j ] += pow(( knownCoords[ j*DIM + d ] 
					- gridCoords[ i + nGrids * d ]) , 2);	
				// gridCoords[ whichGridPt +  nGridPoints * whichDim ]
			}
			distances[ i*nPoints + j ] = pow(distances[ i*nPoints + j ], 0.5);
		}
	}
};

/*
	openCL version
*/
void computeDistances(cl::CommandQueue & cmd, cl::Program &prog, const int DIM, 
	const int nPoints, cl::Buffer &knownCoords, const int nGrids, 
	cl::Buffer &gridCoords, cl::Buffer &distances){

	// cout<<"great8.1"<<endl;

	// cl::Kernel kernel(*gpu.program, "count1");
	cl::Kernel kernel( prog, "computeDistances");

	// void computeDistances(int DIM, int nPoints,
	// 	__global float *knownCoords, int nGrids,
	// 	__global float *gridCoords, __global float* distances){

	kernel.setArg(0, (cl_int) DIM);
	kernel.setArg(1, (cl_int) nPoints);
	kernel.setArg(2, knownCoords);
	kernel.setArg(3, (cl_int) nGrids);
	kernel.setArg(4, gridCoords);
	kernel.setArg(5, distances);

	// cout<<"great8.2"<<endl;

	size_t lSize=16;
	cl::NDRange local( lSize, lSize );
	cl::NDRange global( ((cl_int)nGrids + lSize-1) / 
		lSize * lSize, ((cl_int)nPoints + lSize-1) / lSize * lSize );

	cl::Event event;
	MY_DATA_TYPE t = 0.0;
	// cl_uint count = 0;

	// cout<<"great8.3"<<endl;

	cmd.enqueueNDRangeKernel(kernel, cl::NullRange, global, local, 0, &event);

	// cout<<"great8.4"<<endl;

	event.wait();

	// cout<<"great8.5"<<endl;
	
	t += (event.getProfilingInfo<CL_PROFILING_COMMAND_END>() 
		- event.getProfilingInfo<CL_PROFILING_COMMAND_START>()) / 1.0e9;

	cout<<"cl t2: "<<t<<endl;
};


void computeWeights(const int nGrids, const int nPoints, 
	MY_DATA_TYPE *distances, MY_DATA_TYPE *weightSum, const MY_DATA_TYPE p ){

	// g1 = d(g1,d1)+d(g1,d2)+....
	for( int i=0; i< nGrids; i++){
		weightSum[i] = 0;
		for(int j=0;j< nPoints; j++){
			distances[ i*nPoints + j ] = pow(1/distances[ i*nPoints + j ], p);
			weightSum[i] += distances[ i*nPoints + j ];
			// weightSum[i] += pow(1/distances[ i*nPoints + j ], p);

		}
	}
};

/*
	openCL version
*/
void computeWeights(cl::CommandQueue &cmd, cl::Program &prog, const int nGrids, 
	const int nPoints, cl::Buffer &distances, cl::Buffer &weightSum, 
	const MY_DATA_TYPE p){

	cl::Kernel kernel( prog, "computeWeights");

	// void computeWeights(const int nGrids, const int nPoints, 
	// 	float *distances, float *weightSum, const float p ){

	kernel.setArg(0, (cl_int) nGrids);
	kernel.setArg(1, (cl_int) nPoints);
	kernel.setArg(2, distances);
	kernel.setArg(3, weightSum);
	kernel.setArg(4, (cl_float)p);

	// cout<<"great8.2"<<endl;

	size_t lSize= 256 ;
	cl::NDRange local( lSize );
	cl::NDRange global( ( (cl_int)nGrids + lSize-1) / lSize * lSize );

	cl::Event event;
	MY_DATA_TYPE t = 0.0;

	// cout<<"great8.3"<<endl;

	cmd.enqueueNDRangeKernel(kernel, cl::NullRange, global, local, 0, &event);

	// cout<<"great8.4"<<endl;

	event.wait();

	// cout<<"great8.5"<<endl;
	
	t += (event.getProfilingInfo<CL_PROFILING_COMMAND_END>() 
		- event.getProfilingInfo<CL_PROFILING_COMMAND_START>()) / 1.0e9;

	cout<<"cl t3: "<<t<<endl;	


};

void computeInterpolation(const int nValues, const int nGrids, 
	const int nPoints, const MY_DATA_TYPE *distances, 
	const MY_DATA_TYPE *weightSum, const MY_DATA_TYPE *knownValues, 
	MY_DATA_TYPE *gridValues){

	// MY_DATA_TYPE *knownCoords = new MY_DATA_TYPE[DIM * nPoints];
	// MY_DATA_TYPE *knownValues = new MY_DATA_TYPE[nPoints * nValues];
	// MY_DATA_TYPE (*bounds)[2] = new MY_DATA_TYPE[DIM][2];					// size of DIMS * 2
	// const int nGrids = (int)pow(nDivisions, DIM);				// # of grid points
	// MY_DATA_TYPE *gridCoords = new MY_DATA_TYPE[(size_t) pow(nDivisions, DIM) * DIM]; // multiply DIM since the dimension coordinate
	// MY_DATA_TYPE *distances = new MY_DATA_TYPE[nGrids * nPoints];
	// MY_DATA_TYPE *weightSum = new MY_DATA_TYPE[nGrids]; // why is grid? not nPoints? 
	// MY_DATA_TYPE *gridValues = new MY_DATA_TYPE[nGrids * nValues];


	MY_DATA_TYPE sum_w_attr = 0;
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

/*
	openCL version
*/

void computeInterpolation(cl::CommandQueue &cmd, cl::Program &prog, 
	const int nValues, const int nGrids, const int nPoints, 
	cl::Buffer &distances, cl::Buffer &weightSum, cl::Buffer &knownValues, 
	cl::Buffer &gridValues){

	// void computeInterpolation(const int nValues, const int nGrids, 
	// const int nPoints, const __global float *distances, 
	// const __global float *weightSum, const __global float *knownValues, 
	// __global float *gridValues){

	cl::Kernel kernel( prog, "computeInterpolation");

	kernel.setArg(0, (cl_int) nValues);
	kernel.setArg(1, (cl_int) nGrids);
	kernel.setArg(2, (cl_int) nPoints);
	kernel.setArg(3, distances);
	kernel.setArg(4, weightSum);
	kernel.setArg(5, knownValues);
	kernel.setArg(6, gridValues);

	// cout<<"great8.2"<<endl;

	size_t lSize=16;
	cl::NDRange local( lSize, lSize );
	cl::NDRange global( ((cl_int)nGrids + lSize-1) / 
		lSize * lSize, ((cl_int)nValues + lSize-1) / lSize * lSize );

	cl::Event event;
	MY_DATA_TYPE t = 0.0;
	// cl_uint count = 0;

	// cout<<"great8.3"<<endl;

	cmd.enqueueNDRangeKernel(kernel, cl::NullRange, global, local, 0, &event);

	// cout<<"great8.4"<<endl;

	event.wait();

	// cout<<"great8.5"<<endl;
	
	t += (event.getProfilingInfo<CL_PROFILING_COMMAND_END>() 
		- event.getProfilingInfo<CL_PROFILING_COMMAND_START>()) / 1.0e9;

	cout<<"cl t4: "<<t<<endl;



};