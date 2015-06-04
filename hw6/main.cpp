#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "IDALib.h"

// #include "stopWatch.h"

#include <omp.h> 

#include "YoUtil.hpp"

#define DATA_TYPE float
#define T_MIN 1.0
typedef unsigned int uint;

using namespace std;

static void g1(const int NPOINTS, const int DIM, const DATA_TYPE(*bounds)[2], const int nDiv,
 DATA_TYPE *grids, const int whichDIM, DATA_TYPE *X, int &where) {
	if (whichDIM == DIM) {
		// store the grid point & return
		for (int dim = 0; dim < DIM; ++dim) writeGrid(DIM, NPOINTS, where, dim, grids, X[dim]);
		// for (int dim = 0; dim < DIM; ++dim) cout << "*" << X[dim] << "* "; cout << endl;
		++where;
		return;
	}
	const DATA_TYPE inc = (bounds[whichDIM][1] - bounds[whichDIM][0]) / (nDiv - 1);
	X[whichDIM] = bounds[whichDIM][0];
	for (int i = 0; i < nDiv; ++i) {
		g1(NPOINTS, DIM, bounds, nDiv, grids, whichDIM + 1, X, where);
		X[whichDIM] += inc;
	}
}

// x1, x2, x3, ... y1, y2, y3, ... z1, z1, z3, ...
void computeGridCoordinates(const int DIM, const DATA_TYPE(*bounds)[2], const int nDiv, DATA_TYPE *grids) {
	const int NPOINTS = pow(nDiv, DIM);
	DATA_TYPE *X = new DATA_TYPE[DIM]; // x,y,z,...
	int where = 0;
	g1(NPOINTS, DIM, bounds, nDiv, grids, 0, X, where);
	delete[]X;
}

/*
	GPU version define
*/
struct GPU {
	cl::CommandQueue cmdQueue;
	cl::Program *program;
	cl::Context *context;
	GPU() {
		program = nullptr;
		context = nullptr;
	}
	~GPU() {
		delete program;
		delete context;
	}
	cl::Context& operator() () {
		return *context;
	}
};



int main(int argc, char **argv) {
	if (argc < 3) {
		cerr << argv[0] << " [input filename] [output filename]" << endl;
		return 255;
	}
	ifstream inp(argv[1]);
	if (!inp.good()) {
		cerr << "\nError opening file: " << argv[1] << endl;
		return 254;
	}
	int DIM, nPoints, nValues, nDivisions;
	inp >> DIM >> nPoints >> nValues >> nDivisions;

	// Allocate memory for storing all necessary data
	DATA_TYPE *knownCoords = new DATA_TYPE[DIM * nPoints];
	DATA_TYPE *knownValues = new DATA_TYPE[nPoints * nValues];
	DATA_TYPE (*bounds)[2] = new DATA_TYPE[DIM][2];					// size of DIMS * 2
	const int nGrids = (int)pow(nDivisions, DIM);				// # of grid points
	DATA_TYPE *gridCoords = new DATA_TYPE[(size_t) pow(nDivisions, DIM) * DIM]; // multiply DIM since the dimension coordinate
	DATA_TYPE *distances = new DATA_TYPE[nGrids * nPoints];
	DATA_TYPE *weightSum = new DATA_TYPE[nGrids]; // why is grid? not nPoints?
	DATA_TYPE *gridValues = new DATA_TYPE[nGrids * nValues];

	DATA_TYPE t1;
	DATA_TYPE t2;
	DATA_TYPE t3;
	DATA_TYPE t4;
	stopWatch timer;

	// read data from the specified file and store data into appropriate data structures using write & writeAttribute functions
	for (int pt = 0; pt < nPoints; ++pt) {
		for (int dim = 0; dim < DIM; ++dim) {
			DATA_TYPE tmp;
			inp >> tmp;
			write(DIM, nPoints, pt, dim, knownCoords, tmp);
		}
		for (int attr = 0; attr < nValues; ++attr) {
			DATA_TYPE tmp;
			inp >> tmp;
			writeAttribute(nValues, nPoints, pt, attr, knownValues, tmp);
		}
	}
	inp.close();
	// cout<<"great1"<<endl;


	/*
	GPU version define
	*/
	GPU gpu;

	try {
		// 1. Get context and command queues
		vector<cl::CommandQueue> cmdQueues;
		gpu.context = getContext(CL_DEVICE_TYPE_GPU, cmdQueues);
		gpu.cmdQueue = cmdQueues[0];

		// 2. read source code from file.
		string source = readSourceCode("kernel.cl");
		if (source.length() == 0) {
			return 255;
		}

		// 3. Compile code
		gpu.program = compile(gpu(), source);

		// 4. copy data from host to device
		
		// find bounds of known data points
		// non opencl
		timer.start();
		computeBounds(DIM, nPoints, knownCoords, bounds);
		// cout<<"great3"<<endl;
		timer.stop();
		t1 = timer.elapsedTime();
		// cout<<DIM<<" "<<nPoints<<" "<<nValues<<" "<<nDivisions<<endl;
		// cout<<"t1 "<<t1<<endl;

		cout<<DIM<<" "<<nPoints<<" "<<nValues<<" "<<nDivisions<<" t1 "<<t1<<endl;

		// create grid points in unknownCoords
		// non opencl
		computeGridCoordinates(DIM, bounds, nDivisions, gridCoords);
		// for (int i = 0; i < nGrids*DIM; ++i) cout << gridCoords[i] << " "; cout << endl<<endl;

		// copy data from host to device
		cl::Buffer knownCoords_(gpu(), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
			DIM * nPoints * sizeof(cl_float), knownCoords);
		
		cl::Buffer knownValues_(gpu(), CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, 
			nPoints * nValues * sizeof(cl_float), knownValues);
		
		cl::Buffer gridCoords_(gpu(), CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, 
			pow(nDivisions, DIM) * DIM * sizeof(cl_float), gridCoords);

		cl::Buffer distances_(gpu(), CL_MEM_READ_WRITE , 
			nGrids * nPoints * sizeof(cl_float));
		// cl::Buffer distances_(gpu(), CL_MEM_READ_WRITE , 
		// 	nGrids * nPoints * sizeof(cl_float), distances);

		// step 3. compute distances between all grids points and all known data points
		// opencl version
		// computeDistances(DIM, nPoints, knownCoords_, nGrids, gridCoords_, distances_);
		timer.start();
		computeDistances( gpu.cmdQueue , *gpu.program, DIM, nPoints, 
			knownCoords_, nGrids, gridCoords_, distances_);
		timer.stop();
		t2 = timer.elapsedTime();
		// cout<<"t2 "<<t2<<endl;
		cout<<DIM<<" "<<nPoints<<" "<<nValues<<" "<<nDivisions<<" t2 "<<t2<<endl;
		// cout<<"great9"<<endl;

		// // move data from device to host
		// gpu.cmdQueue.enqueueReadBuffer(distances_, CL_TRUE, 0, 
		// 	nGrids * nPoints * sizeof(cl_float), distances);

		// cout<<"great10"<<endl;

		// for (int i = 0; i < nPoints*nGrids; ++i) cout << distances[i] << " "; cout << endl;

		// cl::Buffer weightSum_(gpu(), CL_MEM_READ_WRITE , 
		// 	nGrids * sizeof(cl_float));
		for (int i = 0; i < nGrids; ++i) weightSum[i] = 0 ;
		cl::Buffer weightSum_(gpu(), CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, 
			nGrids * sizeof(cl_float), weightSum);

		// step 4. turn the distance array into weight array, and compute the total weight
		timer.start();
		computeWeights( gpu.cmdQueue , *gpu.program, nGrids, nPoints, 
			distances_, weightSum_);
		timer.stop();
		t3 = timer.elapsedTime();
		// cout<<"t3 "<<t3<<endl;
		cout<<DIM<<" "<<nPoints<<" "<<nValues<<" "<<nDivisions<<" t3 "<<t3<<endl;

		
		// for (int i = 0; i < nPoints*nGrids; ++i) cout << distances[i] << " "; cout << endl;
		// for (int i = 0; i < nGrids; ++i) cout << weightSum[i] << " "; cout << endl;

		cl::Buffer gridValues_(gpu(), CL_MEM_READ_WRITE , 
			nGrids * nValues * sizeof(cl_float));
		// cl::Buffer gridValues_(gpu(), CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, 
		// 	nGrids * nValues * sizeof(cl_float), gridValues);

		// // step 5 & 6. compute sum (weights * known Values) / totalWeight
		timer.start();
		computeInterpolation(gpu.cmdQueue , *gpu.program, nValues, nGrids, 
			nPoints, distances_, 
			weightSum_, knownValues_, gridValues_);
		timer.stop();
		t4 = timer.elapsedTime();
		// cout<<"t4 "<<t4<<endl;
		cout<<DIM<<" "<<nPoints<<" "<<nValues<<" "<<nDivisions<<" t4 "<<t4<<endl;

		
		// move data from device to host
		gpu.cmdQueue.enqueueReadBuffer(distances_, CL_TRUE, 0, 
			nGrids * nPoints * sizeof(cl_float), distances);
		gpu.cmdQueue.enqueueReadBuffer(weightSum_, CL_TRUE, 0, 
			nGrids  * sizeof(cl_float), weightSum);
		gpu.cmdQueue.enqueueReadBuffer(gridValues_, CL_TRUE, 0, 
			nGrids * nValues  * sizeof(cl_float), gridValues);

		// for (int i = 0; i < nGrids * nValues ; ++i) cout << gridValues[i] << " "; cout << endl;
		
		// All calculations are finished.  Write grid data to the specified output file.
	
		// ofstream outp(argv[2]);
		// // Output points
		// if (outp.good()) {
			
		// 	// This part writes scattered data points
		// 	// for (int i = 0; i < nPoints; ++i) {
		// 	// 	for (int dim = 0; dim < DIM; ++dim) {
		// 	// 		outp << read(DIM, nPoints, i, dim, knownCoords) << " ";
		// 	// 	}
		// 	// 	for (int attr = 0; attr < nValues; ++attr) {
		// 	// 		outp << readAttribute(nValues, nPoints, i, attr, knownValues) << " ";
		// 	// 	}
		// 	// 	outp << knownValues[i] << "\n";
		// 	// }	
		// 	// outp << "\n\n";
			
		// 	for (int i = 0; i < nGrids; ++i) {
		// 		for (int dim = 0; dim < DIM; ++dim) {
		// 			outp << readGrid(DIM, nGrids, i, dim, gridCoords) << " ";
		// 		}
		// 		for (int attr = 0; attr < nValues; ++attr) {
		// 			outp << readGridAttribute(nValues, nGrids, i, attr, gridValues) << " ";
		// 		}
		// 		outp << "\n";
		// 	}
			
		// 	outp.close();
		// }

		delete[] bounds;
		delete[] knownCoords;
		delete[] gridCoords;
		delete[] weightSum;
		delete[] distances;
		delete[] knownValues;
		delete[] gridValues;

		return 0;
	} 
	catch (cl::Error &e) {
		cerr << "\nMain: " << e.what();
		cerr << "\nError no: " << e.err() << endl;
		return 255;
	}


}

