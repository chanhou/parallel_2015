#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "IDALib.h"

#include "stopWatch.h"

#include <omp.h> 

using namespace std;

static void g1(const int NPOINTS, const int DIM, const double(*bounds)[2], const int nDiv,
 double *grids, const int whichDIM, double *X, int &where) {
	if (whichDIM == DIM) {
		// store the grid point & return
		for (int dim = 0; dim < DIM; ++dim) writeGrid(DIM, NPOINTS, where, dim, grids, X[dim]);
		// for (int dim = 0; dim < DIM; ++dim) cout << "*" << X[dim] << "* "; cout << endl;
		++where;
		return;
	}
	const double inc = (bounds[whichDIM][1] - bounds[whichDIM][0]) / (nDiv - 1);
	X[whichDIM] = bounds[whichDIM][0];
	for (int i = 0; i < nDiv; ++i) {
		g1(NPOINTS, DIM, bounds, nDiv, grids, whichDIM + 1, X, where);
		X[whichDIM] += inc;
	}
}

// x1, x2, x3, ... y1, y2, y3, ... z1, z1, z3, ...
void computeGridCoordinates(const int DIM, const double(*bounds)[2], const int nDiv, double *grids) {
	const int NPOINTS = pow(nDiv, DIM);
	double *X = new double[DIM]; // x,y,z,...
	int where = 0;
	g1(NPOINTS, DIM, bounds, nDiv, grids, 0, X, where);
	delete[]X;
}

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
	double *knownCoords = new double[DIM * nPoints];
	double *knownValues = new double[nPoints * nValues];
	double (*bounds)[2] = new double[DIM][2];					// size of DIMS * 2
	const int nGrids = (int)pow(nDivisions, DIM);				// # of grid points
	double *gridCoords = new double[(size_t) pow(nDivisions, DIM) * DIM]; // multiply DIM since the dimension coordinate
	double *distances = new double[nGrids * nPoints];
	double *weightSum = new double[nGrids]; // why is grid? not nPoints?
	double *gridValues = new double[nGrids * nValues];

	double t1;
	double t2;
	double t3;
	double t4;
	stopWatch timer;

	// read data from the specified file and store data into appropriate data structures using write & writeAttribute functions
	for (int pt = 0; pt < nPoints; ++pt) {
		for (int dim = 0; dim < DIM; ++dim) {
			double tmp;
			inp >> tmp;
			write(DIM, nPoints, pt, dim, knownCoords, tmp);
		}
		for (int attr = 0; attr < nValues; ++attr) {
			double tmp;
			inp >> tmp;
			writeAttribute(nValues, nPoints, pt, attr, knownValues, tmp);
		}
	}
	inp.close();
	// show data
	// for (int i = 0; i < nPoints * DIM; ++i) cout << knownCoords[i] << " "; cout << endl;
	// for (int i = 0; i < nPoints * nValues; ++i) cout << knownValues[i] << " "; cout << endl;

	// #pragma omp parallel 
	// cout<<"Hello from thread, nthreads \n"<<omp_get_thread_num()<<", "<<omp_get_num_threads(); 
	// cout<<endl<<omp_get_max_threads()<<endl;

	// while(true)	{
	// 	#pragma omp parallel for
	// 	for(int i=0; i<10; ++i){
	// 		cout<<i<<endl;
	// 	}
	// 	// printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads()); 
	// 	cout<<endl<<omp_get_max_threads();
	// }

	// int where = 0;


	// find bounds of known data points
	timer.start();
	computeBounds(DIM, nPoints, knownCoords, bounds);
	timer.stop();
	t1 = timer.elapsedTime();
	cout<<DIM<<" "<<nPoints<<" "<<nValues<<" "<<nDivisions<<" t1 "<<t1<<endl;
	// for (int i = 0; i < DIM; ++i) cout << bounds[i][0] << ":" << bounds[i][1] << endl; cout << endl;

	// // create grid points in unknownCoords
	computeGridCoordinates(DIM, bounds, nDivisions, gridCoords);
	// for (int i = 0; i < nGrids*DIM; ++i) cout << gridCoords[i] << " "; cout << endl<<endl;

	// // step 3. compute distances between all grids points and all known data points
	timer.start();
	computeDistances(DIM, nPoints, knownCoords, nGrids, gridCoords, distances);
	timer.stop();
	t2 = timer.elapsedTime();
	// cout<<"DIM "<<DIM<<" t2 "<<t2<<endl;
	cout<<DIM<<" "<<nPoints<<" "<<nValues<<" "<<nDivisions<<" t2 "<<t2<<endl;
	// for (int i = 0; i < nPoints*nGrids; ++i) cout << distances[i] << " "; cout << endl;

	// // step 4. turn the distance array into weight array, and compute the total weight
	timer.start();
	computeWeights(nGrids, nPoints, distances, weightSum);
	timer.stop();
	t3 = timer.elapsedTime();
	// cout<<"DIM "<<DIM<<" t3 "<<t3<<endl;
	cout<<DIM<<" "<<nPoints<<" "<<nValues<<" "<<nDivisions<<" t3 "<<t3<<endl;
	// for (int i = 0; i < nPoints*nGrids; ++i) cout << distances[i] << " "; cout << endl;

	// // step 5 & 6. compute sum (weights * known Values) / totalWeight
	timer.start();
	computeInterpolation(nValues, nGrids, nPoints, distances, 
		weightSum, knownValues, gridValues);
	timer.stop();
	t4 = timer.elapsedTime();
	// cout<<"DIM "<<DIM<<" t4 "<<t4<<endl;
	cout<<DIM<<" "<<nPoints<<" "<<nValues<<" "<<nDivisions<<" t4 "<<t4<<endl;
	// All calculations are finished.  Write grid data to the specified output file.
	

	/*  for testting speed  */
	
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

