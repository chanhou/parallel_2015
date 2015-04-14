#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "IDALib.h"
using namespace std;

static void g1(const int NPOINTS, const int DIM, const double(*bounds)[2], const int nDiv, double *grids, const int whichDIM, double *X, int &where) {
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
	double *X = new double[DIM];
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
	double *gridCoords = new double[(size_t) pow(nDivisions, DIM) * DIM];
	double *distances = new double[nGrids * nPoints];
	double *weightSum = new double[nGrids];
	double *gridValues = new double[nGrids * nValues];

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

	// int where = 0;
	// find bounds of known data points
	computeBounds(DIM, nPoints, knownCoords, bounds);
	// for (int i = 0; i < DIM; ++i) cout << bounds[i][0] << ":" << bounds[i][1] << endl; cout << endl;

	// create grid points in unknownCoords
	computeGridCoordinates(DIM, bounds, nDivisions, gridCoords);
	// for (int i = 0; i < nGrids*DIM; ++i) cout << gridCoords[i] << " "; cout << endl;

	// step 3. compute distances between all grids points and all known data points
	computeDistances(DIM, nPoints, knownCoords, nGrids, gridCoords, distances);
	// for (int i = 0; i < nPoints*nGrids; ++i) cout << distances[i] << " "; cout << endl;

	// step 4. turn the distance array into weight array, and compute the total weight
	computeWeights(nGrids, nPoints, distances, weightSum);
	// for (int i = 0; i < nPoints*nGrids; ++i) cout << distances[i] << " "; cout << endl;

	// step 5 & 6. compute sum (weights * known Values) / totalWeight
	computeInterpolation(nValues, nGrids, nPoints, distances, weightSum, knownValues, gridValues);

	// All calculations are finished.  Write grid data to the specified output file.
	ofstream outp(argv[2]);
	// Output points
	if (outp.good()) {
		/*
		// This part writes scattered data points
		for (int i = 0; i < nPoints; ++i) {
			for (int dim = 0; dim < DIM; ++dim) {
				outp << read(DIM, nPoints, i, dim, knownCoords) << " ";
			}
			for (int attr = 0; attr < nValues; ++attr) {
				outp << readAttribute(nValues, nPoints, i, attr, knownValues);
			}
			outp << knownValues[i] << "\n";
		}
		outp << "\n\n";
		*/
		
		for (int i = 0; i < nGrids; ++i) {
			for (int dim = 0; dim < DIM; ++dim) {
				outp << readGrid(DIM, nGrids, i, dim, gridCoords) << " ";
			}
			for (int attr = 0; attr < nValues; ++attr) {
				outp << readGridAttribute(nValues, nGrids, i, attr, gridValues) << " ";
			}
			outp << "\n";
		}
		
		outp.close();
	}
	
	delete[] bounds;
	delete[] knownCoords;
	delete[] gridCoords;
	delete[] weightSum;
	delete[] distances;
	delete[] knownValues;
	delete[] gridValues;

	return 0;
}
