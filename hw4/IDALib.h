#ifndef IDALIB
#define IDALIB

// read/write to known coordinates
double read(const int DIM, const int nPoints, const int whichPt, const int whichDim, const double *coords);
void  write(const int DIM, const int nPoints, const int whichPt, const int whichDim, double *coords, const double val);

// read/write to grid coordinates
double readGrid(const int DIM, const int nGridPoints, const int whichGridPt, const int whichDim, const double *gridCoords);
void  writeGrid(const int DIM, const int nGridPoints, const int whichGridPt, const int whichDim, double *gridCoords, const double val);

// read/write to known values
double readAttribute(const int noAttr, const int nPoints, const int whichPt, const int whichAttr, const double *values);
void   writeAttribute(const int noAttr, const int nPoints, const int whichPt, const int whichAttr, double *values, const double val);

// read/write to grid values
double readGridAttribute(const int noAttr, const int nPoints, const int whichPt, const int whichAttr, const double *values);
void   writeGridAttribute(const int noAttr, const int nPoints, const int whichPt, const int whichAttr, double *values, const double val);

void computeBounds(const int DIM, const int nPoints, const double *knownCoords, double(*bounds)[2]);
void computeDistances(const int DIM, const int nPoints, const double *knownCoords, 
	const int nGrids, const double *gridCoords, double *distances);
void computeWeights(const int nGrids, const int nPoints, double *distances, 
	double *weightSum, const double p = 2.0);
void computeInterpolation(const int nValues, const int nGrids, 
	const int nPoints, const double *distances, 
	const double *weightSum, const double *knownValues, double *gridValues);


#endif
