#ifndef IDALIB
#define IDALIB

#define __CL_ENABLE_EXCEPTIONS
#include "CL/cl.hpp"

typedef float MY_DATA_TYPE;

// read/write to known coordinates
MY_DATA_TYPE read(const int DIM, const int nPoints, 
	const int whichPt, const int whichDim, 
	const MY_DATA_TYPE *coords);
void  write(const int DIM, const int nPoints, const int whichPt, 
	const int whichDim, MY_DATA_TYPE *coords, const MY_DATA_TYPE val);

// read/write to grid coordinates
MY_DATA_TYPE readGrid(const int DIM, const int nGridPoints, 
	const int whichGridPt, const int whichDim, 
	const MY_DATA_TYPE *gridCoords);
void  writeGrid(const int DIM, const int nGridPoints, 
	const int whichGridPt, 
	const int whichDim, MY_DATA_TYPE *gridCoords, 
	const MY_DATA_TYPE val);

// read/write to known values
MY_DATA_TYPE readAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, 
	const MY_DATA_TYPE *values);
void   writeAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, MY_DATA_TYPE *values, 
	const MY_DATA_TYPE val);

// read/write to grid values
MY_DATA_TYPE readGridAttribute(const int noAttr, const int nPoints,
 const int whichPt, const int whichAttr, const MY_DATA_TYPE *values);
void   writeGridAttribute(const int noAttr, const int nPoints, 
	const int whichPt, const int whichAttr, MY_DATA_TYPE *values, 
	const MY_DATA_TYPE val);

void computeBounds(const int DIM, const int nPoints, 
	const MY_DATA_TYPE *knownCoords, MY_DATA_TYPE(*bounds)[2]);
void computeDistances(const int DIM, const int nPoints, 
	const MY_DATA_TYPE *knownCoords, const int nGrids, 
	const MY_DATA_TYPE *gridCoords, MY_DATA_TYPE *distances);
void computeWeights(const int nGrids, const int nPoints, 
	MY_DATA_TYPE *distances, MY_DATA_TYPE *weightSum, 
	const MY_DATA_TYPE p = 2.0);
void computeInterpolation(const int nValues, const int nGrids, 
	const int nPoints, const MY_DATA_TYPE *distances, 
	const MY_DATA_TYPE *weightSum, 
	const MY_DATA_TYPE *knownValues, MY_DATA_TYPE *gridValues);

// For OpenCL version
void computeDistances(cl::CommandQueue &, cl::Program &prog, 
	const int DIM, const int nPoints, cl::Buffer &knownCoords, 
	const int nGrids, cl::Buffer &gridCoords, 
	cl::Buffer &distances);

void computeWeights(cl::CommandQueue &, cl::Program &prog, 
	const int nGrids, const int nPoints, cl::Buffer &distances, 
	cl::Buffer &weightSum, const MY_DATA_TYPE p = 2.0);

void computeInterpolation(cl::CommandQueue &, cl::Program &prog, 
	const int nValues, const int nGrids, const int nPoints, 
	cl::Buffer &distances, cl::Buffer &weightSum, 
	cl::Buffer &knownValues, cl::Buffer &gridValues);

#endif
