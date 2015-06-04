

__kernel
void computeDistances_AOS(const int DIM,const int nPoints,
	const __global float *knownCoords,const int nGrids,
	const __global float *gridCoords, __global float* distances){

	int k;
	int i = get_global_id(0);
	int j = get_global_id(1);

	if( i >= nGrids || j >= nPoints ) return;

	float tmp = 0.0f;
	for (k = 0; k < DIM; k++){
		// tmp += A[i*Ndim+k] * B[k*Pdim+j];
		tmp += pow( ( knownCoords[ j*DIM + k ] - gridCoords[ i + nGrids * k ]) , 2.0f);
	}
	distances[ i*nPoints + j] = pow( tmp , 0.5f );

	// (g1,d1), (g1,d2),...,(g2,d1),(g2,d2),...
	// for(int i=0;i<nGrids;++i){
	// 	for(int j=0;j< nPoints; ++j){
	// 		distances[ i*nPoints + j ] = 0;
	// 		for(int d=0; d< DIM; d++){
	// 			distances[ i*nPoints + j ] += pow(( knownCoords[ j*DIM + d ] 
	// 				- gridCoords[ i + nGrids * d ]) , 2);	
	// 			// gridCoords[ whichGridPt +  nGridPoints * whichDim ]
	// 		}
	// 		distances[ i*nPoints + j ] = pow(distances[ i*nPoints + j ], 0.5);
	// 	}
	// }

}

__kernel
void computeWeights_AOS(const int nGrids, const int nPoints, 
	__global float *distances, __global float *weightSum, const float p ){

	int j;
	int i = get_global_id(0);

	if( i >= nGrids ) return;

	for(j=0; j<nPoints; j++){
		distances[ i*nPoints + j ] = pow( 1.0f/distances[ i*nPoints + j ], p );	
		weightSum[i] += distances[ i*nPoints + j ];
	}
	
	// g1 = d(g1,d1)+d(g1,d2)+....
	// for( int i=0; i< nGrids; i++){
	// 	weightSum[i] = 0;
	// 	for(int j=0;j< nPoints; j++){
	// 		distances[ i*nPoints + j ] = pow(1/distances[ i*nPoints + j ], p);
	// 		weightSum[i] += distances[ i*nPoints + j ];
	// 		// weightSum[i] += pow(1/distances[ i*nPoints + j ], p);

	// 	}
	// }

};

__kernel 
void computeInterpolation_AOS(const int nValues, const int nGrids, 
	const int nPoints, const __global float *distances, 
	const __global float *weightSum, const __global float *knownValues, 
	__global float *gridValues){


	int p;
	int i = get_global_id(0);
	int v = get_global_id(1);

	if( i < nGrids  ) {
		bool distan_zero = false;
		if( v < nValues){
			float tmp = 0.0f;
			for (p = 0; p < nPoints; p++){
				tmp +=  distances[ i*nPoints + p ] * knownValues[ p*nValues + v];

				if (isinf(distances[ i*nPoints + p ])){ //isinf
					tmp = knownValues[ p*nValues + v];
					gridValues[ i*nValues + v ] = tmp;
					distan_zero = true;
					break;
				}
			}

			if(!distan_zero) {
				gridValues[ i*nValues + v ] = tmp / weightSum[i];
			}		
		}
	}
	// MY_DATA_TYPE sum_w_attr = 0;
	// for(int i=0; i<nGrids; ++i){
	// 	bool distan_zero = false;
	// 	for(int v=0; v<nValues; ++v){
	// 		// g1_1,g1_2,g1_3,...,g2_1,g2_2,g2_3
	// 		gridValues[ i*nValues + v ] = 0;

	// 		for(int p=0; p<nPoints; ++p){
	// 			gridValues[ i*nValues + v ] +=  distances[ i*nPoints + p ] * knownValues[ p*nValues + v];
				
	// 			if (isinf(distances[ i*nPoints + p ])){ //isinf
	// 				gridValues[ i*nValues + v ] = knownValues[ p*nValues + v];
	// 				distan_zero = true;
	// 				break;
	// 			}
	// 		}

	// 		if(!distan_zero) {
	// 			gridValues[ i*nValues + v ] /= weightSum[i];
	// 		}
	// 	}
	// }


};

__kernel
void computeDistances_SOA(const int DIM,const int nPoints,
	const __global float *knownCoords,const int nGrids,
	const __global float *gridCoords, __global float* distances){

	int d;
	int i = get_global_id(0);
	int j = get_global_id(1);

	if( i >= nGrids || j >= nPoints ) return;

	float tmp = 0.0f;
	for (d = 0; d < DIM; d++){
		// tmp += A[i*Ndim+k] * B[k*Pdim+j];
		tmp += pow( ( knownCoords[ j + d * nPoints ] - gridCoords[ i + nGrids * d ]) , 2.0f);
	}
	distances[ i + j * nGrids ] = pow( tmp , 0.5f );

	// d(g1,d1), d(g2,d1),...,d(g1,d2),d(g2,d2),...

	// for(int i=0;i<nGrids;++i){
	// 	for(int j=0;j< nPoints; ++j){

	// 		distances[ i + j * nGrids ] = 0;

	// 		for(int d=0; d< DIM; d++){
	// 			// change known coords and grid coords
	// 			distances[ i + j * nGrids ] += pow(( knownCoords[ j + d * nPoints ] 
	// 				- gridCoords[ i + nGrids * d ]) , 2);	
	// 		}

	// 		distances[ i + j * nGrids ] = pow(distances[ i + j * nGrids ], 0.5);
	// 	}
	// }

};

__kernel
void computeWeights_SOA(const int nGrids, const int nPoints, 
	__global float *distances, __global float *weightSum, const float p ){

	int j;
	int i = get_global_id(0);

	if( i >= nGrids ) return;

	for(j=0; j<nPoints; j++){
		distances[ i + j * nGrids ] = pow( 1.0f/distances[ i + j * nGrids ], p );	
		weightSum[i] += distances[ i + j * nGrids ];
	}
	
	// g1 = d(g1,d1)+d(g1,d2)+....

	// for( int i=0; i< nGrids; i++){
	// 	weightSum[i] = 0;
	// 	for(int j=0;j< nPoints; j++){

	// 		distances[ i + j * nGrids ] = pow(1/distances[ i + j * nGrids ], p);
	// 		weightSum[i] += distances[ i + j * nGrids ];
	// 		// weightSum[i] += pow(1/distances[ i*nPoints + j ], p);
	// 	}
	// }

};

__kernel 
void computeInterpolation_SOA(const int nValues, const int nGrids, 
	const int nPoints, const __global float *distances, 
	const __global float *weightSum, const __global float *knownValues, 
	__global float *gridValues){


	int p;
	int i = get_global_id(0);
	int v = get_global_id(1);

	if( i < nGrids  ) {
		bool distan_zero = false;
		if( v < nValues){
			float tmp = 0.0f;
			for (p = 0; p < nPoints; p++){
				tmp +=  distances[ i + p * nGrids ] * knownValues[ p + v * nPoints];

				if (isinf(distances[ i*nPoints + p ])){ //isinf
					tmp = knownValues[ p + v * nPoints ];
					gridValues[ i + v * nGrids ] = tmp;
					distan_zero = true;
					break;
				}
			}

			if(!distan_zero) {
				gridValues[ i + v * nGrids ] = tmp / weightSum[i];
			}		
		}
	}
	// MY_DATA_TYPE sum_w_attr = 0;
	// for(int i=0; i<nGrids; ++i){
	// 	bool distan_zero = false;
	// 	for(int v=0; v<nValues; ++v){
	// 		// g1_1,g2_1,g3_1,..., g1_2,g2_2,g3_2,...
			
	// 		gridValues[ i + v * nGrids  ] = 0;

	// 		for(int p=0; p<nPoints; ++p){
	// 			gridValues[ i + v * nGrids ] +=  distances[ i + p * nGrids ] * knownValues[ p + v * nPoints];
	// 			if (isinf(distances[ i*nPoints + p ]) ){
	// 				gridValues[ i + v * nGrids ] = knownValues[ p + v * nPoints ];
	// 				distan_zero = true;
	// 				break;
	// 			}
	// 		}
	// 		if(!distan_zero) {
	// 			gridValues[ i + v * nGrids ] /= weightSum[i];
	// 		}
	// 	}
	// }


};
