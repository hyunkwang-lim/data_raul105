// Return linear array of two dimension
// example: 
// int v[I_MAX=3][J_MAX=5];
// v[1][3]=3
// v[ i=8 ] => 3
//  return i=8, index where you can find element 1,3
int index2D(int i, int j, int max_j);

// Return linear array of three dimension. See index2D
int index3D(int i, int j, int k, int max_j, int max_k);


// Return linear index of array of ndimiension where the indexes are specified in indexes argument and the number of element of each dimension is given in size_dimensions
int indexND(int ndimensions, int *indexes, int *size_dimensions);

// Return indexes of a multidimensional array given the linear index position and specifying the shape of the array: ndimiension is the number of dimensions and size_dimensions is the number of elements of each dimension. This function is the reverse of indexND
void indexesofND(int index, int *indexes, int ndimensions,  int *size_dimensions);


// This function check if the indexes are in the bounds of the array, returning 0. Otherwise it will return -1. This function should be called before indexND or indexesofND
int indexNDisvalid(int ndimensions, int *indexes, int *size_dimensions);
