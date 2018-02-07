

int index2D(int i, int j, int max_j){
    return (i*max_j) + j;
}

int index3D(int i, int j, int k, int max_j, int max_k){
    return (i*max_j*max_k) + (j*max_k) + k ;
}

int indexNDisvalid(int ndimensions, int *indexes, int *size_dimensions){
    int i;

    for(i=0;i<ndimensions;i++){
        if(indexes[i]<0 || indexes[i]>=size_dimensions[i]){ //if the index is out of bound...
            return -1; // return error
        }
    }
    
    return 0;    
}

int indexND(int ndimensions, int *indexes, int *size_dimensions){
    int i,j,sum=0,multi;
    // if indexes > maximum error..
    // ndimiensions >=0
    for(i=0;i<ndimensions;i++){
        multi=1;
        for(j=i+1;j<ndimensions;j++){
            multi=multi*size_dimensions[j];
        }
        sum=sum+(indexes[i]*multi);
    }
    
    return sum;
}

void indexesofND(int index, int *indexes, int ndimensions,  int *size_dimensions){
    int i,j,rest=0,multi;
    // if indexes > maximum error..
    // ndimiensions >=0
    
    rest=index;
    for(i=0;i<ndimensions;i++){
        multi=1;
        for(j=i+1;j<ndimensions;j++){
            multi=multi*size_dimensions[j];
        }
        indexes[i]=(int)(rest/multi);
        rest=rest-(indexes[i]*multi);
    }
}
