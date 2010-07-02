
//#include <assert.h>
//#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <algorithm>
using std::min;
//#include <pthread.h>


#define _ESVCPTR
#include <complex>
#include <essl.h>



int main(){
    int blocksize = 512;
    double *m1, *m2, *m3;

      dgemm( "N", "N",
	     blocksize, blocksize, blocksize,
	     -1.0, m1,
	     blocksize, m2, blocksize,
	     1.0, m3, blocksize);

      return 0;
}

