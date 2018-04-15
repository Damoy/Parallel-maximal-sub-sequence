/*
----------------------------------------------------

	-- Parallel maximal sub-sequence project --
	
	Author: FORNALI Damien
	Grade: Master I - IFI
	University: Nice-Sophia-Antipolis
	Year: 2017-2018
	Project subject link: https://sites.google.com/site/fabricehuet/teaching/parallelisme-et-distribution/sous-sequence-maximale

----------------------------------------------------
*/ 

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <omp.h>


/* - Structures -*/
struct array {
	long long* data;
	unsigned long long size;
};


/* 

	-- Declarations --

*/

/* - Core - */
void ascent(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long), long long bInit);
void pre_downhill(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long), long long bInit);
void suf_downhill(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long), long long bInit);
void ultimate(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long));
struct array* prefix(struct array* src, unsigned long long size, long long (*ibinary_fun)(long long, long long), long long bInit, unsigned long long organiseOffset);
struct array* suffix(struct array* src, unsigned long long size, long long (*ibinary_fun)(long long, long long), long long bInit, unsigned long long organiseOffset);
struct array* computeM(struct array* src, struct array* psum, struct array* ssum, struct array* smax, struct array* pmax);
long long findMax(struct array* M);
void findMSS(struct array* src, struct array* M);
void computeMSS(const char* file_name);

/* - Tools - */
struct array* parse(const char* file_name);
struct array* allocateArray(unsigned long long size);
long long plus(long long left, long long right);
long long max(long long left, long long right);
int isPow2(unsigned long long value);
void organise(struct array* src, unsigned long long offset);
void restore(struct array* src, unsigned long long noffset);
void destroy(struct array* arr);
void printSolution(long long max, struct array* src, unsigned long long startIndex, unsigned long long endIndex);
void lightPrintSolution(long long max);


/* 

	-- Definitions --

*/


/*
	- Core -
*/

/**
	Classic prefix / suffix ascent phase.

	a: the source array
	b: the destination array
	ibinary_fun: the int binary function called
	bInit: initially, the destination array is filled with this value

	- [Parallel] -
*/
inline void ascent(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long), long long bInit) {
	// first, fill the destination with bInit
	#pragma omp parallel for
	for(unsigned long long i = 0; i < a->size; ++i)
		b->data[i] = bInit;

	// then with the source content
	#pragma omp parallel for
  	for(unsigned long long i = a->size; i < b->size; ++i)
    	b->data[i] = a->data[i - a->size];

    /*  performs the tree ascent,
    	l varies from (m - 1) to 1 included */
  	for(unsigned long long l = (log2(a->size)) - 1; l > 0; --l){
  		// j varies from 2^l to (2^(l + 1) - 1) included
    	#pragma omp parallel for
    	for(unsigned long long j = 1 << l; j < (1 << (l + 1)); ++j)
      		b->data[j] = (*ibinary_fun)(b->data[j << 1], b->data[(j << 1) + 1]);
	}
}

/**
	Prefix algorithm downhill phase.

	a: the source array
	b: the destination array
	ibinary_fun: the int binary function called
	bInit: initially, the destination array is filled with this value

	e.g. of result tree:
						0
					 0	   15
				   0  10 15  16

	- [Parallel] -
*/
inline void pre_downhill(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long), long long bInit) {
	// fill the destination with bInit
	#pragma omp parallel for
	for(unsigned long long i = 0; i < b->size; ++i)
		b->data[i] = bInit;

	/*  performs the tree downhill
		l varies from 1 to (m- 1) included */
	for(unsigned long long l = 1; l < log2(a->size); ++l){
  		// j varies from 2^l to (2^(l + 1) - 1) included
		#pragma omp parallel for
		for(unsigned long long j = 1 << l; j < 1 << (l + 1); ++j){
			if(j % 2 == 0)
				b->data[j] = b->data[j >> 1];
      		else
				b->data[j] = (*ibinary_fun)(b->data[j >> 1], a->data[j - 1]);
		}
	}
}

/**
	Suffix algorithm downhill phase.
	The algorithm is similar to the prefix one but the resulted tree
	is the prefix's one mirror by the centered axe.

	a: the source array
	b: the destination array
	ibinary_fun: the int binary function called
	bInit: initially, the destination array is filled with this value

	e.g. of result tree:
					0
				15     0
			  16  15 10  0

	- [Parallel] -
*/
inline void suf_downhill(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long), long long bInit) {
	// fill the destination with bInit
	#pragma omp parallel for
	for(unsigned long long i = 0; i < b->size; ++i)
		b->data[i] = bInit;

	/*  performs the tree downhill
		l varies from 1 to (m- 1) included */
	for(unsigned long long l = 1; l < log2(a->size); ++l){
		 // j varies from 2^l to (2^(l + 1) - 1) included
		#pragma omp parallel for
		for(unsigned long long j = 1 << l; j < 1 << (l + 1); ++j){
			// odd part
			if(j % 2 == 1)
				b->data[j] = b->data[j >> 1];
			// even part
      		else   
				b->data[j] = (*ibinary_fun)(b->data[j >> 1], a->data[j + 1]);
		}
	}
}

/**
	Classic prefix / suffix ultimate phase.

	a: the source array
	b: the destination array
	ibinary_fun: the int binary function called
	
	- [Parallel] -
*/
inline void ultimate(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long)) {
	unsigned long long m = log2(a->size);
	// i varies from 2^(m - 1) to (2^m) - 1 included
	#pragma omp parallel for
	for(unsigned long long i = 1 << (m - 1); i < 1 << m; ++i)
		b->data[i] = (*ibinary_fun)(b->data[i], a->data[i]);
}

/**
	The prefix algorithm.

	src: the prefix source
	size: the new allocations size
	ibinary_fun: the int binary function used in the different phases
	bInit: value used to initially fill some arrays
	organiseOffset: value used to organise some arrays

	- [Parallel] -
*/
inline struct array* prefix(struct array* src, unsigned long long size, long long (*ibinary_fun)(long long, long long),
	long long bInit, unsigned long long organiseOffset){

	// create a
	struct array* a = allocateArray(size);
	// phase 1 with 'src' as the source and 'a' as the destination
	ascent(src, a, ibinary_fun, bInit);
	// create b
	struct array* b = allocateArray(size);
	// phase 2 with 'a' as the source and 'b' as the destination, prefix-specific
	pre_downhill(a, b, ibinary_fun, bInit);
	// phase 3 with 'a' as the source and 'b' as the destination
	ultimate(a, b, ibinary_fun);
	// free the 'a' array memory, not used anymore
	destroy(a);
	// re-organise the destination array before return
	organise(b, organiseOffset);

	return b;
}

/**
	The suffix algorithm.

	src: the prefix source
	size: the new allocations size
	ibinary_fun: the int binary function used in the different phases
	bInit: value used to initially fill some arrays
	organiseOffset: value used to organise some arrays

	- [Parallel] -
*/
inline struct array* suffix(struct array* src, unsigned long long size, long long (*ibinary_fun)(long long, long long),
	long long bInit, unsigned long long organiseOffset){

	// create a
	struct array* a = allocateArray(size);
	// phase 1 with 'src' as the source and 'a' as the destination
	ascent(src, a, ibinary_fun, bInit);
	// create b
	struct array* b = allocateArray(size);
	// phase 2 with 'a' as the source and 'b' as the destination, suffix-specific
	suf_downhill(a, b, ibinary_fun, bInit);
	// phase 3 with 'a' as the source and 'b' as the destination
	ultimate(a, b, ibinary_fun);
	// free the 'a' array memory, not used anymore
	destroy(a);
	// re-organise the destination array before return
	organise(b, organiseOffset);

	return b;
}

/**
	The MSS algorithn fifth phase.
	Assumes that the inputs has been re-organised, i.e. all have same size.

	src: the MSS algorithm input data
	psum: the src's sum-prefix
	ssum: the src's sum-suffix
	smax: the psum's max-suffix
	pmax: the ssum's max-prefix

	- [Parallel] -
*/
inline struct array* computeM(struct array* src, struct array* psum, struct array* ssum,
	struct array* smax, struct array* pmax){
	
	// the resulted M array
	struct array* M = allocateArray(src->size);

	// the M computing
	#pragma omp parallel for
	for(unsigned long long i = 0; i < src->size; ++i){
		M->data[i] = (pmax->data[i] - ssum->data[i] + src->data[i]) +
					(smax->data[i] - psum->data[i] + src->data[i]) -
					src->data[i];
	}

	return M;
}

/**
	The MSS algorithm sixth phase.
	Finds the maximum value in the inquired array.
	Uses an omp parallel reduction.

	src: the source array

	- [Parallel] -
*/
#pragma omp declare reduction(maximum : long long : omp_out = omp_in > omp_out ? omp_in : omp_out)
inline long long findMax(struct array* src){
	long long max = LONG_MIN;

	#pragma omp parallel for reduction(maximum:max)
	for(unsigned long long i = 0; i < src->size; ++i){
		if(src->data[i] > max)
			max = src->data[i];
	}

	// omp parallelizing returns 0 for
	// all negative data input
	if(max == 0){
		max = LONG_MIN;
		for(unsigned long long i = 0; i < src->size; ++i){
			if(src->data[i] > max)
				max = src->data[i];
		}
	}

	return max;
}

/**
	The last MSS algorithm phase.
	Finds the maximum sub-sequence in src according to the M array.

	src: the MSS algorithm input data
	M: the sixth phase result

	- [Parallel] -
*/
inline void findMSS(struct array* src, struct array* M){
	// get the maximum value of M
	long long maxValue = findMax(M);
	// the two indices used to locate the MSS
	long long startIndex = -1;
	long long endIndex = -1;

	if(src->size < 1){
		printf("Unknown source data.\n");
		return;
	}

	// handle the alone-element case
	if(src->size == 1){
		startIndex = 0;
		endIndex = 0;
	} else if(src->size == 2){
		// handle all possibilities of the two-elements case
		if(src->data[0] < 0){
			if(src->data[1] < 0){
				startIndex = (src->data[0] < src->data[1]) ? 1 : 0;
				endIndex = startIndex;
			} else {
				startIndex = 1;
				endIndex = 1;
			}
		} else {
			if(src->data[1] < 0){
				startIndex = 0;
				endIndex = 0;
			} else {
				startIndex = 0;
				endIndex = 1;
			}
		}
	} else {
		/*  handles the n elements
			where n > 2 and a power of two */

		// if the first element is the beginning of the sentence
		if(M->data[0] == maxValue)
			startIndex = 0;
		// if the last element is the end of the sentence
		if(M->data[M->size - 1] == maxValue)
			endIndex = M->size - 1;

		// if one of the indices is unset we iterate over M
		if(startIndex == -1 || endIndex == -1){
			#pragma omp parallel for
			for(unsigned long long i = 1; i < M->size - 1; ++i){
				// if the current M value is the maximum value
				if(M->data[i] == maxValue){

					// if the startIndex is unset and
					// the previous element is not the max value,
					// it means we are at the beginning of the MSS
					if(startIndex == -1 && M->data[i - 1] != maxValue)
						startIndex = i;

					// if the endIndex is unset and
					// the next element is not the max value,
					// it means we are at the end of the MSS
					if(endIndex == -1 && M->data[i + 1] != maxValue)
						endIndex = i;
				}
			}
		} else if(startIndex != -1 && endIndex != -1){
			/* 
				- Sequential loop -

				Allows to handle the "tricky mirror" case.
				See archive-root/resources/trickyMirror/trickyMirrorTest.txt
			*/
			long long sumCheck = 0;

			// updates the end index
			for(unsigned long long i = startIndex; i <= endIndex; ++i){
				sumCheck += src->data[i];
				if(sumCheck == maxValue){
					endIndex = i;
					break;
				}
			}
		}
	}

	/*  The following is for unhandled cases.
		For instance the case where the maximum value is only
		the source's first element:
			e.g. (src) 3 2 -7 11 */
	if(startIndex < 0 || endIndex < 0){
		// if the first element is equal to the max value
		if(M->data[0] == maxValue){
			startIndex = 0;
			endIndex = 0;
		} else if(M->data[M->size - 1] == maxValue){
			// if the last element is equal to the max value
			startIndex = M->size - 1;
			endIndex = startIndex;
		} else {
			printf("Error found while looking for a maximal sub-sequence.\n");
			return;
		}
	}

	// this is the final solution print !
	printSolution(maxValue, src, startIndex, endIndex);
	// lightPrintSolution(maxValue);
}

/**
	-- The MSS algorithm entry point --
	Computes the six phases, print the solution
	and free the memory used.

	file_name: the data source file name

	- [Parallel] -
*/
inline void computeMSS(const char* file_name){
	// get the source data
	struct array* src = parse(file_name);

	/*
		The data file will be considered as always well formed, i.e. composed of a sequence of relative integers separated by spaces.
		The table will be n = 2^m in size.

		if(!isPow2(src->size)){
			printf("Input file size must be a power of 2.\n");
			return;
		}
	*/

	/* Computes the prefixed sum (PSUM) - phase 1 */
	struct array* prefixedSum = prefix(src, src->size << 1, &plus, 0, src->size);

	/* Computes the suffixed sum (SSUM) - phase 2 */
	struct array* suffixedSum = suffix(src, src->size << 1, &plus, 0, src->size);

	/* Computes the suffixed max (SMAX) - phase 3 */
	long long suffixedMaxOffset = prefixedSum->size;
	struct array* suffixedMax = suffix(prefixedSum, prefixedSum->size << 1, &max, LONG_MIN, prefixedSum->size);

	/* Computes the prefixed max (PMAX) - phase 4 */
	long long prefixedMaxOffset = suffixedSum->size;
	struct array* prefixedMax = prefix(suffixedSum, suffixedSum->size << 1, &max, LONG_MIN, suffixedSum->size);

	/* Computes the M array - phase 5 */
	struct array* M = computeM(src, prefixedSum, suffixedSum, suffixedMax, prefixedMax);

	/* Finds the maximal subsequence sum - last phases (also prints the solution) */
	findMSS(src, M);

	// restore before destruction
	restore(prefixedSum, src->size);
	restore(suffixedSum, src->size);
	restore(suffixedMax, suffixedMaxOffset);
	restore(prefixedMax, prefixedMaxOffset);

	// free the memory
	destroy(src);
	destroy(prefixedSum);
	destroy(suffixedSum);
	destroy(suffixedMax);
	destroy(prefixedMax);
	destroy(M);
}

/*
	--Tools -
*/

/**
	Parses a file and produces an array from the int data collected.

	file_name: the data source file name

	- [Parallel] -
*/
inline struct array* parse(const char* file_name){
	FILE* file = fopen(file_name, "r");

	long long size = 0;
	long long ptr = 0;

	// iterates one first time to obtain the array size
	while(fscanf(file, "%lld", &ptr) == 1)
		size++;

	struct array* src;

	src = allocateArray(size);
	ptr = 0;
	rewind(file);

	long long tmp;

	// iterates a second time to obtain the array data
	while(ptr < size){
		fscanf(file, "%lld", &tmp);
		src->data[ptr++] = tmp;
	}

	fclose(file);
	return src;
}

/**
	Allocates a struct array.
	
	size: the number of elements of the generated array
*/
inline struct array* allocateArray(unsigned long long size) {
	struct array* tmp = malloc(sizeof(struct array));
	tmp->size = size;
	tmp->data = malloc(size * sizeof(long long));
	return tmp;
}

/** The plus binary fonction. */
inline long long plus(long long left, long long right){
	return left + right;
}

/** The max binary fonction. */
inline long long max(long long left, long long right){
	return (left > right) ? left : right;
}

/** Is 'value' a power of two ? (to check input data) */
inline int isPow2(unsigned long long value){
	return ceil(log2(value)) == floor(log2(value));
}

/**
	Moves forward the array base adress of an offset value. 
	Reduces its size accordingly.

	src: the array
	offset: the offset value
*/
inline void organise(struct array* src, unsigned long long offset){
	src->data += offset;
	src->size -= offset;
}

/**
	Moves backward the array base adress of an offset value. 
	Reduces its size accordingly.

	src: the array
	noffset: the offset value
*/
inline void restore(struct array* src, unsigned long long noffset){
	src->data -= noffset;
	src->size += noffset;
}

/** Free a struct array. */
inline void destroy(struct array* arr){
	free(arr->data);
	free(arr);
}

/**
	Prints the final solution.

	(May produce a segmentation fault with huge arrays, with such, use the alternative 'lightPrintSolution' function ?)

	max: the MSS algorithm max value computed
	src: the MSS algorithm input data
	startIndex: the start of the MSS
	endIndex: the end of the MSS
*/
inline void printSolution(long long max, struct array* src, unsigned long long startIndex, unsigned long long endIndex){
	printf("%lld ", max);

	for(unsigned long long i = startIndex; i < endIndex; ++i){
		printf("%lld ", src->data[i]);
	}

	printf("%lld\n", src->data[endIndex]);
}

/** Alternative, prints the inquired value. */
inline void lightPrintSolution(long long max){
	printf("%lld\n", max);
}

/* - The MSS algorithn entry point - */
int main(int argc, char** argv){
	if(argc != 2){
		printf("Wrong command arguments.\n");
		return -1;
	}

	// The input data is passed as an argument in the command line
	computeMSS(argv[1]);
	return 0;
}