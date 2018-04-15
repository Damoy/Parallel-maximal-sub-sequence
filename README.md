# Parallel-maximal-sub-sequence
Parallel maximal-sub-sequence using OMP in C

----------------------------------------------------

	-- Parallel maximal sub-sequence project --
	
	Author: FORNALI Damien
	Grade: Master I - IFI
	University: Nice-Sophia-Antipolis
	Year: 2017-2018
	Project subject link: https://sites.google.com/site/fabricehuet/teaching/parallelisme-et-distribution/sous-sequence-maximale

----------------------------------------------------


I. Project archive content
	(-- is root, -> is folder)

	-- fornali.c 														: fully documented max sub-sequence algorithm
	
	-- resources														: some tests that are not covered by the provided ones
			-> trickyMirror -> trickyMirrorTest							: the 'tricky mirror' test illustration
							-> trickyMirrorTest.txt						: tricky mirror pattern definition and explanation

			-> bounds 		-> boundsTest.txt							: solutions on the above tests
							-> lowerBoundTest							: solution is the first data element test
							-> upperBoundTest							: solution is the last data element test
							-> allBoundsTest							: solution contains all the data
							-> noBoundsTest								: solution excludes the lower and upper bounds
							-> sizeBounds -> lightTwoElementsTest		: two elements test
										  -> heavyTwoElementsTest		: two long elements test

	-- README.txt 														: this




II. Parallelization information

	II.1 Parallelized functions
		
		> Functions that are (almost) completely parallelized. (omp parallel for and omp parallel reduction)

		1. void ascent(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long), long long bInit);
		2. void pre_downhill(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long), long long bInit);
		3. void suf_downhill(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long), long long bInit);
		4. void ultimate(struct array* a, struct array* b, long long (*ibinary_fun)(long long, long long));
		5. struct array* prefix(struct array* src, unsigned long long size, long long (*ibinary_fun)(long long, long long), long long bInit, long long organiseOffset);
		6. struct array* suffix(struct array* src, unsigned long long size, long long (*ibinary_fun)(long long, long long), long long bInit, long long organiseOffset);
		7. struct array* computeM(struct array* src, struct array* psum, struct array* ssum, struct array* smax, struct array* pmax);
		8. long long findMax(struct array* M);
		9. void findMSS(struct array* src, struct array* M);
		10. void computeMSS(const char* file_name);

	II.2 Sequential ones

		> Fully sequential functions.

		1. struct array* parse(const char* file_name);
		2. struct array* allocateArray(unsigned long long size);
		3. long long plus(long long left, long long right);
		4. long long max(long long left, long long right);
		5. void destroy(struct array* arr);
		6. void printSolution(long long max, struct array* src, long long startIndex, long long endIndex);
		7. void lightPrintSolution(long long max);


III. Some explanations

	The algorithm is quite classic, however here are some information that might be relevant:

	1.	Calling the prefix (/suffix) algorithm generates 2-times-bigger arrays than the input ones. 
	Thus, calling the prefix (/suffix) algorithm on generated data (by prefix / suffix) increase exponentially the memory space usage. To avoid handling complex and miscellaneous indices according to the data processed, I 'organise' the generated data, that is to say, I change the start address of the array to be equal to the one of the start of the relevant data (pointer offset, organise function). However, this requires to restore the array initial base address and size before releasing memory or it would throw a segmentation fault error.

	2.  The suffix algorithm is implemented following the prefix one.
	The only difference is the downhill, which creates a mirror tree of the prefix one.

	3.	As for the performance, for int 2^n calculations I make use of the shift operator. Even if gcc is smart enough to replace, for instance, heavy operations (division, 125 / 16) by ligther ones (bit shift, 125 >> 4), I prefer to write them as much as possible. Ignoring this improvements assumes that the user will compile with gcc which might not be the case. I also reduce as much as possible the use of 'malloc', useless variables and data and I try to maximize the parallelization.

	4. For data, I use the 'long long' primitive type to be able to accept and handle as much data as possible.

	5. It was quite hard to try to handle a maximum of input patterns. As a result the 'findMSS' method is quite big. Moreover, I created the 'resources' folder which contain tests patterns that are not covered by the three tests provided.
