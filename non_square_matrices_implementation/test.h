#ifndef __TEST__


double test(void (*f)(float **,
	                  float const *const *const,
	                  float const *const *const,
	                  size_t, size_t, size_t), 
	        float **C, float** A, float **B, size_t A_rows, size_t A_columns, size_t B_columns);

#endif // __TEST__
