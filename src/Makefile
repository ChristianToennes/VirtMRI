all:
	cd compiled ; \
	#emcc --clear-cache; \
	emcc -O3 -Wall ../params.c ../simulate.c ../tinycs.c ../noise.c ../fft.c ../fft2.c ../fft3.c ../kissfft/kiss_fftnd.c ../kissfft/kiss_fft.c -I ../kissfft \
		 -s EXPORTED_FUNCTIONS="['_fft', '_ifft', '_fft2d', '_ifft2d', '_fft3d', '_ifft3d', '_fft_stride', \
		 						 '_compressed_sensing', '_make_cs_params', '_simulate', '_make_dataset', '_make_params', '_make_noise_params',\
								 '_malloc', '_free', '_free_dataset', '_free_params', '_free_noise_params', '_free_cs_params']" \
		 -s EXPORTED_RUNTIME_METHODS="['addFunction', 'removeFunction']" -s ALLOW_MEMORY_GROWTH=1 -s RESERVED_FUNCTION_POINTERS=2
