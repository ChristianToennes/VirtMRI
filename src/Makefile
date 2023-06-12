all:
	#emcc --clear-cache; \
	emcc -O3 -Wall params.c simulate.c tinycs.c kspace_filter.c noise.c fft.c fft2.c fft3.c kissfft/kiss_fftnd.c kissfft/kiss_fft.c bart.o \
		-I kissfft \
		 -s EXPORTED_FUNCTIONS="['_kfft', '_kifft', '_kfft2d', '_kifft2d', '_kfft3d', '_kifft3d', '_kfft_stride', \
		 						 '_compressed_sensing', '_make_cs_params', '_simulate', '_make_dataset', '_make_params', '_make_noise_params', '_make_filter_params',\
								 '_malloc', '_free', '_free_dataset', '_free_params', '_free_noise_params', '_free_cs_params', '_bart_version']" \
		 -s EXPORTED_RUNTIME_METHODS="['addFunction', 'removeFunction']" -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s RESERVED_FUNCTION_POINTERS=2 \
	-o ../wwwroot/a.out.js
