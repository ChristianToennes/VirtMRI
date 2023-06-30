all:
	#emcc --clear-cache; \
	emcc -O3 -Wall params.c simulate.c tinycs.c kspace_filter.c noise.c bart.o \
		-I kissfft -I ../../bart/src/ \
		 -s EXPORTED_FUNCTIONS="['_compressed_sensing', '_make_cs_params', '_simulate', '_make_dataset', '_make_params', '_make_noise_params', '_make_filter_params',\
								 '_malloc', '_free', '_free_dataset', '_free_params', '_free_noise_params', '_free_cs_params', \
								 '_bart_version', '_phantom',  '_calc_phantom', '_calc_bart', '_calc_circ', '_fft', '_ifft','_ifftc','_fftc','_num_init', '_pha_opts_defaults']" \
		 -s EXPORTED_RUNTIME_METHODS="['addFunction', 'removeFunction']" -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s RESERVED_FUNCTION_POINTERS=2 \
	-o ../wwwroot/a.out.js

	#emcc -O3 -Wall bart.o \
	#	-I kissfft \
	#	 -s EXPORTED_FUNCTIONS="['_malloc','_free','_bart_version', '_phantom',  '_calc_phantom', '_calc_bart', '_calc_circ', '_fft', '_ifft', '_num_init', '_pha_opts_defaults']" \
	#	 -s EXPORTED_RUNTIME_METHODS="['addFunction', 'removeFunction']" -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s RESERVED_FUNCTION_POINTERS=2 \
	#-o ../wwwroot/bart.js