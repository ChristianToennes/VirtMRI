all:
	#emcc --clear-cache; \
	emcc -O3 -Wall params.c simulate.c tinycs.c kspace_filter.c noise.c ../../../bla/bart.o \
		 -I ../../bart/src/ -fblocks \
		 -s EXPORTED_FUNCTIONS="['_compressed_sensing', '_make_cs_params', '_simulate', '_make_dataset', '_make_params', '_make_noise_params', '_make_filter_params',\
								 '_malloc', '_free', '_free_dataset', '_free_params', '_free_noise_params', '_free_cs_params', '__Block_object_dispose', \
								 '_bart_version', '_phantom',  '_calc_phantom', '_calc_bart', '_calc_circ', '_fft', '_ifft','_ifftc','_fftc', '_fftscale', \
								 '_num_init', '_pha_opts_defaults', '_img_sim', '_kspace_sim', '_sen_sim', '_kspace_cs', '_kspace_cs', \
								 '_img_cs', '_kspace_filt', '_img_filt', '_load_cfl', '_memcfl_unlink', '_memcfl_list_all', '_memcfl_load', \
								 '_memcfl_create', '_memcfl_exists', '_memcfl_unmap', '_main_fft', '_main_ecalib', '_main_pics']" \
		 -s EXPORTED_RUNTIME_METHODS="['addFunction', 'removeFunction']" -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4gb -s RESERVED_FUNCTION_POINTERS=2 -s ASSERTIONS=0 \
		 ../../../bla/libfftw3f.a ../../../bla/libopenblas.a ../../../bla/libBlocksRuntime.a \
	-o ../wwwroot/a.out.js
	