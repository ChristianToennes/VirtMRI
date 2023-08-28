all:
	#emcc --clear-cache; \
	emcc -O3 -Wall params.c simulate.c kspace_filter.c noise.c /mnt/x/webasm/bart/bart.o \
		 -I /mnt/x/webasm/bart/src/ -fblocks \
		 -s EXPORTED_FUNCTIONS="[\
		 					'_simulate',\
		 					'_malloc', '_make_cs_params', '_make_params', '_make_noise_params', '_make_filter_params', '_make_dataset', \
							'_free', '_free_cs_params', '_free_params', '_free_noise_params', '_free_filter_params', '_free_dataset', \
							'__Block_object_dispose', \
							'_num_init', '_pha_opts_defaults',\
							'_load_cfl', '_memcfl_unlink', '_memcfl_list_all', '_memcfl_load', '_memcfl_create', '_memcfl_exists', '_memcfl_unmap',\
							'_main_fft', '_main_ecalib', '_main_pics']" \
		 -s EXPORTED_RUNTIME_METHODS="['addFunction', 'removeFunction']" \
		 -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4gb -s RESERVED_FUNCTION_POINTERS=2 -s ASSERTIONS=0 \
		 /mnt/x/webasm/fftw/lib/libfftw3f.a /mnt/x/webasm/openblas/lib/libopenblas.a /mnt/x/webasm/libblocksruntime/lib/libBlocksRuntime.a \
	-o ../wwwroot/a.out.js
	