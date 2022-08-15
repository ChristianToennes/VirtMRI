all:
	cd compiled ; \
	emcc -O3 -Wall ../simulate.c ../tinycs.c ../fft.c ../fft2.c ../fft3.c ../kissfft/kiss_fftnd.c ../kissfft/kiss_fft.c -I ../kissfft \
		 -s EXPORTED_FUNCTIONS="['_fft', '_ifft', '_fft2d', '_ifft2d', '_fft3d', '_ifft3d', '_fft_stride','_compressed_sensing', '_make_cs_params', '_malloc', '_free', '_simulate', '_make_dataset', '_make_params']" -s ALLOW_MEMORY_GROWTH=1 
