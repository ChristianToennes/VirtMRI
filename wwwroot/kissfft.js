self.importScripts("a.out.js");

function compressed_sensing(data, params) {
    var f_data = allocFromArray(Float64Array.from(data));
    var out = alloc(data.length);

    var xdim = Math.round(params["xdim"]);
    xdim = xdim > 0 ? xdim : k_xdim;
    xdim = xdim > k_xdim ? k_xdim : xdim;
    var ydim = Math.round(params["ydim"]);
    ydim = ydim > 0 ? ydim : k_ydim;
    ydim = ydim > k_ydim ? k_ydim : ydim;
    var zdim = Math.round(params["zdim"]);
    zdim = zdim > 0 ? zdim : k_zdim;
    zdim = zdim > k_zdim ? k_zdim : zdim;

    var ninner = "cs_ninner" in params ? parseInt(params["cs_ninner"]) : 1;
    var nbreg = "cs_nbreg" in params ? parseInt(params["cs_nbreg"]) : 10;
    var lambda = 1.0;
    var lambda2 = 0.3;
    var mu = "cs_mu" in params ? params["cs_mu"] : 1.0;
    var gam = 1;

    var params = _make_cs_params(xdim, ydim, zdim, ninner, nbreg, lambda, lambda2, mu, gam)

    _compressed_sensing(f_data.byteOffset, out.byteOffset, params);

    var img = new Float32Array(data.length);
    var tmp = new Float64Array(Module.HEAPU8.buffer, out.byteOffset, data.length);
    for(var i=0;i<data.length;i++) {
        img[i] = tmp[i];
    }
    
    free(params);
    free(f_data);
    free(out);

    return img;
}

/** Compute the FFT of a real-valued mxn matrix. */
function fft(data, m) {
    /* Allocate input and output arrays on the heap. */
    var heapData = allocFromArray(Float64Array.from(data));
    var heapSpectrum = alloc(2*m*8);

    _fft(heapData.byteOffset, heapSpectrum.byteOffset, m);

    /* Get spectrum from the heap, copy it to local array. */
    var spectrum = new Float64Array(m*2);
    spectrum.set(new Float64Array(Module.HEAPU8.buffer,
                                  heapSpectrum.byteOffset, m*2));

    /* Free heap objects. */
    free(heapData);
    free(heapSpectrum);

    return spectrum;
}

/** Compute the inverse FFT of a real-valued mxn matrix. */
function ifft(spectrum, m) {
    var heapSpectrum = allocFromArray(Float64Array.from(spectrum));
    var heapData = alloc(2*m*8);

    _ifft(heapSpectrum.byteOffset, heapData.byteOffset, m);

    var data = new Float64Array(m*2);
    data.set(new Float64Array(Module.HEAPU8.buffer,
                              heapData.byteOffset, m*2));

    for (i=0;i<m*2;i++) {
        data[i] /= m;
    }

    free(heapSpectrum);
    free(heapData);

    return data;
}

/** Compute the FFT of a real-valued mxn matrix. */
function kfft2d(data, m, n) {
    /* Allocate input and output arrays on the heap. */
    var heapData = allocFromArray(Float64Array.from(data));
    var heapSpectrum = alloc(2*m*n*8);

    _fft2d(heapData.byteOffset, heapSpectrum.byteOffset, m, n);

    /* Get spectrum from the heap, copy it to local array. */
    var spectrum = new Float64Array(m*n*2);
    spectrum.set(new Float64Array(Module.HEAPU8.buffer,
                                  heapSpectrum.byteOffset, m*n*2));

    /* Free heap objects. */
    free(heapData);
    free(heapSpectrum);

    return spectrum;
}

/** Compute the inverse FFT of a real-valued mxn matrix. */
function kifft2d(spectrum, m, n) {
    var heapSpectrum = allocFromArray(Float64Array.from(spectrum));
    var heapData = alloc(2*m*n*8);

    _ifft2d(heapSpectrum.byteOffset, heapData.byteOffset, m, n);

    var data = new Float64Array(m*n*2);
    data.set(new Float64Array(Module.HEAPU8.buffer,
                              heapData.byteOffset, m*n*2));

    for (i=0;i<m*n*2;i++) {
        data[i] /= m*n;
    }

    free(heapSpectrum);
    free(heapData);

    return data;
}

/** Compute the FFT of a real-valued mxn matrix. */
function kfft3d(data, m, n, l) {
    /* Allocate input and output arrays on the heap. */
    var heapData = allocFromArray(Float64Array.from(data));
    var heapSpectrum = alloc(2*m*n*l*8);

    _fft3d(heapData.byteOffset, heapSpectrum.byteOffset, m, n, l);

    /* Get spectrum from the heap, copy it to local array. */
    var spectrum = new Float64Array(m*n*l*2);
    spectrum.set(new Float64Array(Module.HEAPU8.buffer,
                                  heapSpectrum.byteOffset, m*n*l*2));

    /* Free heap objects. */
    free(heapData);
    free(heapSpectrum);

    return spectrum;
}


/** Compute the inverse FFT of a real-valued mxn matrix. */
function kifft3d(spectrum, m, n, l) {
    var heapSpectrum = allocFromArray(Float64Array.from(spectrum));
    var heapData = alloc(m*n*l*8);

    _ifft3d(heapSpectrum.byteOffset, heapData.byteOffset, m, n, l);

    var data = new Float64Array(m*n*l*2);
    data.set(new Float64Array(Module.HEAPU8.buffer,
                              heapData.byteOffset, m*n*l*2));

    for (i=0;i<m*n*l*2;i++) {
        data[i] /= m*n*l;
    }

    free(heapSpectrum);
    free(heapData);

    return data;
}

/** Create a heap array from the array ar. */
function allocFromArray(ar) {
    /* Allocate */
    var nbytes = ar.length * ar.BYTES_PER_ELEMENT;
    var heapArray = alloc(nbytes);

    /* Copy */
    heapArray.set(new Uint8Array(ar.buffer));
    return heapArray;
}


/** Allocate a heap array to be passed to a compiled function. */
function alloc(nbytes) {
    var ptr = Module._malloc(nbytes);
    return new Uint8Array(Module.HEAPU8.buffer, ptr, nbytes);
}


/** Free a heap array. */
function free(heapArray) {
    Module._free(heapArray.byteOffset);
}
