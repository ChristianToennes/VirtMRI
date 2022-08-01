self.importScripts("a.out.js");

const scalar_size = 4;

function compressed_sensing_fast(data, params) {
    
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
    var nbreg = "cs_nbreg" in params ? parseInt(params["cs_nbreg"]) : 80;
    var lambda = "cs_lambda1" in params ? parseFloat(params["cs_lambda1"]):1.0;
    var lambda2 = "cs_lambda2" in params? parseFloat(params["cs_lambda2"]):0.15;
    var mu = "cs_mu" in params ? params["cs_mu"] : 1.0;
    var gam = "cs_gam" in params? parseFloat(params["cs_gam"]):1;
    
    var params = _make_cs_params(xdim, ydim, zdim, ninner, nbreg, lambda, lambda2, mu, gam)
    var f_data = scalar_size == 4 ? allocFromArray(data) : allocFromArray(Float64Array.from(data));
    var f_data_ptr = f_data.byteOffset;
    var out = alloc(xdim*ydim*zdim*scalar_size);
    var out_ptr = out.byteOffset;

    _compressed_sensing(f_data_ptr, out_ptr, params);

    var img = new Float32Array(xdim*ydim*zdim);
    if(scalar_size == 8) {
        var tmp = new Float64Array(img.length);
        tmp.set(new Float64Array(Module.HEAPU8.buffer, out_ptr, img.length));

        for(var i=0;i<img.length;i++) {
            img[i] = tmp[i];
        }
    } else {
        img.set(new Float32Array(Module.HEAPU8.buffer, out_ptr, img.length));
    }

    var kspace = new Float32Array(xdim*ydim*zdim*2);
    if(scalar_size == 8) {
        var tmp = new Float64Array(kspace.length);
        tmp.set(new Float64Array(Module.HEAPU8.buffer, f_data_ptr, kspace.length));

        for(var i=0;i<img.length;i++) {
            kspace[i] = tmp[i];
        }
    } else {
        kspace.set(new Float32Array(Module.HEAPU8.buffer, f_data_ptr, kspace.length));
    }
    
    free(params);
    free(f_data);
    free(out);

    return [img, kspace];
}

/** Compute the FFT of a real-valued mxn matrix. */
function fft(data, m) {
    /* Allocate input and output arrays on the heap. */
    var heapData = scalar_size==4 ? allocFromArray(data) : allocFromArray(Float64Array.from(data));
    var heapSpectrum = alloc(2*m*scalar_size);

    _fft(heapData.byteOffset, heapSpectrum.byteOffset, m);

    /* Get spectrum from the heap, copy it to local array. */
    var spectrum = new Float32Array(m*2);
    if(scalar_size==8) {
        var tmp = new Float64Array(Module.HEAPU8.buffer, heapSpectrum.byteOffset, m*2);
        for(var i=0;i<tmp.length;i++) { spectrum[i] = tmp[i]; }
    } else {
        spectrum.set(new Float32Array(Module.HEAPU8.buffer, heapSpectrum.byteOffset, m*2));
    }   

    /* Free heap objects. */
    free(heapData);
    free(heapSpectrum);

    return spectrum;
}

/** Compute the inverse FFT of a real-valued mxn matrix. */
function ifft(spectrum, m) {
    var heapSpectrum = scalar_size==4 ? allocFromArray(spectrum) : allocFromArray(Float64Array.from(spectrum));
    var heapData = alloc(2*m*scalar_size);

    _ifft(heapSpectrum.byteOffset, heapData.byteOffset, m);

    var data = scalar_size==4 ? Float32Array.from(new Float32Array(Module.HEAPU8.buffer,heapData.byteOffset, m*2)): Float32Array.from(new Float64Array(Module.HEAPU8.buffer,heapData.byteOffset, m*2));

    //for (i=0;i<m*2;i++) {
        //data[i] /= m;
    //}

    free(heapSpectrum);
    free(heapData);

    return data;
}

/** Compute the FFT of a real-valued mxn matrix. */
function kfft2d(data, m, n) {
    /* Allocate input and output arrays on the heap. */
    var heapData = scalar_size==4?allocFromArray(data):allocFromArray(Float64Array.from(data));
    var heapSpectrum = alloc(2*m*n*scalar_size);

    _fft2d(heapData.byteOffset, heapSpectrum.byteOffset, m, n);

    /* Get spectrum from the heap, copy it to local array. */
    var spectrum = scalar_size==4?Float32Array.from(new Float32Array(Module.HEAPU8.buffer, heapSpectrum.byteOffset, m*n*2)):Float32Array.from(new Float64Array(Module.HEAPU8.buffer, heapSpectrum.byteOffset, m*n*2));

    /* Free heap objects. */
    free(heapData);
    free(heapSpectrum);

    return spectrum;
}

/** Compute the inverse FFT of a real-valued mxn matrix. */
function kifft2d(spectrum, m, n) {
    var heapSpectrum = scalar_size==4?allocFromArray(spectrum):allocFromArray(Float64Array.from(spectrum));
    var heapData = alloc(2*m*n*scalar_size);

    _ifft2d(heapSpectrum.byteOffset, heapData.byteOffset, m, n);

    var data = scalar_size==4?Float32Array.from(new Float32Array(Module.HEAPU8.buffer, heapData.byteOffset, m*n*2)):Float32Array.from(new Float64Array(Module.HEAPU8.buffer, heapData.byteOffset, m*n*2));

    //for (i=0;i<m*n*2;i++) {
        //data[i] /= m*n;
    //}

    free(heapSpectrum);
    free(heapData);

    return data;
}

/** Compute the FFT of a real-valued mxn matrix. */
function kfft3d(data, m, n, l) {
    /* Allocate input and output arrays on the heap. */
    var heapData = allocFromArray(Float64Array.from(data));
    var heapSpectrum = alloc(2*m*n*l*scalar_size);

    _fft3d(heapData.byteOffset, heapSpectrum.byteOffset, m, n, l);

    /* Get spectrum from the heap, copy it to local array. */
    var spectrum = scalar_size==4?Float32Array.from(new Float32Array(Module.HEAPU8.buffer, heapSpectrum.byteOffset, m*n*l*2)):Float32Array.from(new Float64Array(Module.HEAPU8.buffer, heapSpectrum.byteOffset, m*n*l*2));

    /* Free heap objects. */
    free(heapData);
    free(heapSpectrum);

    return spectrum;
}


/** Compute the inverse FFT of a real-valued mxn matrix. */
function kifft3d(spectrum, m, n, l) {
    var heapSpectrum = scalar_size==4?allocFromArray(spectrum):allocFromArray(Float64Array.from(spectrum));
    var heapData = alloc(m*n*l*scalar_size);

    _ifft3d(heapSpectrum.byteOffset, heapData.byteOffset, m, n, l);

    var data = scalar_size==4?Float32Array.from(new Float32Array(Module.HEAPU8.buffer,heapData.byteOffset, m*n*l*2)):Float32Array.from(new Float64Array(Module.HEAPU8.buffer,heapData.byteOffset, m*n*l*2));

    //for (i=0;i<m*n*l*2;i++) {
        //data[i] /= m*n*l;
    //}

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
