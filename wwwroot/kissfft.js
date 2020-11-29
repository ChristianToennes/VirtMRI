import Module from "./a.out.js"

let instance = {
    ready: new Promise(resolve => {
        Module({
            onRuntimeInitialized() {
                //instance = Promise.resolve();
                instance = Object.assign(this, {
                    ready: Promise.resolve()
                });
                resolve();
            }
        });
    })
};

/** Compute the FFT of a real-valued mxn matrix. */
export function rfft2d(data, m, n) {
    /* Allocate input and output arrays on the heap. */
    var heapData = allocFromArray(data);
    var heapSpectrum = alloc(2*m*n*4);

    instance._rfft2d(heapData.byteOffset, heapSpectrum.byteOffset, m, n);

    /* Get spectrum from the heap, copy it to local array. */
    var spectrum = new Float32Array(m*n*2);
    spectrum.set(new Float32Array(instance.HEAPU8.buffer,
                                  heapSpectrum.byteOffset, m*n*2));

    /* Free heap objects. */
    free(heapData);
    free(heapSpectrum);

    return spectrum;
}


/** Compute the inverse FFT of a real-valued mxn matrix. */
export function irfft2d(spectrum, m, n) {
    var heapSpectrum = allocFromArray(spectrum);
    var heapData = alloc(m*n*4);

    instance._irfft2d(heapSpectrum.byteOffset, heapData.byteOffset, m, n);

    var data = new Float32Array(m*n);
    data.set(new Float32Array(instance.HEAPU8.buffer,
                              heapData.byteOffset, m*n));

    for (i=0;i<m*n;i++) {
        data[i] /= m*n;
    }

    free(heapSpectrum);
    free(heapData);

    return data;
}


/** Create a heap array from the array ar. */
export function allocFromArray(ar) {
    /* Allocate */
    var nbytes = ar.length * ar.BYTES_PER_ELEMENT;
    var heapArray = alloc(nbytes);

    /* Copy */
    heapArray.set(new Uint8Array(ar.buffer));
    return heapArray;
}


/** Allocate a heap array to be passed to a compiled function. */
export function alloc(nbytes) {
    var ptr = instance._malloc(nbytes);
    return new Uint8Array(instance.HEAPU8.buffer, ptr, nbytes);
}


/** Free a heap array. */
export function free(heapArray) {
    instance._free(heapArray.byteOffset);
}