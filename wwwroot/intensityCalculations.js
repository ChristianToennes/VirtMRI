self.importScripts("./kissfft.js")
self.importScripts("./pako.min.js")

const xdim = 256;
const ydim = 256;
const zdim = 256;

var array_pd = new Float32Array(256 * 256);
var array_t1 = new Uint16Array(256 * 256);
var array_t2 = new Uint16Array(256 * 256);
var k_data_im_re, k_result;

async function loadDataSet() {
    array_pd = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function () {
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var a = new Uint8Array(resp, 8);
            var b = new Float32Array(a.length);
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 256.0 + mm[0]) * (mm[1] - mm[0]);
            }
            resolve(b);
        }
        xhr.onerror = reject
        xhr.open("GET", "pd.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_t1 = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function () {
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var a = new Uint8Array(resp, 8);
            var b = new Float32Array(a.length);
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 256.0 + mm[0]) * (mm[1] - mm[0]);
            }
            resolve(b);
        }
        xhr.onerror = reject
        xhr.open("GET", "t1.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_t2 = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function () {
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var a = new Uint8Array(resp, 8);
            var b = new Float32Array(a.length);
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 256.0 + mm[0]) * (mm[1] - mm[0]);
            }
            resolve(b);
        }
        xhr.onerror = reject
        xhr.open("GET", "t2.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    return [array_pd, array_t1, array_t2];
}

function calcKSpace(result) {
    k_result = new Uint8ClampedArray(xdim * ydim * zdim);
    k_data_im_re = new Float32Array(xdim * ydim * zdim * 2);
    for (var z = 0; z < zdim; z++) {
        var slice_data = new Float32Array(xdim * ydim)
        for (var x = 0; x < xdim * ydim; x++) {
            slice_data[x] = result[x + z * xdim * ydim]
        }

        var fft_res = rfft2d(slice_data, xdim, ydim);
        k_data_im_re.set(fft_res, z * xdim * ydim * 2)

        k_result.set(transformKSpace(fft_res), z * xdim * ydim);
    }
    return k_result
}

function transformKSpace(fft_res) {
    var k_data = new Float32Array(xdim * ydim);
    var k_result = new Float32Array(xdim * ydim);
    var maxval = 0;
    var minval = 999999999;
    for (var i = 0; i < xdim * ydim; i++) {
        k_data[i] = Math.sqrt(fft_res[2 * i] * fft_res[2 * i] + fft_res[2 * i + 1] * fft_res[2 * i + 1]);
        if (k_data[i] == -Infinity) {
            k_data[i] = 0;
        }
        maxval = Math.max(maxval, k_data[i]);
        minval = Math.min(minval, k_data[i]);
    }

    for (var i = 0; i < (xdim / 2 + 1) * ydim; i++) {
        var val = (k_data[i] - minval) * 255 / maxval;

        var y = (Math.floor(i / (xdim / 2 + 1)) + ydim / 2) % ydim;
        var x = i % (xdim / 2 + 1) + xdim / 2;

        var j = y * xdim + x - 1;
        var k = y * xdim + (xdim - x);

        k_result[j] = val;
        k_result[k] = val;
    }
    return k_result;
}

function genMapKSpace(xlines, ylines, fmin, fmax) {
    var dx = xdim / xlines;
    var dy = xdim / ylines;
    return function filterKSpace(value, d_index, array) {
        var index = Math.floor((d_index - 1) / 2);
        var y = (Math.floor(index / (xdim / 2 + 1))) % ydim;
        if(y>ydim/2) { 
            y = ydim-y;
        }
        var x = index % (xdim / 2 + 1);
        if(index > ydim*(xdim+2)/2) {
            x = x + xdim / 2 + 1;
        }
        var res1 = Math.sqrt(x*x+y*y) >= fmin && Math.sqrt(x*x+y*y) <= fmax;
        var res = res1 && (Math.ceil(x / dx) * dx - x) < 1 && (Math.ceil(y / dy) * dy - y) < 1;
        return res ? value : 0;
    }
}

function genFilterKSpace(xlines, ylines, fmin, fmax) {
    var dx = xdim / xlines;
    var dy = xdim / ylines;
    return function filterKSpace(value, d_index, array) {
        var index = Math.floor((d_index - 1) / 2);
        var y = (Math.floor(index / (xdim / 2 + 1))) % ydim;
        if(y>ydim/2) { 
            y = ydim-y;
        }
        var x = index % (xdim / 2 + 1);
        if(index > ydim*(xdim+2)/2) {
            x = x + xdim / 2 + 1;
        }
        var res1 = Math.sqrt(x*x+y*y) >= fmin && Math.sqrt(x*x+y*y) <= fmax;
        var res = res1 && (Math.ceil(x / dx) * dx - x) < 1 && (Math.ceil(y / dy) * dy - y) < 1;
        return res ? true : false;
    }
}

function inverseKSpace(kSpace, xlines, ylines, fmin, fmax) {
    var filterKSpace = genFilterKSpace(xlines, ylines, 0, xdim*xdim)
    var mapKSpace = genMapKSpace(xlines, ylines, fmin, fmax)
    var result = new Float32Array(xdim * ydim * zdim);
    var k_result = new Float32Array(xdim * ydim * zdim);
    for (var z = 0; z < zdim; z++) {
        var slice_data = kSpace.subarray(z * xdim * ydim * 2, (z + 1) * xdim * ydim * 2);

        var input_data = slice_data.map(mapKSpace);

        k_result.set(transformKSpace(input_data), z * xdim * ydim);

        //input_data = input_data.filter(filterKSpace);

        var img_result = irfft2d(input_data, xdim, ydim)
        var maxval = 0;
        var minval = 999999999;
        for (var i = 0; i < img_result.length; i++) {
            maxval = Math.max(maxval, img_result[i]);
            minval = Math.min(minval, img_result[i]);
        }
        for (var i = 0; i < img_result.length; i++) {
            result[i + z * xdim * ydim] = (img_result[i] - minval) / (maxval - minval)
            result[i + 1 + z * xdim * ydim] = (img_result[i] - minval) / (maxval - minval)
        }
    }
    return [result, k_result];
}

function calcSpinEcho(te, tr) {
    var result = new Float32Array(array_t1.length);
    for (var x = 0; x < result.length; x++) {
        var t1 = array_t1[x]
        var t2 = array_t2[x]
        var pd = array_pd[x]
        if (t1 == 0) {
            t1 = 1;
        }
        if (t2 == 0) {
            t2 = 1;
        }
        var val = ((pd) * (Math.exp((-te) / (t2))) * (1 - Math.exp((-tr) / (t1))));

        result[x] = val
    }

    var k_result = calcKSpace(result);
    return [result, k_result];
}

function calcInversionRecovery(ti, tr) {
    var result = new Float32Array(array_t1.length);
    for (var x = 0; x < result.length; x++) {
        // Berechnen der Intensitaetsmatrix
        var t1 = array_t1[x]
        var t2 = array_t2[x]
        var pd = array_pd[x]
        // Divison durch 0 abfangen
        if (t1 == 0) {
            t1 = 1.0;
        }
        var val = Math.abs(pd * (1.0 - 2.0 * Math.exp(-ti / t1) + Math.exp(-tr / ti)));

        result[x] = val
    }
    var k_result = calcKSpace(result);
    return [result, k_result];
}

var queryableFunctions = {
    spinEcho: function (te, tr) {
        reply('result', calcSpinEcho(te, tr));
    },
    inversionRecovery: function (ti, tr) {
        reply('result', calcInversionRecovery(ti, tr));
    },
    reco: function (xlines, ylines, fmin, fmax) {
        reply('result', inverseKSpace(k_data_im_re, xlines, ylines, fmin, fmax));
    },
    loadData: async function () {
        reply('loadData', await loadDataSet());
    }
};

// system functions

function defaultReply(message) {
    // your default PUBLIC function executed only when main page calls the queryableWorker.postMessage() method directly
    // do something
}

function reply() {
    if (arguments.length < 1) {
        throw new TypeError('reply - not enough arguments');
        return;
    }
    postMessage({
        'queryMethodListener': arguments[0],
        'queryMethodArguments': Array.prototype.slice.call(arguments, 1)
    });
}

onmessage = async function (oEvent) {
    if (oEvent.data instanceof Object && oEvent.data.hasOwnProperty('queryMethod') && oEvent.data.hasOwnProperty('queryMethodArguments')) {
        queryableFunctions[oEvent.data.queryMethod].apply(self, oEvent.data.queryMethodArguments);
    } else {
        defaultReply(oEvent.data);
    }
};