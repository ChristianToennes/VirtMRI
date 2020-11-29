import * as pako from "./pako.esm.js"
import * as kissfft from "./kissfft.js"

const xdim = 256;
const ydim = 256;
const zdim = 256;

var array_pd = new Float32Array(256 * 256);
var array_t1 = new Uint16Array(256 * 256);
var array_t2 = new Uint16Array(256 * 256);

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
    return result;
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

    var k_result = new Uint8ClampedArray(xdim * ydim * zdim);
    for (var z = 0; z < zdim; z++) {
        var slice_data = new Float32Array(xdim * ydim)
        for (var x = 0; x < xdim * ydim; x++) {
            slice_data[x] = result[x + z * xdim * ydim]
        }

        var k_data_im_re = kissfft.rfft2d(slice_data, xdim, ydim)
        var k_data = new Float32Array(xdim * ydim);
        var maxval = 0;
        var minval = 999999999;
        for (var i = 0; i < xdim * ydim; i++) {
            k_data[i] = Math.sqrt(k_data_im_re[2 * i] * k_data_im_re[2 * i] + k_data_im_re[2 * i + 1] * k_data_im_re[2 * i + 1]);
            if (k_data[i] == -Infinity) {
                k_data[i] = 0;
            }
            maxval = Math.max(maxval, k_data[i]);
            minval = Math.min(minval, k_data[i]);
        }
        
        var x = xdim / 2, y = ydim / 2
        for (var i = 0; i < (xdim / 2 + 1) * ydim; i++) {
            var val = (k_data[i] - minval) * 255 / maxval;

            var j = y * xdim + x;
            var k = (y + 1) * xdim - (x + 1);

            k_result[j + z * xdim * ydim] = val;
            k_result[k + z * xdim * ydim] = val;

            x++;
            if (x == xdim + 1) {
                x = xdim / 2;
                y++;
                if (y == ydim) {
                    y = 0;
                }
            }
        }
    }



    return [result, k_result];
}

var queryableFunctions = {
    spinEcho: function (te, tr) {
        reply('result', calcSpinEcho(te, tr));
    },
    inversionRecovery: function (ti, tr) {
        reply('result', calcInversionRecovery(ti, tr));
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
    if (arguments.length < 1) { throw new TypeError('reply - not enough arguments'); return; }
    postMessage({ 'queryMethodListener': arguments[0], 'queryMethodArguments': Array.prototype.slice.call(arguments, 1) });
}

onmessage = async function (oEvent) {
    if (oEvent.data instanceof Object && oEvent.data.hasOwnProperty('queryMethod') && oEvent.data.hasOwnProperty('queryMethodArguments')) {
        queryableFunctions[oEvent.data.queryMethod].apply(self, oEvent.data.queryMethodArguments);
    } else {
        defaultReply(oEvent.data);
    }
};