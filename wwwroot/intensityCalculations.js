self.importScripts("./kissfft.js")
self.importScripts("./pako.min.js")

var xdim = 256;
var ydim = 256;
var zdim = 256;

const na_tau1 = 50;
const na_tau2 = 0.1;
const na_t1 = 39.2;
const na_t2s = 50;
const na_t2f = 4;
const na_t2fr = 60;
const na_mm = 140;
const na_vol = 0.7;

var array_pd = new Float32Array(256 * 256);
var array_t1 = new Uint16Array(256 * 256);
var array_t2 = new Uint16Array(256 * 256);
var array_t2s = new Uint16Array(256 * 256);
var array_na_mm = new Uint16Array(256 * 256);
var array_na_t1 = new Uint16Array(256 * 256);
var array_na_ex_frac = new Uint16Array(256 * 256);
var array_na_t2s = new Uint16Array(256 * 256);
var array_na_t2f = new Uint16Array(256 * 256);
var k_data_im_re, k_result;

async function loadDataSet(path) {
    array_pd = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function () {
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            //console.log("pd", mm[0],mm[1], a.reduce((a,b)=>a+b)/a.length)
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0) * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
        xhr.onerror = reject
        xhr.open("GET", path+"/pd.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_t1 = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function () {
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            //console.log("t1", mm[0],mm[1], a.reduce((a,b)=>a+b)/a.length)
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0) * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
        xhr.onerror = reject
        xhr.open("GET", path+"/t1.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_t2 = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function () {
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            //console.log("t2", mm[0],mm[1], a.reduce((a,b)=>a+b)/a.length)
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0) * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
        xhr.onerror = reject
        xhr.open("GET", path+"/t2.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_t2s = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function () {
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            //console.log("t2", mm[0],mm[1], a.reduce((a,b)=>a+b)/a.length)
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0) * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
        xhr.onerror = reject
        xhr.open("GET", path+"/t2s.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_mm = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function (e) {
            if (e.target.status != 200) {
                resolve(null);
            }
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            //console.log("t2", mm[0],mm[1], a.reduce((a,b)=>a+b)/a.length)
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0) * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_mm.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_ex_frac = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function (e) {
            if (e.target.status != 200) {
                resolve(null);
            }
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            //console.log("t2", mm[0],mm[1], a.reduce((a,b)=>a+b)/a.length)
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0)// * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_ex_frac.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_t1 = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function (e) {
            if (e.target.status != 200) {
                resolve(null);
            }
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            //console.log("t2", mm[0],mm[1], a.reduce((a,b)=>a+b)/a.length)
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0) * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_t1.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_t2f = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function (e) {
            if (e.target.status != 200) {
                resolve(null);
            }
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            //console.log("t2", mm[0],mm[1], a.reduce((a,b)=>a+b)/a.length)
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0) * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_t2f.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_t2s = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = function (e) {
            if (e.target.status != 200) {
                resolve(null);
            }
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            //console.log("t2", mm[0],mm[1], a.reduce((a,b)=>a+b)/a.length)
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0) * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_t2s.bin.gz", true)
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    return [array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_t2f, array_na_t2s, zdim, ydim, xdim];
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
    for (var i = 0; i < k_data.length; i++) {
        if (fft_res[2 * i] == -1) {
            k_data[i] = -1;
        } else {
            k_data[i] = Math.sqrt(fft_res[2 * i] * fft_res[2 * i] + fft_res[2 * i + 1] * fft_res[2 * i + 1]);
            if (k_data[i] == -Infinity) {
                k_data[i] = 0;
            }
            maxval = Math.max(maxval, k_data[i]);
            minval = Math.min(minval, k_data[i]);
        }
    }

    for (var i = 0; i < (xdim / 2 + 1) * ydim; i++) {

        var val = (k_data[i] - minval) * 255 / maxval;

        var y = (Math.floor(i / (xdim / 2 + 1)) + ydim / 2) % ydim;
        var x = i % (xdim / 2 + 1) + xdim / 2;

        var j = y * xdim + x - 1;
        var k = y * xdim + (xdim - x);
        if (k_data[i] == -1) {
            k_result[j] = -1;
            k_result[k] = -1;
        } else {
            k_result[j] = val;
            k_result[k] = val;
        }
    }
    return k_result;
}

function genMapKSpace(xlines, ylines, fmin, fmax) {
    var dx = xdim / xlines;
    var dy = ydim / ylines;
    return function filterKSpace(value, d_index, array) {
        var index = Math.floor((d_index) / 2);
        var y = (Math.floor(index / (xdim/2 + 1))) % ydim;
        if (y > ydim / 2) {
            y = ydim - y;
        }
        var x = index % (xdim / 2 + 1);
        if (index > ydim * (xdim + 2) / 2) {
            x = x + xdim / 2 + 1;
        }
        var f = Math.sqrt((x) * (x) + (y) * (y));

        var y = (Math.floor(index / (xdim / 2 + 1)) + ydim / 2);
        var x = index % (xdim / 2 + 1) + xdim / 2;
        
        var res1 = f >= fmin && f <= fmax;
        var res = res1 && ((Math.floor(x / dx) * dx - x) < 1 && (Math.floor(x / dx) * dx - x) > -1) && ((Math.floor(y / dy) * dy - y) < 1 && (Math.floor(y / dy) * dy - y) > -1);


        return res ? value : 0;
    }
}

function genFilterKSpace(xlines, ylines, fmin, fmax) {
    var dx = xdim / xlines;
    var dy = ydim / ylines;
    return function filterKSpace(value, d_index, array) {
        //var index = Math.floor((d_index - 1) / 2);
        var index = d_index;
        var y = (Math.floor(index / (ydim / 2 + 1))) % ydim;
        if (y > ydim / 2) {
            y = ydim - y;
        }
        var x = index % (xdim / 2 + 1);
        if (index > ydim * (xdim + 2) / 2) {
            x = x + xdim / 2 + 1;
        }
        var res1 = Math.sqrt(x * x + y * y) >= fmin && Math.sqrt(x * x + y * y) <= fmax;
        var res = res1 && (Math.ceil(x / dx) * dx - x) < 1 && (Math.ceil(y / dy) * dy - y) < 1;
        return res ? true : false;
    }
}

function inverseKSpace(kSpace, xlines, ylines, fmin, fmax, noIfft = false) {
    //var filterKSpace = genFilterKSpace(xlines, ylines, 0, xdim * xdim)
    //console.log(xlines, ylines, fmin, fmax);
    //var mapKSpace = genMapKSpace(xlines, ylines, fmin, fmax)
    var result = new Float32Array(xdim * ydim * zdim);
    var k_result = new Float32Array(xdim * ydim * zdim);
    var zend = 10;
    for (var z = 0; z < zdim; z++) {
        var slice_data = kSpace.subarray(z * xdim * ydim * 2, (z + 1) * xdim * ydim * 2);

        var _fmax = fmax;
        if (z<zend) {
            _fmax = z/zend*z/zend * fmax;
        }
        if (z > zdim-zend) {
            _fmax = (zdim-z)/zend*(zdim-z)/zend * fmax;
        }
        var mapKSpace = genMapKSpace(xlines, ylines, fmin, _fmax)
        var input_data = slice_data.map(mapKSpace);
        k_result.set(transformKSpace(input_data), z * xdim * ydim);

        if (!noIfft) {
            var img_result = irfft2d(input_data, xdim, ydim)
            var maxval = 0;
            var minval = 999999999;
            for (var i = 0; i < img_result.length; i++) {
                maxval = Math.max(maxval, img_result[i]);
                minval = Math.min(minval, img_result[i]);
            }
            for (var i = 0; i < img_result.length; i++) {
                //result[i + z * xdim * ydim] = (img_result[i] - minval) / (maxval - minval)
                //result[i + 1 + z * xdim * ydim] = (img_result[i] - minval) / (maxval - minval)
                result[i + z * xdim * ydim] = img_result[i]
                result[i + 1 + z * xdim * ydim] = img_result[i]
            }
        }
    }
    if (noIfft) {
        return [undefined, k_result, [xlines, ylines, fmin, fmax]];
    }
    return [result, k_result, [xlines, ylines, fmin, fmax]];
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
        var val = pd * Math.exp(-te / t2) * (1 - Math.exp(-tr / t1));

        result[x] = Math.abs(val);
    }

    var k_result = calcKSpace(result);
    return [result, k_result];
}

function calcInversionRecovery(te, tr, ti) {
    var result = new Float32Array(array_t1.length);
    for (var x = 0; x < result.length; x++) {
        var t1 = array_t1[x]
        var t2 = array_t2[x]
        var pd = array_pd[x]
        if (t1 == 0) {
            t1 = 1.0;
        }
        var val = pd * (1.0 - 2.0 * Math.exp(-ti / t1) + Math.exp(-tr / ti)) * Math.exp(-te / t2);

        result[x] = Math.abs(val);
    }
    var k_result = calcKSpace(result);
    return [result, k_result];
}

const DELTAB = 19.5E-9;
const GYRO = 42.58;

function calcFlash(te, tr, fa) {
    var result = new Float32Array(array_t1.length);
    fa = fa * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    for (var x = 0; x < result.length; x++) {
        var t1 = array_t1[x]
        var t2 = array_t2[x]
        var pd = array_pd[x]
        if (t1 == 0) {
            t1 = 1.0;
        }
        var E1 = Math.exp(-tr/t1);
        var t2s = 1 / (1/t2 + GYRO*DELTAB);
        var val = pd * (1-E1)*sfa/((1-E1)*cfa) * Math.exp(-te/t2s)

        result[x] = Math.abs(val);
    }
    var k_result = calcKSpace(result);
    return [result, k_result];
}

function calcPSIF(te, tr, fa) {
    var result = new Float32Array(array_t1.length);
    fa = fa * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    var tfa = Math.tan(fa / 2);
    for (var x = 0; x < result.length; x++) {
        var t1 = array_t1[x]
        var t2 = array_t2[x]
        var pd = array_pd[x]
        if (t1 == 0) {
            t1 = 1.0;
        }
        var e1 = Math.exp(-tr/t1)
        var e2 = Math.exp(-tr/t2)
        var val = pd * sfa/(1+cfa)*(1-(1-e1*cfa)*Math.sqrt((1-e2*e2) /( (1-e1*cfa)*(1-e1*cfa)-e2*e2*(e1-cfa)*(e1-cfa) ))) * Math.exp(-te/t2);

        result[x] = Math.abs(val);
    }
    var k_result = calcKSpace(result);
    return [result, k_result];
}

function calcFISP(te, tr, fa) {
    var result = new Float32Array(array_t1.length);
    fa = fa * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    for (var x = 0; x < result.length; x++) {
        var t1 = array_t1[x]
        var t2 = array_t2[x]
        var t2s = array_t2s[x]
        var pd = array_pd[x]
        if (t1 == 0) {
            t1 = 1.0;
        }

        var e1 = Math.exp(-tr/t1)
        var e2 = Math.exp(-tr/t2)
        var val = pd * sfa/(1+cfa) * (1 - (e1-cfa)*Math.sqrt((1-e2*e2)/( (1-e1*cfa)*(1-e1*cfa)-e2*e2*(e1-cfa)*(e1-cfa) ) ) ) * Math.exp(-te/t2s);

        result[x] = Math.abs(val);
    }
    var k_result = calcKSpace(result);
    return [result, k_result];
}

function calcBalancedSSFP(te, tr, fa) {
    var result = new Float32Array(array_t1.length);
    fa = fa * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    for (var x = 0; x < result.length; x++) {
        var t1 = array_t1[x]
        var t2 = array_t2[x]
        var pd = array_pd[x]
        if (t1 == 0) {
            t1 = 1.0;
        }
        var e_tr_t1 = Math.exp(-tr/t1)
        var e_tr_t2 = Math.exp(-tr/t2)
        var val = pd * sfa * (1-Math.exp(-tr/t1)) / (1 - (e_tr_t1+e_tr_t2)*cfa - e_tr_t1*e_tr_t2 ) * Math.exp(-te/t2);

        result[x] = Math.abs(val);
    }
    var k_result = calcKSpace(result);
    return [result, k_result];
}

function calcSGRE(te, tr, fa) {
    var result = new Float32Array(array_t1.length);
    fa = fa * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    for (var x = 0; x < result.length; x++) {
        var t1 = array_t1[x]
        var t2s = array_t2s[x]
        var pd = array_pd[x]
        if (t1 == 0) {
            t1 = 1.0;
        }
        var E1 = Math.exp(-tr/t1);
        var val = pd * (1-E1)*sfa/(1-E1*cfa) * Math.exp(-te/t2s)

        result[x] = Math.abs(val);
    }
    var k_result = calcKSpace(result);
    return [result, k_result];
}

function downscale(arr, factor=2) {
    var r = new Float32Array(arr.length / (factor*factor));
    for(var z = 0; z<zdim; z++) {
        for(var y = 0;y<ydim; y+=factor) {
            for(var x=0;x<xdim;x+=factor) {
                opos = z*xdim*ydim + y*xdim + x;
                npos = z*xdim*ydim/(factor*factor) + y*xdim/(factor*factor) + x / factor;
                var sum = 0
                for(var i=0;i<factor;i++) {
                    for(var j=0;j<factor;j++) {
                        sum += arr[z*xdim*ydim + (y+i)*xdim + x + j];
                    }
                }
                r[npos] = sum / (factor*factor);
            }
        }
    }
    return r;
}

function upscale(arr, factor=2) {
    var r = new Float32Array(arr.length * (factor*factor));
    for(var z = 0; z<zdim; z++) {
        for(var y = 0;y<ydim; y+=factor) {
            for(var x=0;x<xdim;x+=factor) {
                opos = z*xdim*ydim + y*xdim + x;
                npos = z*xdim*ydim/(factor*factor) + y*xdim/(factor*factor) + x / factor;
                for(var i=0;i<factor;i++) {
                    for(var j=0;j<factor;j++) {
                        r[z*xdim*ydim + (y+i)*xdim + x + j] = arr[npos];
                    }
                }
            }
        }
    }
    return r;
}

function calcNa(te, tr, downscale=4) {
    //Na: s(t) = VbCb (1-exp(-tr/t1))[0.6*exp(-t/T2s) + 0.4*exp(-t/T2L)] + VfrCfrafr(1-exp(-tr/t1))exp(-t/T2fr) + n(t). 
    var result = new Float32Array(array_t1.length);
    //for (var x = 0; x < result.length; x++) {
    for(var z = 0; z<zdim; z++) {
        for(var y = Math.floor(downscale/2);y<ydim; y+=downscale) {
            for(var x = Math.floor(downscale/2);x<xdim;x+=downscale) {
                pos = z*xdim*ydim + y*xdim + x;
                var t2f = array_na_t2f[pos];
                var t2s = array_na_t2s[pos];
                var mm = array_na_mm[pos];
                var vol = array_na_ex_frac[pos];
                if (t2s == 0) {
                    t2s = 1;
                }
                if (t2f == 0) {
                    t2f = 1;
                }
                var val = (na_vol-vol)*mm * (1-Math.exp(-tr/na_t1)) * (0.6*Math.exp(-te/t2f) + 0.4*Math.exp(-te/t2s)) + vol*na_mm * (1-Math.exp(-tr/na_t1))*Math.exp(-te/na_t2fr)
                if (z == 10 && x == Math.floor(downscale/2) && y == Math.floor(downscale/2) ) {
                    console.log(vol, mm, val, x, y, z);
                }

                val = Math.abs(val) / 140;
                result[pos] = val;
                for(var yi =-Math.floor(downscale/2);yi<Math.floor(downscale/2); yi++) {
                    for(var xi =-Math.floor(downscale/2);xi<Math.floor(downscale/2);xi++) {
                        if (yi!=0 || xi!=0) {
                            pos = z*xdim*ydim + (y+yi)*xdim + x+xi;
                            result[pos] = val;
                        }
                    }
                }
            }
        }
    }
        
    var k_result = calcKSpace(result);

    //var r = inverseKSpace(k_data_im_re, 256, 256, 0, 64);
    return [result, k_result, [256, 256, 0, 256]];
}

function calcSQ(te_start, te_end, te_step, tau1=10, downscale=4) {
    //SQ: [mM].*(exp(-(TE_vec(kk)+t1)/T2s)+exp(-(TE_vec(kk)+t1)/T2f)) Simulation  von einer multi-echo Akquisition mit TE_vec=TE=[1,2,3,…]. 
    if (downscale < 1) { downscale = 1;}
    var result = new Float32Array(array_t1.length);
    var fa = 0.5*Math.PI;
    for (var x = 0; x < result.length; x++) {
        result[x] = 0;
    }
    var te_count = 0;
    for (var te=te_start; te<=te_end; te+=te_step) {
        te_count += 1;
        for(var z = 0; z<zdim; z++) {
            for(var y = Math.floor(downscale/2);y<ydim; y+=downscale) {
                for(var x = Math.floor(downscale/2);x<xdim;x+=downscale) {
                    pos = z*xdim*ydim + y*xdim + x;
                    if (pos > array_na_t2f.length) { continue; }
                    var t2f = array_na_t2f[pos];
                    var t2s = array_na_t2s[pos];
                    var mm = array_na_mm[pos];
                    if (t2s == 0) {
                        t2s = 1;
                    }
                    if (t2f == 0) {
                        t2f = 1;
                    }
                    var val = mm * ( Math.exp(-(te+tau1)/t2s) + Math.exp(-(te+tau1)/t2f) ) * Math.sin(fa);

                    //result[i] += Math.abs(val);
                    result[pos] += Math.abs(val);
                    for(var yi =-Math.floor(downscale/2);yi<Math.floor(downscale/2); yi++) {
                        for(var xi =-Math.floor(downscale/2);xi<Math.floor(downscale/2);xi++) {
                            if (yi!=0 || xi!=0) {
                                pos = z*xdim*ydim + (y+yi)*xdim + x+xi;
                                result[pos] += Math.abs(val);
                            }
                        }
                    }
                }
            }
        }
    }
    for (var x = 0; x < result.length; x++) {
        result[x] /= (te_count*140);
    }
    
    var k_result = calcKSpace(result);

    //var r = inverseKSpace(k_data_im_re, 256, 256, 0, 128);
    return [result, k_result, [256, 256, 0, 256]];
}

function calcTQ(te_start, te_end, te_step, tau1=10, tau2=0.1, downscale=10) {
    //TQ:  [mM].*((exp(-TE_vec(kk)/T2s)-exp(-TE_vec(kk)/T2f))*(exp(-t1/T2s)-exp(-t1/T2f))*exp(t2/T2s));
    // s = (e(-te/t2s)-e(-te/t2f))*(e(-τ1/t2s)-e(-τ1/t2f))*e(-τ2/t2s)  
    var result = new Float32Array(array_t1.length);
    for (var x = 0; x < result.length; x++) {
        result[x] = 0;
    }
    var te_count = 0;
    for (var te=te_start; te<=te_end; te+=te_step) {
        te_count += 1;
        for(var z = 0; z<zdim; z++) {
            for(var y = Math.floor(downscale/2);y<ydim; y+=downscale) {
                for(var x = Math.floor(downscale/2);x<xdim;x+=downscale) {
                    pos = z*xdim*ydim + y*xdim + x;
                    if (pos > array_na_t2f.length) { continue; }
                    var t2f = array_na_t2f[pos];
                    var t2s = array_na_t2s[pos];
                    var mm = array_na_mm[pos];
                    if (t2s == 0) {
                        t2s = 1;
                    }
                    if (t2f == 0) {
                        t2f = 1;
                    }
                    var val = mm * ( (Math.exp(-te/t2s) - Math.exp(-te/t2f)) * (Math.exp(-tau1/t2s)-Math.exp(-tau1/t2f)) * Math.exp(-tau2/t2s) );

                    //result[i] += Math.abs(val);
                    result[pos] += Math.abs(val);
                    for(var yi =-Math.floor(downscale/2);yi<Math.floor(downscale/2); yi++) {
                        for(var xi =-Math.floor(downscale/2);xi<Math.floor(downscale/2);xi++) {
                            if (yi!=0 || xi!=0) {
                                pos = z*xdim*ydim + (y+yi)*xdim + x+xi;
                                result[pos] += Math.abs(val);
                            }
                        }
                    }
                }
            }
        }
    }
    for (var x = 0; x < result.length; x++) {
        result[x] /= (te_count*140);
    }

    var k_result = calcKSpace(result);

    //var r = inverseKSpace(k_data_im_re, 256, 256, 0, 128);
    return [result, k_result, [256, 256, 0, 256]];
}

function calcTQF(te,tau1=10e-3,tau2=0.1e-3,fa=90,downscale=4) { // https://doi.org/10.1002/nbm.1548
    var s = new Float32Array(array_t1.length);
    fa = fa * Math.PI / 180;
    //for(var x=0;x<s.length;x++) {
    for(var z = 0; z<zdim; z++) {
        for(var y = Math.floor(downscale/2);y<ydim; y+=downscale) {
            for(var x = Math.floor(downscale/2);x<xdim;x+=downscale) {
                pos = z*xdim*ydim + y*xdim + x;
                var mm = array_na_mm[pos];
                var t2s = array_na_t2s[pos];
                var t2f = array_na_t2f[pos];
                val = Math.abs(mm * (Math.exp(-te/t2s)-Math.exp(-te/t2f)) * (Math.exp(-tau1/t2s)-Math.exp(-tau1/t2f)) * Math.exp(-tau2/t2s) * (Math.sin(fa)**5));
                result[pos] = val;
                for(var yi =-Math.floor(downscale/2);yi<Math.floor(downscale/2); yi++) {
                    for(var xi =-Math.floor(downscale/2);xi<Math.floor(downscale/2);xi++) {
                        if (yi!=0 || xi!=0) {
                            pos = z*xdim*ydim + (y+yi)*xdim + x+xi;
                            result[pos] = val;
                        }
                    }
                }
            }
        }
    }
    
    var k_result = calcKSpace(s);

    var r = inverseKSpace(k_data_im_re, 256, 256, 0, 64);
    return [r[0], k_result, [256, 256, 0, 64]];
}

function calcTQ_SQ_Ratio(tq_params, sq_params) {
    var tq = calcTQ(...tq_params)[0];
    var sq = calcSQ(...sq_params)[0];
    var result = new Float32Array(array_t1.length);
    for (var x = 0; x < result.length; x++) {
        result[x] = tq[x]/sq[x];
    }
    var k_result = calcKSpace(result);

    return [result, k_result, [256, 256, 0, 256]];
}

var queryableFunctions = {
    spinEcho: function (te, tr) {
        reply('result', calcSpinEcho(te, tr));
    },
    inversionRecovery: function (te, tr, ti) {
        reply('result', calcInversionRecovery(te, tr, ti));
    },
    flash: function(te, tr, fa) {
        reply('result', calcFlash(te, tr, fa));
    },
    bSSFP: function(te, tr, fa) {
        reply('result', calcBalancedSSFP(te, tr, fa));
    },
    psif: function(te, tr, fa) {
        reply('result', calcPSIF(te, tr, fa));
    },
    fisp: function(te, tr, fa) {
        reply('result', calcFISP(te, tr, fa));
    },
    sgre: function(te, tr, fa) {
        reply('result', calcSGRE(te, tr, fa));
    },
    na: function(te, tr, spacing) {
        reply('result', calcNa(te, tr, spacing));
    },
    sq: function(te_start, te_end, te_step, tau1, spacing) {
        reply('result', calcSQ(te_start, te_end, te_step, tau1, spacing));
    },
    tq: function(te_start, te_end, te_step, tau1, tau2, spacing) {
        reply('result', calcTQ(te_start, te_end, te_step, tau1, tau2, spacing));
    },
    tqf: function(te, tau1, tau2, fa, spacing) {
        reply('result', calcTQF(te, tau1, tau2, fa, spacing));
    },
    tqsqr: function(tq_params, sq_params) {
        reply('result', calcTQ_SQ_Ratio(tq_params, sq_params));
    },
    reco: function (xlines, ylines, fmin, fmax, noIfft) {
        reply('result', inverseKSpace(k_data_im_re, xlines, ylines, fmin, fmax, noIfft));
    },
    filterKSpace: function (xlines, ylines, fmin, fmax) {
        reply('result', inverseKSpace(k_data_im_re, xlines, ylines, fmin, fmax, true));
    },
    loadData: async function (path) {
        reply('loadData', await loadDataSet(path));
    },
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