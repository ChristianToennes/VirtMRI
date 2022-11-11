self.importScripts("./kissfft.js")
self.importScripts("./pako.min.js")
//self.importScripts("./discrete-wavelets.umd.min.js")
self.importScripts("./Mask_CS_Accel2.txt.js");
self.importScripts("./Mask_CS_Accel4.txt.js");

var k_xdim = 256;
var k_ydim = 256;
var k_zdim = 256;

const na_tau1 = 50;
const na_tau2 = 0.1;
const na_t1 = 39.2;
const na_t2s = 50;
const na_t2f = 4;
const na_t2fr = 60;
const na_mm = 140;
const na_vol = 0.7;

var ds = undefined;

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

class MRImage {
    constructor(xdim, ydim, zdim, kspace_filt, data, k_data, params, raw_data) {
        this.xdim = xdim;
        this.ydim = ydim;
        this.zdim = zdim;
        this.data = data;
        this.kSpace = k_data;
        this.kspace_filt = kspace_filt;
        this.params = params;
        this.raw_data = raw_data;
    }
}

function onloadDataSet(xhr,resolve) {
    return function(e) {
        if (e.target.status != 200) {
            resolve(null);
        } else {
            var resp = pako.inflate(xhr.response).buffer;
            var mm = new Float32Array(resp, 0, 2);
            var shape = new Uint16Array(resp, 8, 3);
            zdim = shape[0];
            ydim = shape[1];
            xdim = shape[2];
            var a = new Uint8Array(resp, 14);
            var b = new Float32Array(a.length);
            for (var x = 0; x < a.length; x++) {
                b[x] = (a[x] / 255.0) * (mm[1] - mm[0]) + mm[0];
            }
            resolve(b);
        }
    };
}

async function loadDataSet(path) {
    array_pd = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = onloadDataSet(xhr,resolve);
        xhr.onerror = reject;
        xhr.open("GET", path+"/pd.bin.gz", true);
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_t1 = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = onloadDataSet(xhr,resolve);
        xhr.onerror = reject;
        xhr.open("GET", path+"/t1.bin.gz", true);
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_t2 = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = onloadDataSet(xhr,resolve);
        xhr.onerror = reject;
        xhr.open("GET", path+"/t2.bin.gz", true);
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_t2s = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = onloadDataSet(xhr,resolve);
        xhr.onerror = reject;
        xhr.open("GET", path+"/t2s.bin.gz", true);
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_mm = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = onloadDataSet(xhr,resolve);
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_mm.bin.gz", true);
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_ex_frac = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = onloadDataSet(xhr,resolve);
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_ex_frac.bin.gz", true);
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_t1 = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = onloadDataSet(xhr,resolve);
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_t1.bin.gz", true);
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_t2f = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = onloadDataSet(xhr,resolve);
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_t2f.bin.gz", true);
        xhr.responseType = "arraybuffer";
        xhr.send()
    });
    array_na_t2s = await new Promise((resolve, reject) => {
        const xhr = new XMLHttpRequest();
        xhr.onload = onloadDataSet(xhr,resolve);
        xhr.onerror = reject;
        xhr.open("GET", path+"/na_t2s.bin.gz", true);
        xhr.responseType = "arraybuffer";
        xhr.send();
    });
    if (ds == undefined) {
        ds = make_dataset(array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_ex_frac, array_na_t2s, array_na_t2f);
    } else {
        free_dataset(ds);
        ds = make_dataset(array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_ex_frac, array_na_t2s, array_na_t2f);
    }
    return [array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_t2f, array_na_t2s, zdim, ydim, xdim];
}

function my_fft_stride(data, offset, stride, length) {
    var slice = new Float32Array(2*length);
    for(var i=0;i<slice.length;i+=2) {
        slice[i] = data[offset+i*stride];
        slice[i+1] = data[offset+i*stride+1];
        if(slice[i] == NaN || slice[i] == undefined) {
            console.log("nan before", i, offset, i*stride);
        }
    }
    //console.log("data", data.reduce((p,v) => Math.min(p,v)), data.reduce((p,v) => Math.max(p,v)))
    //console.log("input", slice.reduce((p,v) => Math.min(p,v)), slice.reduce((p,v) => Math.max(p,v)))
    var result = fft(slice, length);
    //console.log(length, slice.length, result.length);
    //console.log("result", result.reduce((p,v) => Math.min(p,v)), result.reduce((p,v) => Math.max(p,v)));
    for(var i=0;i<slice.length;i+=2) {
        data[offset+i*stride] = result[i];
        data[offset+i*stride+1] = result[i+1];
        if(data[offset+i*stride] == undefined) {
            console.log("nan after", i, offset+i*stride-data.length);
        }
    }
}

function my_ifft_stride(data, offset, stride, length) {
    var slice = new Float32Array(2*length);
    for(var i=0;i<slice.length;i+=2) {
        slice[i] = data[offset+i*stride];
        slice[i+1] = data[offset+i*stride+1];
        if(slice[i] == NaN || slice[i] == undefined) {
            console.log("nan before", i, offset, i*stride);
        }
    }
    //console.log("data", data.reduce((p,v) => Math.min(p,v)), data.reduce((p,v) => Math.max(p,v)))
    //console.log("input", slice.reduce((p,v) => Math.min(p,v)), slice.reduce((p,v) => Math.max(p,v)))
    var result = ifft(slice, length);
    //console.log(length, slice.length, result.length);
    //console.log("result", result.reduce((p,v) => Math.min(p,v)), result.reduce((p,v) => Math.max(p,v)));
    for(var i=0;i<slice.length;i+=2) {
        data[offset+i*stride] = result[i];
        data[offset+i*stride+1] = result[i+1];
        if(data[offset+i*stride] == undefined) {
            console.log("nan after", i, offset+i*stride-data.length);
        }
    }
}

function fft2d(data, xdim, ydim) {
    return kfft2d(data, xdim, ydim);
}

function ifft2d(data, xdim, ydim) {
    return kifft2d(data, xdim, ydim);
}

function fft3d(data, xdim, ydim, zdim) {
    return kfft3d(data, xdim, ydim, zdim);
}

function ifft3d(data, xdim, ydim, zdim) {
    return kifft3d(data, xdim, ydim, zdim);
}

function calcKSpace(result) {
    var xdim = result.xdim;
    var ydim = result.ydim;
    var zdim = result.zdim;
    var xoff = Math.floor((xdim-xdim) / 2);
    var yoff = Math.floor((ydim-ydim) / 2);
    xoff = 0;
    yoff = 0;
    k_data_im_re = new Float32Array(xdim * ydim * zdim * 2);
    k_data_im_re.fill(0);
    for (var z = 0; z < zdim; z++) {
        var slice_data = new Float32Array(xdim * ydim * 2);
        slice_data.fill(0);
        for (var x = 0; x < xdim; x++) {
            for (var y=0; y < ydim; y++) {
                slice_data[2*(xoff+x + (yoff+y)*xdim)] = result.data[x + y*xdim + z*xdim*ydim];
            }
        }
        var fft_res = fft2d(slice_data, xdim, ydim);
        k_data_im_re.set(fft_res, z * xdim * ydim * 2)
    }
    //console.log("calc", k_data_im_re.length, k_data_im_re.reduce((p,v)=>v!=v?p+1:p, 0))
    result.kSpace = transformKSpace3d(k_data_im_re, false, xdim, ydim, zdim);
    return result.kSpace;
}

function calcKSpace3d(result) {
    var xdim = result.xdim;
    var ydim = result.ydim;
    var zdim = result.zdim;
    var c_data = new Float32Array(xdim * ydim * zdim * 2);
    c_data.fill(0);
    for (var x = 0; x < xdim; x++) {
        for (var y=0; y < ydim; y++) {
            for(var z=0; z< zdim; z++) {
                c_data[2*(x + y*xdim + z*xdim*ydim)] = result.data[x + y*xdim + z*xdim*ydim];
            }
        }
    }
    k_data_im_re = fft3d(c_data, xdim, ydim, zdim);
    result.kSpace = transformKSpace3d(k_data_im_re, true, xdim, ydim, zdim);
    return result.kSpace;
}

function transformKSpace3d(k_data_im_re, fft3d, xdim, ydim, zdim) {
    var k_result = new Float32Array(xdim * ydim * zdim);
    try {
    for (var z = 0; z < zdim; z++) {
        var slice_data = new Float32Array(xdim * ydim*2);
        for (var x = 0; x < xdim; x++) {
            for (var y=0; y < ydim; y++) {
                slice_data[2*(x + y*xdim)] = k_data_im_re[2*(x + y*xdim + z*xdim*ydim)];
                slice_data[2*(x + y*xdim)+1] = k_data_im_re[2*(x + y*xdim + z*xdim*ydim) + 1];
                if(slice_data[2*(x + y*xdim)] == NaN || slice_data[2*(x + y*xdim)] == undefined || slice_data[2*(x + y*xdim)+1] == undefined || slice_data[2*(x + y*xdim)+1] == NaN) {
                    console.log("fft output undefined", x,y,z);
                }
                if(slice_data[2*(x + y*xdim)] == Infinity || slice_data[2*(x + y*xdim)] == -Infinity) {
                    slice_data[2*(x + y*xdim)] = 0;
                }
                if(slice_data[2*(x + y*xdim)+1] == Infinity || slice_data[2*(x + y*xdim)+1] == -Infinity) {
                    slice_data[2*(x + y*xdim)+1] = 0;
                }
            }
        }
        if(fft3d) {
            k_result.set(transformKSpace(slice_data, xdim, ydim), ((z+zdim/2)%zdim) * xdim * ydim);
            //k_result.set(transformKSpace(slice_data, xdim, ydim), z * xdim * ydim);
        } else {
            k_result.set(transformKSpace(slice_data, xdim, ydim), z * xdim * ydim);
        }
    } } catch (e) {
        console.log(e);
    }
    var max = k_result.reduce((p,v,i) => Math.max(p,v));
    var min = k_result.reduce((p,v) => Math.min(p,v));
    for(var i=0;i<k_result.length;i++) {
        k_result[i] = 255*(k_result[i]-min) / (max-min);
    }
    return k_result
}

function transformKSpace(fft_res, xdim, ydim) {
    var k_data = new Float32Array(xdim * ydim);
    var k_result = new Float32Array(xdim * ydim);
    var maxval = 0;
    var minval = 999999999;
    for (var i = 0; i < k_data.length; i++) {
        k_data[i] = Math.sqrt(fft_res[2 * i] * fft_res[2 * i] + fft_res[2 * i + 1] * fft_res[2 * i + 1]);
        maxval = Math.max(maxval, k_data[i]);
        minval = Math.min(minval, k_data[i]);
    }
    for(var x = 0;x<xdim;x++) {
        for(var y=0;y<ydim;y++) {
            var i = x + y*xdim;
            var val = (k_data[i] - minval) * 255 / (maxval-minval);

            var j = ((x+xdim/2)%xdim)+((y+ydim/2)%ydim)*xdim;
            k_result[j] = k_data[i];
        }
    }

    return k_result;
}

function genMapKSpace(xdim, ydim, xlines, ylines, fmin, fmax) {
    var dx = xdim / xlines;
    var dy = ydim / ylines;
    return function filterKSpace(value, d_index) {
        var index = Math.floor(d_index / 2);
        var y = Math.floor(index / xdim);
        var x = index-y*xdim;
        if(y>ydim/2) {
            y = ydim-y;
        }
        if(x>xdim/2) {
            x = xdim-x;
        }
        var f = Math.sqrt(x*x+y*y);
        
        var res1 = f >= fmin && f <= fmax;
        var res = res1 && ((Math.floor(x / dx) * dx - x) < 1 && (Math.floor(x / dx) * dx - x) > -1) && ((Math.floor(y / dy) * dy - y) < 1 && (Math.floor(y / dy) * dy - y) > -1);


        return res ? value : 0;
    }
}

function genMapKSpace3d(xdim, ydim, xlines, ylines, fmin, fmax) {
    var dx = xdim / xlines;
    var dy = ydim / ylines;
    return function filterKSpace(value, d_index) {
        var index = Math.floor(d_index / 2);
        var z = Math.floor(index / (xdim*ydim));
        var y = Math.floor((index-z*xdim*ydim) / xdim);
        var x = index-y*xdim-z*xdim*ydim;
        if(y>ydim/2) {
            y = ydim-y;
        }
        if(x>xdim/2) {
            x = xdim-x;
        }
        var f = Math.sqrt(x*x+y*y);
        
        var res1 = f >= fmin && f <= fmax;
        var res = res1 && ((Math.floor(x / dx) * dx - x) < 1 && (Math.floor(x / dx) * dx - x) > -1) && ((Math.floor(y / dy) * dy - y) < 1 && (Math.floor(y / dy) * dy - y) > -1);


        return res ? value : 0;
    }
}

function inverseKSpace(kSpace, xdim, ydim, zdim, xlines, ylines, fmin, fmax, noIfft = false) {
    //var mapKSpace = genMapKSpace(xlines, ylines, fmin, fmax)
    var result = new Float32Array(xdim * ydim * zdim);
    var k_result = new Float32Array(xdim * ydim * zdim);
    k_result.fill(0);
    var zend = Math.max(Math.floor(zdim/25), 1);
    for (var z = 0; z < zdim; z++) {
        if (2*z*xdim*ydim > kSpace.length) { zdim = z; break; }
        var slice_data = new Float32Array(xdim*ydim*2);
        for(var i=0;i<slice_data.length;i++) {
            slice_data[i] = kSpace[i+z*xdim*ydim*2];
        }

        var _fmax = fmax;
        if (z<zend) {
            _fmax = z/zend*z/zend * fmax;
        }
        if (z > zdim-zend) {
            _fmax = (zdim-z)/zend*(zdim-z)/zend * fmax;
        }
        var mapKSpace = genMapKSpace(xdim, ydim, xlines, ylines, fmin, _fmax)
        var input_data = slice_data.map(mapKSpace);
        k_result.set(transformKSpace(input_data, xdim, ydim), z * xdim * ydim);

        if (!noIfft) {
            var fft_result = ifft2d(input_data, xdim, ydim);

            var xoff = Math.floor((xdim-xdim) / 2);
            var yoff = Math.floor((ydim-ydim) / 2);
            for (var x = 0; x < xdim; x++) {
                for (var y=0; y < ydim; y++) {
                    result[x+y*xdim+z*xdim*ydim] = Math.sqrt(fft_result[2*(xoff+x + (yoff+y)*xdim)]*fft_result[2*(xoff+x + (yoff+y)*xdim)]+fft_result[2*(xoff+x + (yoff+y)*xdim)+1]*fft_result[2*(xoff+x + (yoff+y)*xdim)+1]);
                }
            }
        }
    }
    if (noIfft) {
        return [undefined, k_result, [xlines, ylines, fmin, fmax]];
    }
    //result = new MRImage(xdim, ydim, zdim, result);
    return [result, k_result, [xlines, ylines, fmin, fmax]];
}

function inverseKSpace3d(kSpace, xdim, ydim, zdim, xlines, ylines, fmin, fmax, noIfft = false) {
    var result = new Float32Array(xdim * ydim * zdim);
    var zend = Math.max(Math.floor(zdim/25), 1);
    
    var _fmax = fmax;
    if (z<zend) {
        _fmax = z/zend*z/zend * fmax;
    }
    if (z > zdim-zend) {
        _fmax = (zdim-z)/zend*(zdim-z)/zend * fmax;
    }
    var mapKSpace = genMapKSpace3d(xdim, ydim, xlines, ylines, fmin, fmax)
    var input_data = kSpace.map(mapKSpace);
    var k_result = transformKSpace3d(kSpace, true, xdim, ydim, zdim);

    if (!noIfft) {
        var fft_result = ifft3d(kSpace, xdim, ydim, zdim);
        var xoff = Math.floor((xdim-xdim) / 2);
        var yoff = Math.floor((ydim-ydim) / 2);
        for (var x = 0; x < xdim; x++) {
            for (var y=0; y < ydim; y++) {
                for(var z=0;z<zdim;z++) {
                    var pos = 2*(xoff+x + (yoff+y)*xdim+z*xdim*k_ydim);
                    result[x+y*xdim+z*xdim*ydim] = Math.sqrt(fft_result[pos]*fft_result[pos]+fft_result[pos+1]*fft_result[pos+1]);
                }
            }
        }
    }
    if (noIfft) {
        return [undefined, k_result, [xlines, ylines, fmin, fmax]];
    }
    //result = new MRImage(xdim, ydim, zdim, result);
    return [result, k_result, [xlines, ylines, fmin, fmax]];
}

var call_count = 0;
function calcSpinEcho(pos, params) {
    call_count += 1;
    var te = params["te"];
    var tr = params["tr"];
    //var fa = params["fa"] * Math.PI / 180;
    var t1 = array_t1[pos]
    var t2 = array_t2[pos]
    var pd = array_pd[pos]
    if (t1 == 0) { t1 = 1.0; }
    if (t2 == 0) { t2 = 1.0; }
    var val = pd * Math.exp(-te / t2) * (1 - Math.exp(-tr / t1));
    val = Math.abs(val);
    return val;
}

function calcInversionRecovery(pos, params) {
    call_count += 1;
    var te = params["te"];
    var tr = params["tr"];
    var ti = params["ti"];
    //var fa = params["fa"] * Math.PI / 180;
    var t1 = array_t1[pos]
    var t2 = array_t2[pos]
    var pd = array_pd[pos]
    if (t1 == 0) { t1 = 1.0; }
    if (t2 == 0) { t2 = 1.0; }    
    var val = pd * (1.0 - 2.0 * Math.exp(-ti / t1) + Math.exp(-tr / ti)) * Math.exp(-te / t2);
    val = Math.abs(val);
    return val;
}

const DELTAB = 19.5E-9;
const GYRO = 42.58;

function calcFlash(pos, params) {
    var te = params["te"];
    var tr = params["tr"];
    var fa = params["fa"] * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    var t1 = array_t1[pos]
    var t2 = array_t2[pos]
    var pd = array_pd[pos]
    if (t1 == 0) { t1 = 1.0; }
    if (t2 == 0) { t2 = 1.0; }
    var E1 = Math.exp(-tr/t1);
    var t2s = 1 / (1/t2 + GYRO*DELTAB);
    var val = pd * (1-E1)*sfa/((1-E1)*cfa) * Math.exp(-te/t2s)
    val = Math.abs(val);
    return val;
}

function calcPSIF(pos, params) {
    var te = params["te"];
    var tr = params["tr"];
    var fa = params["fa"] * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    var t1 = array_t1[pos]
    var t2 = array_t2[pos]
    var pd = array_pd[pos]
    if (t1 == 0) { t1 = 1.0; }
    if (t2 == 0) { t2 = 1.0; }
    var e1 = Math.exp(-tr/t1)
    var e2 = Math.exp(-tr/t2)
    var val = pd * sfa/(1+cfa)*(1-(1-e1*cfa)*Math.sqrt((1-e2*e2) /( (1-e1*cfa)*(1-e1*cfa)-e2*e2*(e1-cfa)*(e1-cfa) ))) * Math.exp(-te/t2);
    val = Math.abs(val);
    return val;
}

function calcFISP(pos, params) {
    var te = params["te"];
    var tr = params["tr"];
    var fa = params["fa"] * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    var t1 = array_t1[pos]
    var t2 = array_t2[pos]
    var t2s = array_t2s[pos]
    var pd = array_pd[pos]
    if (t1 == 0) { t1 = 1.0; }
    if (t2 == 0) { t2 = 1.0; }
    var e1 = Math.exp(-tr/t1)
    var e2 = Math.exp(-tr/t2)
    var val = pd * sfa/(1+cfa) * (1 - (e1-cfa)*Math.sqrt((1-e2*e2)/( (1-e1*cfa)*(1-e1*cfa)-e2*e2*(e1-cfa)*(e1-cfa) ) ) ) * Math.exp(-te/t2s);
    val = Math.abs(val);
    return val;
}

function calcBalancedSSFP(pos, params) {
    var te = params["te"];
    var tr = params["tr"];
    var fa = params["fa"] * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    var t1 = array_t1[pos]
    var t2 = array_t2[pos]
    var pd = array_pd[pos]
    if (t1 == 0) { t1 = 1.0; }
    if (t2 == 0) { t2 = 1.0; }
    var e_tr_t1 = Math.exp(-tr/t1)
    var e_tr_t2 = Math.exp(-tr/t2)
    var val = pd * sfa * (1-Math.exp(-tr/t1)) / (1 - (e_tr_t1+e_tr_t2)*cfa - e_tr_t1*e_tr_t2 ) * Math.exp(-te/t2);
    val = Math.abs(val);
    return val;
}

function calcSGRE(pos, params) {
    var te = params["te"];
    var tr = params["tr"];
    var fa = params["fa"] * Math.PI / 180;
    var sfa = Math.sin(fa);
    var cfa = Math.cos(fa);
    
    var t1 = array_t1[pos]
    var t2s = array_t2s[pos]
    var pd = array_pd[pos]
    if (t1 == 0) {
        t1 = 1.0;
    }
    var E1 = Math.exp(-tr/t1);
    var val = pd * (1-E1)*sfa/(1-E1*cfa) * Math.exp(-te/t2s)
    val = Math.abs(val);
    return val;
}

function calcNa(pos, params) {
    call_count += 1;
    var te = params["te"];
    var tr = params["tr"];
    //Na: s(t) = VbCb (1-exp(-tr/t1))[0.6*exp(-t/T2s) + 0.4*exp(-t/T2L)] + VfrCfrafr(1-exp(-tr/t1))exp(-t/T2fr) + n(t). 
    var t1 = array_na_t1[pos];
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
    var val = (na_vol-vol)*mm * (1-Math.exp(-tr/t1)) * (0.6*Math.exp(-te/t2f) + 0.4*Math.exp(-te/t2s)) + vol*na_mm * (1-Math.exp(-tr/na_t1))*Math.exp(-te/na_t2fr)
    val = Math.abs(val) / 140;
    return val;
}

function calcSQ(pos, params) {
    var te_start = params["te_start"];
    var te_end = params["te_end"];
    var te_step = params["te_step"];
    var tau1 = "tau1" in params ? params["tau1"] : 10;
    var fa = ("fa" in params ? params["fa"] : 90) * Math.PI / 180;
    //SQ: [mM].*(exp(-(TE_vec(kk)+t1)/T2s)+exp(-(TE_vec(kk)+t1)/T2f)) Simulation  von einer multi-echo Akquisition mit TE_vec=TE=[1,2,3,…]. 
    var result = 0;
    var te_count = 0;
    for (var te=te_start; te<=te_end; te+=te_step) {
        te_count += 1;
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
        result += Math.abs(val);
    }
    result /= (te_count*140);
    return result
}

function calcTQ(pos, params) {
    var te_start = params["te_start"];
    var te_end = params["te_end"];
    var te_step = params["te_step"];
    var tau1 = "tau1" in params ? params["tau1"] : 10;
    var tau2 = "tau2" in params ? params["tau2"] : 0.1;
    var result = 0;
    var te_count = 0;
    for (var te=te_start; te<=te_end; te+=te_step) {
        te_count += 1;
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

        result += Math.abs(val);
    }
    
    result /= (te_count*140);
    
    return result;
}

function calcTQF(pos, params) { // https://doi.org/10.1002/nbm.1548
    var te = params["te"]
    var tau1 = "tau1" in params ? params["tau1"] : 10e-3
    var tau2 = "tau2" in params ? params["tau2"] : 0.1e-3
    var fa = ("fa" in params ? params["fa"] : 90 ) * Math.PI / 180;
    var mm = array_na_mm[pos];
    var t2s = array_na_t2s[pos];
    var t2f = array_na_t2f[pos];
    val = Math.abs(mm * (Math.exp(-te/t2s)-Math.exp(-te/t2f)) * (Math.exp(-tau1/t2s)-Math.exp(-tau1/t2f)) * Math.exp(-tau2/t2s) * (Math.sin(fa)**5));
    return val;
}

function calcTQ_SQ_Ratio(pos, params) {
    var tq_params = params["tq_params"];
    var sq_params = params["sq_params"];
    var tq = calcTQ(pos, tq_params);
    var sq = calcSQ(pos, sq_params);
    var val = tq / sq;
    return val;
}

function sparsify(k_space) {
    return k_space;
}

function sparsify_transform(img) {
    return img;
}

function i_sparsify_transform(img) {
    return img;
}

function idwt2d(slice_data, wavelet, xdim, ydim) {
    //console.log(slice_data.length);
    var line = new Float32Array(xdim);
    var result = new Float32Array(slice_data.length);
    for(var y=0;y<ydim;y++) {
        for(var x=0;x<xdim;x++) {
            line[x] = slice_data[x+y*xdim];
        }
        //console.log(line.length, Array.from(line).length);
        result.set(wt.idwt(Array.from(line), wavelet, "zero"), y*xdim);
    }
    var tmp_result;
    for(var x=0;x<xdim;x++) {
        for(var y=0;y<ydim;y++) {
            line[y] = result[x+y*xdim];
        }
        tmp_result = wt.idwtc(line, wavelet);
        for(var y=0;y<ydim;y++) {
            result[x+y*xdim] = tmp_result[y];
        }
    }
    return result;
}

function dwt2d(slice_data, wavelet, xdim, ydim) {
    var line = new Float32Array(xdim);
    var result = new Float32Array(slice_data.length);
    for(var y=0;y<ydim;y++) {
        for(var x=0;x<xdim;x++) {
            line[x] = slice_data[x+y*xdim];
        }
        result.set(wt.waverec(line, wavelet), y*xdim);
    }
    var tmp_result;
    for(var x=0;x<xdim;x++) {
        for(var y=0;y<ydim;y++) {
            line[y] = result[x+y*xdim];
        }
        tmp_result = wt.waverec(line, wavelet);
        for(var y=0;y<ydim;y++) {
            result[x+y*xdim] = tmp_result[y];
        }
    }
    return result;
}

function compressed_sensing_mriquestions(k_space, params) {
    var threshold = "cs_threshold" in params ? params["cs_threshold"] : 0.3;
    k_space = sparsify(k_space);
    var initial_img = idwt2d(k_space, "db1", xdim, ydim);
    var update_img = sparsify_transform(initial_img);
    for(var i=0;i<update_img.length;i++) {
        if(update_img[i] < threshold) { 
            update_img[i] = 0;
        }
    }
    update_img = i_sparsify_transform(update_img);
    var k_space_d = dwt2d(update_img, "db1", xdim, ydim);
    for(var i=0;i<k_space_d.length;i++) {
        k_space_d[i] -= k_space[i];
    }
    for(var i=0;i<initial_img.length;i++) {
        initial_img[i] += update_img[i];
    }
    return initial_img;
}

function filter_kspace_cs(k_data_im_re, mask_cs, zdim) {
    var k_data_im_re_f = new Float32Array(k_data_im_re.length);
    for(var z=0;z<zdim;z++) {
        for(var y=0;y<ydim;y++) {
            for(var x=0;x<xdim;x++) {
                var pos = 2*(x+y*xdim+z*xdim*ydim);
                if(mask_cs[y][z%100] == 1) {
                    k_data_im_re_f[pos] = k_data_im_re[pos]; 
                    k_data_im_re_f[pos+1] = k_data_im_re[pos+1];
                } else {
                    k_data_im_re_f[pos] = 0
                    k_data_im_re_f[pos+1] = 0
                }
            }
        }
    }
    return k_data_im_re_f;
}

function print_stats(name, array) {
    var max = array.reduce((p,c) => Math.max(p,c), Math.max());
    var min = array.reduce((p,c) => Math.min(p,c), Math.min());
    var sum = array.reduce((p,c) => p+c, 0);
    var nans = array.reduce((p,c) => { if(c == NaN){return p+1;}else{return p;}}, 0);
    var ones = array.reduce((p,c) => { if(c != 0){return p+1;}else{return p;}}, 0);
    var zeros = array.reduce((p,c) => { if(c == 0){return p+1;}else{return p;}}, 0);
    //console.log(name, array.length, min, max, sum/array.length, nans, ones, zeros);
    console.log(name, min, max, sum);
}

function compressed_sensing(f_data, params) {
    //%%CS_RECON_CARTESIAN  Spatial-only 2D and 3D CS recon.
    //%
    //% Total variation (L1-norm of gradient) constrained compressed
    //% sensing reconstruction of undersampled Cartesian data.
    //% Samples not acquired are expected to be zeros in the arrays.
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
    var lambda = 1.0;
    var lambda2 = 0.3;
    var mu = "cs_mu" in params ? params["cs_mu"] : 1.0;
    var gam = 1;
    
    var data_ndims = zdim == 1 ? 2 : 3;
    //mask = data ~= 0.0;   //% not perfect, but good enough
    //print_stats("data", f_data);
    var f_mask = new Uint8Array(f_data.length);
    for(var i=0;i<f_data.length;i+=2) {
        f_mask[i] = (f_data[i]+f_data[i+1]) != 0;
        f_mask[i+1] = f_mask[i];
    }
    //print_stats("mask", f_mask);
    //% normalize the data so that standard parameter values work
    var norm_factor = get_norm_factor(k_ydim, f_data);
    for(var i=0;i<f_data.length;i++) {
        f_data[i] = norm_factor * f_data[i];
    }
    //print_stats("data", f_data);
    //% Reserve memory for the auxillary variables
    var f_data0 = new Float32Array(f_data.length);
    f_data0.set(f_data);
    var img = new Float32Array(f_data.length);
    img.fill(0);
    var X = new Float32Array(f_data.length*3);
    X.fill(0);
    var B = new Float32Array(f_data.length*3);
    B.fill(0);
    var i_murf = new Float32Array(f_data.length);
    var i_rhs = new Float32Array(f_data.length);
    //% Build Kernels
    var scale = Math.sqrt(f_data.length/2);
    var f_uker = new Float32Array(f_data.length);
    f_uker.fill(0);
    if (data_ndims == 2) {
      f_uker[0] = 4;
      f_uker[1*2] = -1;
      f_uker[xdim*2] = -1;
      f_uker[(xdim-1)*2] = -1;
      f_uker[xdim*(ydim-1)*2] = -1;
      //console.log("index", 0, 1, xdim, xdim-1, xdim*(ydim-1)); 
    } else {// data_ndims == 3
      f_uker[0] = 8;
      f_uker[1*2] = -1;
      f_uker[xdim*2] = -1;
      f_uker[xdim*ydim*2] = -1;
      f_uker[(xdim-1)*2] = -1;
      f_uker[xdim*(ydim-1)*2] = -1;
      f_uker[xdim*ydim*(zdim-1)*2] = -1;
    }
    //print_stats("uker in", f_uker);
    f_uker = data_ndims==2? fft2d(f_uker, xdim, ydim) : fft3d(f_uker, xdim, ydim, zdim);
    //print_stats("uker mid", f_uker);
    //console.log(mu, lambda, gam);
    var real, imag, d;
    for(var i=0;i<f_uker.length;i+=2) {
        real = (mu*f_mask[i]+lambda*f_uker[i]+gam);
        imag = lambda*f_uker[i+1];
        d = real*real+imag*imag;
        f_uker[i] = real/d;
        f_uker[i+1] = -imag/d;
    }
    //print_stats("uker out", f_uker);
    //%  Do the reconstruction
    for(var outer = 0;outer < nbreg; outer++) {
      for(var i=0;i<f_data.length;i++) {
        i_murf[i] = mu*f_mask[i]*f_data[i];
      }
      //print_stats("murf in", i_murf);
      i_murf = data_ndims==2? ifft2d(i_murf, xdim, ydim) : ifft3d(i_murf, xdim, ydim, zdim);
      //print_stats("murf out", i_murf);
      for(var inner = 0;inner < ninner;inner++) {
        //% update u
        for(var x=0;x<xdim;x++) {
            for(var y=0;y<ydim;y++) {
                for(var z=0;z<zdim;z++) {
                    for(var j=0;j<2;j++) {
                        var i = 2*3*(x + y*xdim + z*xdim*ydim) + j;
                        var d = x<xdim-1 ? X[i]-B[i] - X[2*3*(x+1 + y*xdim + z*xdim*ydim)+j]+B[2*3*(x+1 + y*xdim + z*xdim*ydim)+j] : X[i]-B[i] - X[2*3*(y*xdim + z*xdim*ydim)+j]+B[2*3*(y*xdim + z*xdim*ydim)+j];
                        d += y<ydim-1 ? X[i+2]-B[i+2] - X[2*3*(x + (y+1)*xdim + z*xdim*ydim)+j+2]+B[2*3*(x + (y+1)*xdim + z*xdim*ydim)+j+2] : X[i+2]-B[i+2] - X[2*3*(x + z*xdim*ydim)+j+2]+B[2*3*(x + z*xdim*ydim)+j+2];
                        d += z<zdim-1 ? X[i+4]-B[i+4] - X[2*3*(x + y*xdim + (z+1)*xdim*ydim)+j+4]+B[2*3*(x + y*xdim + (z+1)*xdim*ydim)+j+4] : X[i+4]-B[i+4] - X[2*3*(x + y*xdim) +j+4]+B[2*3*(x + y*xdim) + j+4];
                        i_rhs[2*(x + y*xdim + z*xdim*ydim) + j] = i_murf[2*(x + y*xdim + z*xdim*ydim) + j]*scale + lambda*d + gam*img[2*(x + y*xdim + z*xdim*ydim) + j];
                    }
                }
            }
        }
        //console.log("scale", scale);
        //print_stats("rhs in", i_rhs);
        var f_rhs = data_ndims==2? fft2d(i_rhs, xdim, ydim) : fft3d(i_rhs, xdim, ydim, zdim);
        //print_stats("rhs out", f_rhs);
        for(var i=0;i<f_rhs.length;i+=2) {
            var x = f_rhs[i]; var y = f_rhs[i+1]; var u = f_uker[i]; var v = f_uker[i+1];
            f_rhs[i] = x*u-y*v;
            f_rhs[i+1] = x*v+y*u;
        }
        //print_stats("img_in", f_rhs);
        img = data_ndims==2? ifft2d(f_rhs, xdim, ydim) : ifft3d(f_rhs, xdim, ydim, zdim);
        //print_stats("img", img);
        //% update x and y
        for(var x=0;x<xdim;x++) {
            for(var y=0;y<ydim;y++) {
                for(var z=0;z<zdim;z++) {
                    for(var j=0;j<2;j++) {
                        var i = 2*3*(x + y*xdim + z*xdim*ydim) + j;
                        B[i] = x>0 ? img[2*(x + y*xdim + z*xdim*ydim)+j]-img[2*(x-1 + y*xdim + z*xdim*ydim)+j] + B[i] : img[2*(x + y*xdim + z*xdim*ydim)+j]-img[2*(xdim-1 + y*xdim + z*xdim*ydim)+j] + B[i];
                        B[i+2] = y>0 ? img[2*(x + y*xdim + z*xdim*ydim)+j]-img[2*(x + (y-1)*xdim + z*xdim*ydim)+j] + B[i+2] : img[2*(x + y*xdim + z*xdim*ydim)+j]-img[2*(x + (ydim-1)*xdim + z*xdim*ydim)+j] + B[i+2];
                        B[i+4] = z>0 ? img[2*(x + y*xdim + z*xdim*ydim)+j]-img[2*(x + y*xdim + (z-1)*xdim*ydim)+j] + B[i+4] : img[2*(x + y*xdim + z*xdim*ydim)+j]-img[2*(x + y*xdim + (zdim-1)*xdim*ydim)+j] + B[i+4];
                    }
                }
            }
        }
        for(var i=0;i<X.length;i+=2) {
            var s = Math.sqrt(X[i]*X[i]+X[i+1]*X[i+1]);
            var ss = s - 1/lambda2;
            ss = ss > 0 ? ss : 0;
            s = s + (s < 1/lambda2 ? 1 : 0);
            ss = ss / s;
            X[i] = ss * X[i];
            //% update bregman parameters
            B[i] = B[i] - X[i];
            X[i+1] = ss * X[i+1];
            //% update bregman parameters
            B[i+1] = B[i+1] - X[i+1];
        }
      }
      var k_img = data_ndims==2? fft2d(img, xdim, ydim) : fft3d(img, xdim, ydim, zdim);
      for(var i=0;i<f_data.length;i++) {
        f_data[i] = f_data[i] + f_data0[i] - f_mask[i]*k_img[i] / scale;
      }
    }
    //% undo the normalization so that results are scaled properly
    //print_stats("img in", img);
    var result_img = new Float32Array(f_data.length/2);
    for(var i=0;i<result_img.length;i++) {
        result_img[i] = Math.sqrt(img[2*i]*img[2*i]+img[2*i+1]*img[2*i+1]) / norm_factor / scale;
    }
    //print_stats("img out", result_img);
    return result_img;
}

function get_norm_factor(ydim, data) {
    var norm = 0;
    for(var i=0;i<data.length;i++) {
        norm += data[i]*data[i];
    }
    norm = Math.sqrt(norm) / ydim;
    //print_stats("norm", norm);
    //console.log("norm", norm, 1/norm);
    return 1 / norm;
}

function dist(z0,y0,x0,z1,y1,x1) {
    var d = Math.sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))
    return [z1*xdim*ydim+y1*xdim+x1, d];
}

function getNearest(z,y,x) {
    return [[Math.round(z-0.5)*k_xdim*k_ydim+Math.round(y-0.5)*k_xdim+Math.round(x-0.5)], [1]];
    /*var [pos, d] = getNeighbours(z,y,x);
    var [max,mi] = d.reduce(function(o, v, i) { 
        if(o[0] < v) {
            return [v, i];
        } else { return o; }
    }, [d[0], 0]);
    return [[pos[mi]], [max]];*/
}

function getNeighbours(z,y,x) {
    var [dpos1, d1] = dist(z,y,x, Math.floor(z), Math.floor(y), Math.floor(x));
    var [dpos2, d2] = dist(z,y,x, Math.floor(z), Math.floor(y), Math.ceil(x));
    var [dpos3, d3] = dist(z,y,x, Math.floor(z), Math.ceil(y),  Math.floor(x));
    var [dpos4, d4] = dist(z,y,x, Math.floor(z), Math.ceil(y),  Math.ceil(x));
    var [dpos5, d5] = dist(z,y,x, Math.ceil(z),  Math.floor(y), Math.floor(x));
    var [dpos6, d6] = dist(z,y,x, Math.ceil(z),  Math.floor(y), Math.ceil(x));
    var [dpos7, d7] = dist(z,y,x, Math.ceil(z),  Math.ceil(y),  Math.floor(x));
    var [dpos8, d8] = dist(z,y,x, Math.ceil(z),  Math.ceil(y),  Math.ceil(x));
    var pos = [];
    var d = [];
    pos.push(dpos1);
    d.push(1-d1);
    if (d2 > 0 && !pos.includes(dpos2)) { pos.push(dpos2); d.push(1-d2); }
    if (d3 > 0 && !pos.includes(dpos3)) { pos.push(dpos3); d.push(1-d3); }
    if (d4 > 0 && !pos.includes(dpos4)) { pos.push(dpos4); d.push(1-d4); }
    if (d5 > 0 && !pos.includes(dpos5)) { pos.push(dpos5); d.push(1-d5); }
    if (d6 > 0 && !pos.includes(dpos6)) { pos.push(dpos6); d.push(1-d6); }
    if (d7 > 0 && !pos.includes(dpos7)) { pos.push(dpos7); d.push(1-d7); }
    if (d8 > 0 && !pos.includes(dpos8)) { pos.push(dpos8); d.push(1-d8); }
    return [pos, d];
}

function getVolumeVoxel(z,y,x, zdim, ydim, xdim) {
    var xs = k_xdim/xdim/2;
    var ys = k_ydim/ydim/2;
    var zs = k_zdim/zdim/2;

    var poss = [];
    var ds = [];
    x_start=Math.floor(x-xs);
    x_end=Math.ceil(x_start+2*xs);
    y_start=Math.floor(y-ys);
    y_end=Math.ceil(y_start+2*ys);
    z_start=Math.floor(z-zs);
    z_end=Math.ceil(z_start+2*zs);
    for(xi=x_start;xi<x_end;xi+=1) {
        for(yi=y_start;yi<y_end;yi+=1) {
            for(zi=z_start;zi<z_end;zi+=1) {
                var pos = zi*k_xdim*k_ydim+yi*k_xdim+xi;
                if(pos >= 0 && pos < k_xdim*k_ydim*k_zdim && xi>=0 && xi<k_xdim && yi>=0 && yi<k_ydim && zi>=0 && zi<k_zdim) {
                    poss.push(pos);
                    ds.push(1);
                }
            }
        }
    }
    /*if(poss.length == 0) {
        console.log("x", k_xdim, xdim, x, xs, Math.round(x-xs), Math.ceil(x+xs), x_start, x_end);
        console.log("y", k_ydim, ydim, y, ys, Math.round(y-ys), Math.ceil(y+ys), y_start, y_end);
        console.log("z", k_zdim, zdim, z, zs, Math.round(z-zs), Math.ceil(z+zs), z_start, z_end);
    }*/
    return [poss, ds];
}

function addGaussianNoiseRatio(img, mean, sigma) {
    var result = new Float32Array(img.length);
    var e1 = Math.sqrt(2.0/Math.E);
    var r1, r2;
    for(var i=0;i<result.length;i++) {
        {
            r1 = 2*Math.random()-1;
            r2 = 2*Math.random()-1;
            r2 = (2.0*r2-1.0)*e1;
        } while(r1== 0 && -4*r1*r1*Math.log(r1) < r2*r2);
        var n = sigma * r2 / r1 + mean;
        result[i] = img[i] + n;
    }
    return result
}

function addGaussianNoiseBoxMuller(img, mean, sigma) {
    var result = new Float32Array(img.length);
    var u,v,s,x,y;
    for(var i=0;i<result.length;i+=2) {
        u = Math.random();
        v = Math.random();
        s = Math.sqrt(-2.0*Math.log(u));
        x = s * Math.cos(2.0*Math.PI*v);
        y = s * Math.sin(2.0*Math.PI*v);
        result[i] = img[i] + (sigma*x+mean);
        result[i+1] = img[i+1] + (sigma*y+mean);
    }
    return result
}

function addImageNoise(img, params) {
    var type = "img_noise_type" in params ? params["img_noise_type"] : 0;
    switch (type) {
        case "0":
            return img;
        case "1":
            var mean = "img_noise_mean" in params ? parseFloat(params["img_noise_mean"]) : 0;
            var sigma = "img_noise_sigma" in params ? parseFloat(params["img_noise_sigma"]) : 0.001;
            return addGaussianNoiseBoxMuller(img, mean, sigma);
    }
    return img;
}

function addKSpaceNoise(kSpace, img, params) {
    var type = "img_noise_type" in params ? params["img_noise_type"] : 0;
    switch (type) {
        case "0":
            return [kSpace, img];
        case "2":
            var mean = "img_noise_mean" in params ? parseFloat(params["img_noise_mean"]) : 0;
            var sigma = "img_noise_sigma" in params ? parseFloat(params["img_noise_sigma"]) : 0.001;
            kSpace = addGaussianNoiseBoxMuller(kSpace, mean, sigma);
            return [kSpace, undefined]
    }
    return [kSpace, img];
}

function simulateImage(params) {
    try {
        call_count = 0;
        var sequence = "sequence" in params ? params["sequence"] : undefined;
        var S = sequence in imageFunctions ? imageFunctions[sequence] : undefined;
        var xdim = Math.round(params["xdim"]);
        xdim = xdim > 0 ? xdim : k_xdim;
        xdim = xdim > k_xdim ? k_xdim : xdim;
        var ydim = Math.round(params["ydim"]);
        ydim = ydim > 0 ? ydim : k_ydim;
        ydim = ydim > k_ydim ? k_ydim : ydim;
        var zdim = Math.round(params["zdim"]);
        zdim = zdim > 0 ? zdim : k_zdim;
        zdim = zdim > k_zdim ? k_zdim : zdim;
        var nearest = "nearest" in params ? parseInt(params["nearest"]) : 2;
        var cs = "cs" in params ? parseInt(params["cs"]) : 0;
        var cs_filt = "cs_filter" in params ? parseInt(params["cs_filter"]) : 0;
        var cs_mask_t = cs_mask2_t;
        switch(cs_filt) {
            case 0:
                cs_mask_t = cs_mask2_t;
                break;
            case 1:
                cs_mask_t = cs_mask4_t;
                break;
        }
        var fft3d = 'fft' in params ? params['fft'] == '3d' : true;

        if (S == undefined) {
            return undefined;
        }
        var result = new Float32Array(xdim*ydim*zdim);
        //var all_poss = new Set();
        //var double_poss = [];
        for(var z = 0; z<zdim; z++) {
            for(var y = 0;y<ydim; y++) {
                for(var x = 0;x<xdim;x++) {
                    var ipos = z*xdim*ydim + y*xdim + x;
                    var pos = [], d = [];
                    var nz = z*k_zdim/zdim + (k_zdim/zdim/2.0);
                    var ny = y*k_ydim/ydim + (k_ydim/ydim/2.0);
                    var nx = x*k_xdim/xdim + (k_xdim/xdim/2.0);
                    switch (nearest) {
                        case 0:
                            [pos, d] = getNearest(nz, ny, nx);
                            break;
                        case 1:
                            [pos, d] = getNeighbours(nz, ny, nx);
                            break;
                        default:
                            [pos, d] = getVolumeVoxel(nz, ny, nx, zdim, ydim, xdim);
                            break;
                    }

                    /*pos.map(function(p) {
                        if(all_poss.has(p)) {
                            double_poss.push(p);
                        } else {
                            all_poss.add(p);
                        }
                    });*/
                    var ds = d.reduce(function(s, v) { return s+v;});
                    if(ds == 0) { ds = 1/8;}
                    var vals = pos.map(function(v) { return S(v, params) });
                    var val = vals.reduce(function(s, v, i) { return s + v*d[i]; });
                    result[ipos] = val / ds;
                }
            }
        }
        result = addImageNoise(result, params);
        result = new MRImage(xdim, ydim, zdim, [256,256,0,256], result, undefined, params);
        var kspace;
        if (fft3d) {
            kspace = calcKSpace3d(result);
        } else {
            kspace = calcKSpace(result);
        }
        result.kSpace = kspace;

        [k_data_im_re, result] = addKSpaceNoise(k_data_im_re, result, params);
        if (result == undefined) {
            if(fft3d) {
                [result, k_result, p] = inverseKSpace3d(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
            } else {
                [result, k_result, p] = inverseKSpace(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
            }
            result = new MRImage(xdim, ydim, zdim, [256, 256, 0, 256], result, k_result, params);
        }

        var cs_result = undefined;
        var filt_result = undefined;
        var k_result = undefined;
        switch(cs) {
            case 0:
                break;
            case 1:
                k_data_im_re = filter_kspace_cs(k_data_im_re, cs_mask_t, zdim);
                if(fft3d) {
                    [filt_result, k_result, p] = inverseKSpace3d(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
                } else {
                    [filt_result, k_result, p] = inverseKSpace(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
                }
                filt_result = new MRImage(xdim, ydim, zdim, [256,256,0,256], filt_result, k_result, params);
                break;
            case 2:
                k_data_im_re = filter_kspace_cs(k_data_im_re, cs_mask_t, zdim);
                if(fft3d) {
                    [filt_result, k_result, p] = inverseKSpace3d(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
                } else {
                    [filt_result, k_result, p] = inverseKSpace(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
                }
                filt_result = new MRImage(xdim, ydim, zdim, [256,256,0,256], filt_result, k_result, params);
                //k_result = transformKSpace3d(k_data_im_re, fft3d);
                cs_result = new Float32Array(xdim*ydim*zdim);
                if(fft3d) {
                    cs_result = compressed_sensing(k_data_im_re, params);
                } else {
                    for(var z=0;z<zdim;z++) {
                        //console.log("z", z);
                        reply('progress', 100*z/zdim);
                        var k_data_slice = new Float32Array(k_xdim*k_ydim*2);
                        for(var i=0;i<k_data_slice.length;i++) {
                            k_data_slice[i] = k_data_im_re[i+k_xdim*k_ydim*2*z];
                        }
                        var nparams = Object.assign({}, params);
                        nparams["zdim"] = 1;
                        var slice_result = compressed_sensing(k_data_slice, nparams);
                        cs_result.set(slice_result, z*xdim*ydim);
                    }
                    reply('progress', 100);
                }
                cs_result = new MRImage(xdim, ydim, zdim, [256,256,0,256], cs_result, k_result, params);
                break;
            case 3:
                k_data_im_re = filter_kspace_cs(k_data_im_re, cs_mask_t, zdim);
                if(fft3d) {
                    [filt_result, k_result, p] = inverseKSpace3d(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
                } else {
                    [filt_result, k_result, p] = inverseKSpace(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
                }
                filt_result = new MRImage(xdim, ydim, zdim, [256,256,0,256], filt_result, k_result, params);
                //k_result = transformKSpace3d(k_data_im_re, fft3d);            
                var k_cs_re_img;
                if(fft3d) {
                    [cs_result, k_cs_re_img] = compressed_sensing_fast(k_data_im_re, params);
                } else {
                    cs_result = new Float32Array(xdim*ydim*zdim);
                    k_cs_re_img = new Float32Array(xdim*ydim*zdim*2);
                    for(var z=0;z<zdim;z++) {
                        //console.log("z", z);
                        reply('progress', 100*z/zdim);
                        var k_data_slice = new Float32Array(xdim*ydim*2);                    
                        for(var i=0;i<k_data_slice.length;i++) {
                            k_data_slice[i] = k_data_im_re[i+xdim*ydim*2*z];
                        }
                        var nparams = Object.assign({}, params);
                        nparams["zdim"] = 1;
                        var [slice_result, k_slice_re_img] = compressed_sensing_fast(k_data_slice, nparams);
                        cs_result.set(slice_result, z*xdim*ydim);
                        k_cs_re_img.set(k_slice_re_img, z*xdim*ydim*2);
                    }
                    reply('progress', 100);
                }
                var k_cs_result = transformKSpace3d(k_cs_re_img, fft3d, xdim, ydim, zdim);
                //console.log(k_cs_result.reduce((p,v)=>v!=v?p+1:p, 0));
                cs_result = new MRImage(xdim, ydim, zdim, [256,256,0,256], cs_result, k_cs_result, params);
                break;
            case 4:
                break;
            case 5:
                k_result = magnitude_image(k_data_im_re);
                break;
        }
        //console.log("sim voxel function calls:", call_count, "overlap:", call_count/(256*256*256));
        /*var missing_pos = [];
        for(var x=0;x<256;x++){
            for(var y=0;y<256;y++){
                for(var z=0;z<256;z++) {
                    if(!all_poss.has(z*256*256+y*256+x)) {
                        missing_pos.push([x,y,z]);
                    }
                }
            }
        }
        console.log("missing:", missing_pos);
        console.log("double:", double_poss);*/

        if(cs_result != undefined) {
            return [result, filt_result, cs_result];
        } else {
            return result;
        }
    } catch (e) {
        return e;
    }
}

function simulateImageFast(params) {
    try{
        if (ds == undefined) {
            ds = make_dataset(array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_ex_frac, array_na_t2s, array_na_t2f);
        }
        var fft3d = 'fft' in params ? params['fft'] == '3d' : true;
        var xdim = Math.round(params["xdim"]);
        xdim = xdim > 0 ? xdim : k_xdim;
        xdim = xdim > k_xdim ? k_xdim : xdim;
        var ydim = Math.round(params["ydim"]);
        ydim = ydim > 0 ? ydim : k_ydim;
        ydim = ydim > k_ydim ? k_ydim : ydim;
        var zdim = Math.round(params["zdim"]);
        zdim = zdim > 0 ? zdim : k_zdim;
        zdim = zdim > k_zdim ? k_zdim : zdim;
        var nbreg = params["cs_nbreg"];

        if (fft3d) {
            params["cs_callback"] = function (outer) {reply('progress', 100*outer/nbreg)};
        } else {
            params["cs_callback"] = function (z) {reply('progress', 100*z/zdim)};
        }
        [image, kspace, filt_image, filt_kspace, cs_image, cs_kspace] = simulate_fast(ds, params);
        delete params["cs_callback"];

        var k_result = transformKSpace3d(kspace, fft3d, xdim, ydim, zdim);
        result = new MRImage(xdim, ydim, zdim, [256, 256, 0, 256], image, k_result, params, kspace);
        if(cs_image != undefined) {
            var cs_k_result = transformKSpace3d(cs_kspace, fft3d, xdim, ydim, zdim);
            var filt_k_result = transformKSpace3d(filt_kspace, fft3d, xdim, ydim, zdim);
            cs_result = new MRImage(xdim, ydim, zdim, [256, 256, 0, 256], cs_image, cs_k_result, params, cs_kspace);
            filt_result = new MRImage(xdim, ydim, zdim, [256, 256, 0, 256], filt_image, filt_k_result, params, filt_kspace);
            return [result, filt_result, cs_result];
        } else {
            return result;
        }
    } catch (e) {
        console.log(e);
        return e;
    }
}

function profile_list(params_list) {
    var times = [];
    for(var j=0;j<params_list[0];j++) {
        for(var i=1;i<params_list.length;i++) {
            while(times.length<=i) {
                times.push([]);
            }
            reply('profile', profile(params_list[i], times[i]));
        }
    }
    return 0;
}

function profile(params, times) {
    //var times = [];
    var res = 0;
    for(var i=0;i<params["count"];i++) {
        var timer = performance.now();
        if(params["compute"] == "WebASM") {
            res = simulateImageFast(params);
            //console.debug(typeof(res));
        } else {
            res = simulateImage(params);
            //console.debug(typeof(res));
        }
        times.push((performance.now()-timer)/1000);
        if(params["cs"]=="2") {
            console.debug(i, times[times.length-1], res[2].data.reduce((p, c) => p+c));
        } else {
            console.debug(i, times[times.length-1], res.data.reduce((p, c) => p+c));
        }
    }
    return [params, times];
}

var imageFunctions = {
    SE: calcSpinEcho,
    IR: calcInversionRecovery,
    FLASH: calcFlash,
    bSSFP: calcBalancedSSFP,
    PSIF: calcPSIF,
    FISP: calcFISP,
    SGRE: calcSGRE,
    Na: calcNa,
    SQ: calcSQ,
    TQ: calcTQ,
    TQF: calcTQF,
    TQSQR: calcTQ_SQ_Ratio,
};

var queryableFunctions = {
    simulateImage: function(params) {
        reply('result', simulateImage(params));
    },
    simulateImageFast: function(params) {
        reply('result', simulateImageFast(params));
    },
    profile: function(params) {
        reply('endprofile', profile_list(params));
    },
    reco: function(params) {
        reply('result', inverseKSpace3d(k_data_im_re, ...params));
    },
    loadData: async function (path) {
        reply('loadData', await loadDataSet(path));
    },
    kspace: async function() {
        reply('kspace', k_data_im_re);
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

var running = false;
onmessage = async function (oEvent) {
    if (running) {
        console.log("function already running", oEvent);
        defaultReply(oEvent.data);
    }
    if (oEvent.data instanceof Object && oEvent.data.hasOwnProperty('queryMethod') && oEvent.data.hasOwnProperty('queryMethodArguments')) {
        running = true;
        try {
            queryableFunctions[oEvent.data.queryMethod].apply(self, oEvent.data.queryMethodArguments);
        } catch (e) {
            throw e;
        }
        running = false;
    } else {
        defaultReply(oEvent.data);
    }
};