self.importScripts("./kissfft.js")
self.importScripts("./pako.min.js")
self.importScripts("./math.js")
self.importScripts("./discrete-wavelets.umd.min.js")
self.importScripts("./Mask_CS_Accel2.txt.js");
self.importScripts("./Mask_CS_Accel4.txt.js");

//var xdim = 256;
//var ydim = 256;
//var zdim = 256;

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
    constructor(xdim, ydim, zdim, kspace_filt, data, k_data, params) {
        this.xdim = xdim;
        this.ydim = ydim;
        this.zdim = zdim;
        this.data = data;
        this.kSpace = k_data;
        this.kspace_filt = kspace_filt;
        this.params = params;
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

function _fft2d(data, xdim, ydim) {
    var k_data_im_re = new Float32Array(xdim * ydim * 2);
    k_data_im_re.set(data);
    for(var i=0;i<k_ydim;i++) {
        my_fft_stride(k_data_im_re, 2*i*xdim, 1, xdim);
    }
    for(var i=0;i<k_xdim;i++) {
        my_fft_stride(k_data_im_re, 2*i, xdim, ydim);
    }
    return k_data_im_re;
}

function fft2d(data, xdim, ydim) {
    return kfft2d(data, xdim, ydim);
}

function _ifft2d(data, xdim, ydim) {
    var k_data_im_re = new Float32Array(xdim * ydim * 2);
    k_data_im_re.set(data);
    for(var i=0;i<k_ydim;i++) {
        my_ifft_stride(k_data_im_re, 2*i*xdim, 1, xdim);
    }
    for(var i=0;i<k_xdim;i++) {
        my_ifft_stride(k_data_im_re, 2*i, xdim, ydim);
    }
    return k_data_im_re;
}

function ifft2d(data, xdim, ydim) {
    return kifft2d(data, xdim, ydim);
}

function fft3d(data, xdim, ydim, zdim) {
    var k_data_im_re = new Float32Array(xdim * ydim * zdim * 2);
    if(data.length == k_data_im_re.length) {
        k_data_im_re.set(data);
    } else {
        k_data_im_re.fill(0);
        for(var i=0;i<data.length;i++) {
            k_data_im_re[2*i] = data[i];
        }
    }
    var slice = new Float32Array(xdim*ydim*2);
    for(var z=0;z<zdim;z++) {
        for(var i=0;i<slice.length;i++) {
            slice[i] = k_data_im_re[i+z*xdim*ydim*2];
        }
        k_data_im_re.set(fft2d(slice, xdim, ydim), z*xdim*ydim*2);
    }
    for(var i=0;i<xdim*ydim;i++) {
        my_fft_stride(k_data_im_re, 2*i, xdim*ydim, zdim);
    }
    return k_data_im_re;
}

function _fft3d(data, xdim, ydim, zdim) {
    return kfft3d(data, xdim, ydim, zdim);
}

function ifft3d(data, xdim, ydim, zdim) {
    var img = new Float32Array(xdim * ydim * zdim * 2);
    img.set(data);
    var slice = new Float32Array(xdim*ydim*2);
    for(var z=0;z<zdim;z++) {
        for(var i=0;i<slice.length;i++) {
            slice[i] = img[i+z*xdim*ydim*2];
        }
        img.set(ifft2d(slice, xdim, ydim), z*xdim*ydim*2);
    }
    for(var i=0;i<xdim*ydim;i++) {
        my_ifft_stride(img, 2*i, k_xdim*k_ydim, zdim);
    }
    return img;
}

function _ifft3d(data, xdim, ydim, zdim) {
    return kfft3d(data, xdim, ydim, zdim);
}

function calcKSpace(result) {
    var xdim = result.xdim;
    var ydim = result.ydim;
    var zdim = result.zdim;
    var xoff = Math.floor((k_xdim-xdim) / 2);
    var yoff = Math.floor((k_ydim-ydim) / 2);
    xoff = 0;
    yoff = 0;
    k_data_im_re = new Float32Array(k_xdim * k_ydim * zdim * 2);
    k_data_im_re.fill(0);
    //console.log("calc", k_xdim, k_ydim, zdim, k_data_im_re.length, k_data_im_re.reduce((p,v)=>v!=v?p+1:p, 0))
    for (var z = 0; z < zdim; z++) {
        var slice_data = new Float32Array(k_xdim * k_ydim * 2);
        slice_data.fill(0);
        for (var x = 0; x < xdim; x++) {
            for (var y=0; y < ydim; y++) {
                slice_data[2*(xoff+x + (yoff+y)*k_xdim)] = result.data[x + y*xdim + z*xdim*ydim];
            }
        }
        var fft_res = fft2d(slice_data, k_xdim, k_ydim);
        //console.log(z, z*k_xdim*k_ydim*2, fft_res.length, fft_res.reduce((p,v)=>v!=v?p+1:p,0));
        k_data_im_re.set(fft_res, z * k_xdim * k_ydim * 2)
    }
    //console.log("calc", k_data_im_re.length, k_data_im_re.reduce((p,v)=>v!=v?p+1:p, 0))
    result.kSpace = transformKSpace3d(k_data_im_re, false, zdim);
    return result.kSpace;
}

function calcKSpace3d(result) {
    var xdim = result.xdim;
    var ydim = result.ydim;
    var zdim = result.zdim;
    var c_data = new Float32Array(k_xdim * k_ydim * zdim * 2);
    c_data.fill(0);
    for (var x = 0; x < xdim; x++) {
        for (var y=0; y < ydim; y++) {
            for(var z=0; z< zdim; z++) {
                c_data[2*(x + y*k_xdim + z*k_xdim*k_ydim)] = result.data[x + y*xdim + z*xdim*ydim];
            }
        }
    }
    k_data_im_re = fft3d(c_data, k_xdim, k_ydim, zdim);
    result.kSpace = transformKSpace3d(k_data_im_re, true, zdim);
    return result.kSpace;
}

function transformKSpace3d(k_data_im_re, fft3d, zdim) {
    var k_result = new Float32Array(k_xdim * k_ydim * zdim);
    for (var z = 0; z < zdim; z++) {
        var slice_data = new Float32Array(k_xdim * k_ydim*2);
        for (var x = 0; x < k_xdim; x++) {
            for (var y=0; y < k_ydim; y++) {
                slice_data[2*(x + y*k_xdim)] = k_data_im_re[2*(x + y*k_xdim + z*k_xdim*k_ydim)];
                slice_data[2*(x + y*k_xdim)+1] = k_data_im_re[2*(x + y*k_xdim + z*k_xdim*k_ydim) + 1];
                if(slice_data[2*(x + y*k_xdim)] == NaN || slice_data[2*(x + y*k_xdim)] == undefined || slice_data[2*(x + y*k_xdim)+1] == undefined || slice_data[2*(x + y*k_xdim)+1] == NaN) {
                    console.log("fft output undefined", x,y,z);
                }
            }
        }
        if(fft3d) {
            k_result.set(transformKSpace(slice_data), ((z+zdim/2)%zdim) * k_xdim * k_ydim);
        } else {
            k_result.set(transformKSpace(slice_data), z * k_xdim * k_ydim);
        }
    }
    var max = k_result.reduce((p,v) => Math.max(p,v));
    var min = k_result.reduce((p,v) => Math.min(p,v));
    for(var i=0;i<k_result.length;i++) {
        k_result[i] = 255*(k_result[i]-min) / (max-min);
    }
    return k_result
}

function transformKSpace(fft_res) {
    var k_data = new Float32Array(k_xdim * k_ydim);
    var k_result = new Float32Array(k_xdim * k_ydim);
    var maxval = 0;
    var minval = 999999999;
    for (var i = 0; i < k_data.length; i++) {
        k_data[i] = Math.sqrt(fft_res[2 * i] * fft_res[2 * i] + fft_res[2 * i + 1] * fft_res[2 * i + 1]);
        maxval = Math.max(maxval, k_data[i]);
        minval = Math.min(minval, k_data[i]);
    }
    for(var x = 0;x<k_xdim;x++) {
        for(var y=0;y<k_ydim;y++) {
            var i = x + y*k_xdim;
            var val = (k_data[i] - minval) * 255 / (maxval-minval);

            var j = ((x+k_xdim/2)%k_xdim)+((y+k_ydim/2)%k_ydim)*k_xdim;
            k_result[j] = k_data[i];
            //k_result[j] = val;
            //j = (k_xdim-x-k_xdim/2-1)+((y+k_ydim/2)%k_ydim)*k_xdim;
            //k_result[j] = val;
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
    var k_result = new Float32Array(k_xdim * k_ydim * zdim);
    k_result.fill(0);
    var zend = Math.max(Math.floor(zdim/25), 1);
    for (var z = 0; z < zdim; z++) {
        if (2*z*k_xdim*k_ydim > kSpace.length) { zdim = z; break; }
        //var slice_data = kSpace.subarray(z * k_xdim * k_ydim * 2, (z + 1) * k_xdim * k_ydim * 2);
        var slice_data = new Float32Array(k_xdim*k_ydim*2);
        for(var i=0;i<slice_data.length;i++) {
            slice_data[i] = kSpace[i+z*k_xdim*k_ydim*2];
        }

        var _fmax = fmax;
        if (z<zend) {
            _fmax = z/zend*z/zend * fmax;
        }
        if (z > zdim-zend) {
            _fmax = (zdim-z)/zend*(zdim-z)/zend * fmax;
        }
        var mapKSpace = genMapKSpace(k_xdim, k_ydim, xlines, ylines, fmin, _fmax)
        var input_data = slice_data.map(mapKSpace);
        k_result.set(transformKSpace(input_data), z * k_xdim * k_ydim);

        if (!noIfft) {
            var fft_result = ifft2d(input_data, k_xdim, k_ydim);

            var xoff = Math.floor((k_xdim-xdim) / 2);
            var yoff = Math.floor((k_ydim-ydim) / 2);
            for (var x = 0; x < xdim; x++) {
                for (var y=0; y < ydim; y++) {
                    result[x+y*xdim+z*xdim*ydim] = Math.sqrt(fft_result[2*(xoff+x + (yoff+y)*k_xdim)]*fft_result[2*(xoff+x + (yoff+y)*k_xdim)]+fft_result[2*(xoff+x + (yoff+y)*k_xdim)+1]*fft_result[2*(xoff+x + (yoff+y)*k_xdim)+1]);
                    //result[x+y*xdim+z*xdim*ydim] = mag_result[xoff+x+(yoff+y)*k_xdim]
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
    //var mapKSpace = genMapKSpace(xlines, ylines, fmin, fmax)
    var result = new Float32Array(xdim * ydim * zdim);
    //var k_result = new Float32Array(k_xdim * k_ydim * zdim);
    //k_result.fill(0);
    var zend = Math.max(Math.floor(zdim/25), 1);
    
    var _fmax = fmax;
    if (z<zend) {
        _fmax = z/zend*z/zend * fmax;
    }
    if (z > zdim-zend) {
        _fmax = (zdim-z)/zend*(zdim-z)/zend * fmax;
    }
    var mapKSpace = genMapKSpace3d(k_xdim, k_ydim, xlines, ylines, fmin, fmax)
    var input_data = kSpace.map(mapKSpace);
    var k_result = transformKSpace3d(kSpace, true, zdim);

    if (!noIfft) {
        var fft_result = ifft3d(kSpace, k_xdim, k_ydim, k_zdim);
        /*var maxval = 0;
        var minval = 999999999;
        for (var i = 0; i < img_result.length; i++) {
            maxval = Math.max(maxval, img_result[i]);
            minval = Math.min(minval, img_result[i]);
        }*/
        var xoff = Math.floor((k_xdim-xdim) / 2);
        var yoff = Math.floor((k_ydim-ydim) / 2);
        for (var x = 0; x < xdim; x++) {
            for (var y=0; y < ydim; y++) {
                for(var z=0;z<zdim;z++) {
                    var pos = 2*(xoff+x + (yoff+y)*k_xdim+z*k_xdim*k_ydim);
                    result[x+y*xdim+z*xdim*ydim] = Math.sqrt(fft_result[pos]*fft_result[pos]+fft_result[pos+1]*fft_result[pos+1]);
                //result[x+y*xdim+z*xdim*ydim] = mag_result[xoff+x+(yoff+y)*k_xdim]
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

function calcSpinEcho(pos, params) {
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
    var te = params["te"];
    var tr = params["tr"];
    //Na: s(t) = VbCb (1-exp(-tr/t1))[0.6*exp(-t/T2s) + 0.4*exp(-t/T2L)] + VfrCfrafr(1-exp(-tr/t1))exp(-t/T2fr) + n(t). 
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
    //var slice_data = k_space.subarray(z * k_xdim * k_ydim * 2, (z + 1) * k_xdim * k_ydim * 2);
    //var initial_img = irfft2d(slice_data, k_xdim, k_ydim)
    var initial_img = idwt2d(k_space, "db1", k_xdim, k_ydim);
    var update_img = sparsify_transform(initial_img);
    for(var i=0;i<update_img.length;i++) {
        if(update_img[i] < threshold) { 
            update_img[i] = 0;
        }
    }
    update_img = i_sparsify_transform(update_img);
    //var k_space_d = rfft2d(update_img);
    var k_space_d = dwt2d(update_img, "db1", k_xdim, k_ydim);
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
        for(var y=0;y<k_ydim;y++) {
            for(var x=0;x<k_xdim;x++) {
                var pos = 2*(x+y*k_xdim+z*k_xdim*k_ydim);
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
    f_uker = data_ndims==2? fft2d(f_uker, k_xdim, k_ydim) : fft3d(f_uker, k_xdim, k_ydim, k_zdim);
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
      i_murf = data_ndims==2? ifft2d(i_murf, k_xdim, k_ydim) : ifft3d(i_murf, k_xdim, k_ydim, k_zdim);
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
        var f_rhs = data_ndims==2? fft2d(i_rhs, k_xdim, k_ydim) : fft3d(i_rhs, k_xdim, k_ydim, k_zdim);
        //print_stats("rhs out", f_rhs);
        for(var i=0;i<f_rhs.length;i+=2) {
            var x = f_rhs[i]; var y = f_rhs[i+1]; var u = f_uker[i]; var v = f_uker[i+1];
            f_rhs[i] = x*u-y*v;
            f_rhs[i+1] = x*v+y*u;
        }
        //print_stats("img_in", f_rhs);
        img = data_ndims==2? ifft2d(f_rhs, k_xdim, k_ydim) : ifft3d(f_rhs, k_xdim, k_ydim, k_zdim);
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
      var k_img = data_ndims==2? fft2d(img, k_xdim, k_ydim) : fft3d(img, k_xdim, k_ydim, k_zdim);
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
    return [z1*k_xdim*k_ydim+y1*k_xdim+x1, d];
}

function getNearest(z,y,x) {
    var [pos, d] = getNeighbours(z,y,x);
    var [max,mi] = d.reduce(function(o, v, i) { 
        if(o[0] < v) {
            return [v, i];
        } else { return o; }
    }, [d[0], 0]);
    return [[pos[mi]], [max]];
}

function getNeighbours(z,y,x) {
    var [dpos1, d1] = dist(z,y,x, Math.floor(z*k_zdim/zdim), Math.floor(y*k_ydim/ydim), Math.floor(x*k_xdim/xdim));
    var [dpos2, d2] = dist(z,y,x, Math.floor(z*k_zdim/zdim), Math.floor(y*k_ydim/ydim), Math.ceil(x*k_xdim/xdim));
    var [dpos3, d3] = dist(z,y,x, Math.floor(z*k_zdim/zdim), Math.ceil(y*k_ydim/ydim),  Math.floor(x*k_xdim/xdim));
    var [dpos4, d4] = dist(z,y,x, Math.floor(z*k_zdim/zdim), Math.ceil(y*k_ydim/ydim),  Math.ceil(x*k_xdim/xdim));
    var [dpos5, d5] = dist(z,y,x, Math.ceil(z*k_zdim/zdim),  Math.floor(y*k_ydim/ydim), Math.floor(x*k_xdim/xdim));
    var [dpos6, d6] = dist(z,y,x, Math.ceil(z*k_zdim/zdim),  Math.floor(y*k_ydim/ydim), Math.ceil(x*k_xdim/xdim));
    var [dpos7, d7] = dist(z,y,x, Math.ceil(z*k_zdim/zdim),  Math.ceil(y*k_ydim/ydim),  Math.floor(x*k_xdim/xdim));
    var [dpos8, d8] = dist(z,y,x, Math.ceil(z*k_zdim/zdim),  Math.ceil(y*k_ydim/ydim),  Math.ceil(x*k_xdim/xdim));
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
    for(var xi=Math.floor(x-xs);xi<Math.ceil(x+xs);xi+=1) {
        for(var yi=Math.floor(y-ys);yi<Math.ceil(y+ys);yi+=1) {
            for(var zi=Math.floor(z-zs);zi<Math.ceil(z+zs);zi+=1) {
                var pos = zi*k_xdim*k_ydim+yi*k_xdim+xi;
                if(pos >= 0 && pos < k_xdim*k_ydim*k_zdim && xi>=0 && xi<k_xdim && yi>=0 && yi<k_ydim && zi>=0 && zi<k_zdim) {
                    poss.push(pos);
                    ds.push(1);
                }
            }
        }
    }
    //if(poss.length == 0) {
        //console.log("x", k_xdim, xdim, x, xs, Math.round(x-xs), Math.ceil(x+xs));
        //console.log("y", k_ydim, ydim, y, ys, Math.round(y-ys), Math.ceil(y+ys));
        //console.log("z", k_zdim, zdim, z, zs, Math.round(z-zs), Math.ceil(z+zs));
    //}
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
    var sequence = "sequence" in params ? params["sequence"] : undefined;
    var S = sequence in imageFunctions ? imageFunctions[sequence] : undefined;
    //console.log(params, sequence, S);
    var xdim = Math.round(params["xdim"]);
    xdim = xdim > 0 ? xdim : k_xdim;
    xdim = xdim > k_xdim ? k_xdim : xdim;
    var ydim = Math.round(params["ydim"]);
    ydim = ydim > 0 ? ydim : k_ydim;
    ydim = ydim > k_ydim ? k_ydim : ydim;
    var zdim = Math.round(params["zdim"]);
    zdim = zdim > 0 ? zdim : k_zdim;
    zdim = zdim > k_zdim ? k_zdim : zdim;
    //console.log(xdim, ydim, zdim);
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
    for(var z = 0; z<zdim; z++) {
        for(var y = 0;y<ydim; y++) {
            for(var x = 0;x<xdim;x++) {
                var ipos = z*xdim*ydim + y*xdim + x;
                var pos = [], d = [];
                switch (nearest) {
                    case 0:
                        [pos, d] = getNearest(z*k_zdim/zdim, y*k_ydim/ydim, x*k_xdim/xdim);
                        break;
                    case 1:
                        [pos, d] = getNeighbours(z*k_zdim/zdim, y*k_ydim/ydim, x*k_xdim/xdim);
                        break;
                    default:
                        [pos, d] = getVolumeVoxel(z*k_zdim/zdim, y*k_ydim/ydim, x*k_xdim/xdim, zdim, ydim, xdim);
                        break;
                }
                var ds = d.reduce(function(s, v) { return s+v;});
                if(ds == 0) { ds = 1/8;}
                var vals = pos.map(function(v) { return S(v, params) });
                var val = vals.reduce(function(s, v, i) { return s + v*d[i]/ds; });
                result[ipos] = val;
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
    //console.log(zdim, k_zdim, xdim*ydim*zdim, kspace.length, result.data.length, kspace.reduce((p,v)=>v!=v?p+1:p, 0));
    //console.log(result.kSpace.reduce((p,v)=>v!=v?p+1:p, 0));
    result.kSpace = kspace;
    //console.log(result.kSpace.reduce((p,v)=>v!=v?p+1:p, 0));

    [k_data_im_re, result] = addKSpaceNoise(k_data_im_re, result, params);
    /*for(var z=0;z<zdim;z++) {
        k_data_im_re.set(phant_data, z*k_xdim*k_ydim*2);
        k_result.set(transformKSpace(phant_data),z*k_xdim*k_ydim);
        result = undefined;
    }*/
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
                    var k_data_slice = new Float32Array(k_xdim*k_ydim*2);                    
                    for(var i=0;i<k_data_slice.length;i++) {
                        k_data_slice[i] = k_data_im_re[i+k_xdim*k_ydim*2*z];
                    }
                    var nparams = Object.assign({}, params);
                    nparams["zdim"] = 1;
                    var [slice_result, k_slice_re_img] = compressed_sensing_fast(k_data_slice, nparams);
                    cs_result.set(slice_result, z*xdim*ydim);
                    k_cs_re_img.set(k_slice_re_img, z*xdim*ydim*2);
                    /*for(var i=0;i<slice_result.length;i++) {
                        cs_result[z*xdim*ydim+i] = slice_result[i];
                    }*/
                }
                reply('progress', 100);
            }
            var k_cs_result = transformKSpace3d(k_cs_re_img, fft3d, zdim);
            //console.log(k_cs_result.reduce((p,v)=>v!=v?p+1:p, 0));
            cs_result = new MRImage(xdim, ydim, zdim, [256,256,0,256], cs_result, k_cs_result, params);
            break;
        case 4:
            break;
        case 5:
            k_result = magnitude_image(k_data_im_re);
            break;
    }
    if (result == undefined) {
        if(fft3d) {
            [result, k_result, p] = inverseKSpace3d(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
        } else {
            [result, k_result, p] = inverseKSpace(k_data_im_re, xdim, ydim, zdim, 256, 256, 0, 256, false)
        }
        result = new MRImage(xdim, ydim, zdim, [256, 256, 0, 256], result, k_result, params);
    }
    if(cs_result != undefined) {
        return [result, filt_result, cs_result];
    } else {
        return result;
    }
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
    reco: function (params) {
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