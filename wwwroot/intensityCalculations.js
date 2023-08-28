self.importScripts("./kissfft.js")
self.importScripts("./pako.min.js")
//self.importScripts("./discrete-wavelets.umd.min.js")
self.importScripts("./Mask_CS_Accel2.txt.js");
self.importScripts("./Mask_CS_Accel4.txt.js");

const DIMS = 16;

var ds = undefined;

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
            k_zdim = shape[0];
            k_ydim = shape[1];
            k_xdim = shape[2];
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
    
    array_pd = undefined;
    array_t1 = undefined;
    array_t2 = undefined;
    array_t2s = undefined;
    array_na_mm = undefined;
    array_na_t1 = undefined;
    array_na_ex_frac = undefined;
    array_na_t2s = undefined;
    array_na_t2f = undefined;
    k_data_im_re = undefined;
    k_result = undefined;

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
        ds = make_dataset(k_xdim, k_ydim, k_zdim, array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_ex_frac, array_na_t2s, array_na_t2f);
    } else {
        free_dataset(ds);
        ds = make_dataset(k_xdim, k_ydim, k_zdim, array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_ex_frac, array_na_t2s, array_na_t2f);
    }
    return [array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_t2f, array_na_t2s, k_zdim, k_ydim, k_xdim];
}

function simulateImageFast(params) {
    try{
        if (ds == undefined) {
            ds = make_dataset(k_xdim, k_ydim, k_zdim, array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_ex_frac, array_na_t2s, array_na_t2f);
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

        result = new MRImage(xdim, ydim, zdim, [256, 256, 0, 256], image, kspace, params, kspace);
        if(cs_image != undefined) {
            cs_result = new MRImage(xdim, ydim, zdim, [256, 256, 0, 256], cs_image, cs_kspace, params, cs_kspace);
            filt_result = new MRImage(xdim, ydim, zdim, [256, 256, 0, 256], filt_image, filt_kspace, params, filt_kspace);
            return [result, filt_result, cs_result];
        } else {
            return result;
        }
    } catch (e) {
        //console.log("error", e);
        return e.message;
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
        res = simulateImageFast(params);
        times.push((performance.now()-timer)/1000);
        if(params["cs"]=="2") {
            console.debug(i, times[times.length-1], res[2].data.reduce((p, c) => p+c));
        } else {
            console.debug(i, times[times.length-1], res.data.reduce((p, c) => p+c));
        }
    }
    return [params, times];
}

function debug(msg) {
    switch(msg) {
        default: return ['debug', "not found"];
    }
}

var queryableFunctions = {
    simulateImageFast: function(params) {
        reply('result', simulateImageFast(params));
    },
    profile: function(params) {
        reply('endprofile', profile_list(params));
    },
    loadData: async function (path) {
        reply('loadData', await loadDataSet(path));
    },
    debug: function(msg) {
        result = debug(msg);
        reply(result[0], result[1]);
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