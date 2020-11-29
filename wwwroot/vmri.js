
function QueryableWorker(url) {
    var instance = this,
        worker = new Worker(url, {type:"module"}),
        listeners = {};

    this.defaultListener = function () { };

    this.postMessage = function (message) {
        worker.postMessage(message);
    }

    this.terminate = function () {
        worker.terminate();
    }

    this.addListener = function (name, listener) {
        listeners[name] = listener;
    }

    this.removeListener = function (name) {
        delete listeners[name];
    }

    /* 
      This functions takes at least one argument, the method name we want to query.
      Then we can pass in the arguments that the method needs.
    */
    this.sendQuery = function () {
        if (arguments.length < 1) {
            throw new TypeError('QueryableWorker.sendQuery takes at least one argument');
            return;
        }
        worker.postMessage({
            'queryMethod': arguments[0],
            'queryMethodArguments': Array.prototype.slice.call(arguments, 1)
        });
    }

    worker.onmessage = function (event) {
        if (event.data instanceof Object &&
            event.data.hasOwnProperty('queryMethodListener') &&
            event.data.hasOwnProperty('queryMethodArguments')) {
            listeners[event.data.queryMethodListener].apply(instance, event.data.queryMethodArguments);
        } else {
            this.defaultListener.call(instance, event.data);
        }
    }
}
const w = new QueryableWorker("intensityCalculations.js");

const xdim = 256;
const ydim = 256;
const zdim = 256;
var array_pd, array_t1, array_t2, result, k_data_im_re, slice_data;

const loadDataMessageHandler = function (data) {
    array_pd = data[0];
    array_t1 = data[1];
    array_t2 = data[2];
    r = document.getElementById("content");
    r.classList.remove("hidden");
    spin = document.getElementById("datasetLoading");
    spin.classList.add("hidden");
    displayDataSet();
};
w.addListener('loadData', loadDataMessageHandler);

var imgResult;
var kResult;
const resultMessageHandler = function (data) {
    imgResult = data[0];
    kResult = data[1];
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    r.classList.remove("hidden");
    spin.classList.add("hidden");
    displayAndWindow3DImage();
};
w.addListener('result', resultMessageHandler);



// Name, Label, T1, T2, T2*, PD
const params = {
    0: ["Background", 0, 0, 0, 0, 0, "subject04_bck_v.bin"],
    1: ["CSF", 1, 2569.0, 329, 58, 1, "subject04_csf_v.bin"],
    2: ["Grey Matter", 2.0, 833, 83, 69, 0.86, "subject04_gm_v.bin"],
    3: ["White Matter", 3.0, 500, 70, 61, 0.77, "subject04_wm_v.bin"],
    4: ["Fat", 4, 350.0, 70.0, 58, 1, "subject04_fat_v.bin"],
    5: ["Muscle", 5, 900.0, 47, 30, 1, "subject04_muscles_v.bin"],
    6: ["Muscle / Skin", 6, 569.0, 329, 58, 1, "subject04_muscles_skin_v.bin"],
    7: ["Skull", 7, 0, 0, 0, 0, "subject04_skull_v.bin"],
    8: ["Vessels", 8, 0, 0, 0, 0, "subject04_vessels_v.bin"],
    9: ["Around fat", 9, 500.0, 70, 61, 0.77, "subject04_fat2_v.bin"],
    10: ["Dura Matter", 10, 2569.0, 329, 58, 1, "subject04_dura_v.bin"],
    11: ["Bone Marrow", 11, 500.0, 70, 61, 0.77, "subject04_marrow_v.bin"]
}

function display3DImage(canvas, data) {
    //canvas = document.createElement('canvas');
    ctx = canvas.getContext('2d');
    canvas.width = xdim;
    canvas.height = ydim;
    idata = ctx.createImageData(xdim, ydim);

    slice = parseInt(document.getElementById("slice").value);

    image_result = new Uint8ClampedArray(xdim * ydim * 4);
    min = data[0];
    max = data[0];
    for (var x = 0; x < data.length; x++) {
        if (data[x] < min) min = data[x];
        if (data[x] > max) max = data[x];
    }
    for (var x = 0; x < xdim * ydim; x++) {
        image_result[4 * x] = 255.0 * ((data[x + slice * xdim * ydim] - min) / (max - min))
        image_result[4 * x + 1] = 255.0 * ((data[x + slice * xdim * ydim] - min) / (max - min))
        image_result[4 * x + 2] = 255.0 * ((data[x + slice * xdim * ydim] - min) / (max - min))
        image_result[4 * x + 3] = 255
    }

    idata.data.set(image_result);
    ctx.putImageData(idata, 0, 0);
    //dataUri = canvas.toDataURL();
    //element.src = dataUri;
}

function displayAndWindow3DImage() {
    //canvas = document.createElement('canvas');
    res_element = document.getElementById("result");
    canvas = res_element.getElementsByTagName("canvas")[0];
    ctx = canvas.getContext('2d');
    canvas.width = xdim;
    canvas.height = ydim;
    idata = ctx.createImageData(xdim, ydim);

    inputs = res_element.getElementsByTagName("input");
    in_wc = inputs[0]
    in_ww = inputs[1]
    in_slice = inputs[2]
    slice = parseInt(in_slice.value);

    image_result = new Uint8ClampedArray(xdim * ydim * 4);
    var ww = parseFloat(in_ww.value) * 0.5
    var wc = parseFloat(in_wc.value)

    for (var x = 0; x < xdim * ydim; x++) {
        val = imgResult[x + slice * xdim * ydim] * 4096
        if (val <= (wc - ww)) { val = 0 }
        else if (val >= (wc + ww)) { val = 255; }
        else { val = 255 * (val - (wc - ww)) / (ww) }
        image_result[4 * x] = val
        image_result[4 * x + 1] = val
        image_result[4 * x + 2] = val
        image_result[4 * x + 3] = 255
    }

    idata.data.set(image_result);
    ctx.putImageData(idata, 0, 0);

    k_canvas = res_element.getElementsByTagName("canvas")[1];
    k_canvas.width = xdim;
    k_canvas.height = ydim
    k_ctx = k_canvas.getContext('2d');
    kdata = k_ctx.createImageData(xdim, ydim);

    k_result = new Uint8ClampedArray(xdim * ydim * 4);
    for (var x = 0; x < xdim * ydim; x++) {
        val = kResult[x + slice * xdim * ydim]
        k_result[4 * x] = val
        k_result[4 * x + 1] = val
        k_result[4 * x + 2] = val
        k_result[4 * x + 3] = 255
    }

    kdata.data.set(k_result);
    k_ctx.putImageData(kdata, 0, 0);
}

function loadFuzzyDataSet() {
    content = document.getElementById("content");
    content.classList.add("hidden");
    return new Promise(async (resolve, reject) => {

        slice = document.getElementById("slice");
        slice.max = zdim;
        slice.value = 100;

        spin = document.getElementById("datasetLoading");
        spin.classList.remove("hidden");

        w.sendQuery("loadData");
        resolve();
    });
}

function loadDataSet() {
    GeoTIFF.fromUrl("PD.tif").then(tiff => {
        tiff.getImage().then(image => {
            tiff_pd = image;
            image.readRasters().then(data => {
                array_pd = data[0];
                displayImage(document.getElementById("imgPD"), array_pd);
            })
        })
    })

    GeoTIFF.fromUrl("T1.tif").then(tiff => {
        tiff.getImage().then(image => {
            tiff_t1 = image;
            image.readRasters().then(data => {
                array_t1 = data[0];
                displayImage(document.getElementById("imgT1"), array_t1);
            })
        })
    })

    GeoTIFF.fromUrl("T2.tif").then(tiff => {
        tiff.getImage().then(image => {
            tiff_t2 = image;
            image.readRasters().then(data => {
                array_t2 = data[0];
                displayImage(document.getElementById("imgT2"), array_t2);
            })
        })
    })
}

function setSpinEcho() {
    var seTab = document.getElementById("params-se");
    var irTab = document.getElementById("params-ir");
    var seTabHead = document.getElementById("params-se-tab");
    var irTabHead = document.getElementById("params-ir-tab");
    irTab.classList.remove("active", "show");
    seTab.classList.add("active", "show");
    irTabHead.classList.remove("active");
    seTabHead.classList.add("active");
}

function spinEcho() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var te = parseFloat(document.getElementById("se_ti").value)
    var tr = parseFloat(document.getElementById("se_tr").value)

    w.sendQuery("spinEcho", te, tr);
}

function setInversionRecovery() {
    var seTab = document.getElementById("params-se");
    var irTab = document.getElementById("params-ir");
    var seTabHead = document.getElementById("params-se-tab");
    var irTabHead = document.getElementById("params-ir-tab");
    seTab.classList.remove("active", "show");
    irTab.classList.add("active", "show");
    seTabHead.classList.remove("active");
    irTabHead.classList.add("active");
}

async function inversionRecovery() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var ti = parseFloat(document.getElementById("ir_ti").value)
    var tr = parseFloat(document.getElementById("ir_tr").value)

    w.sendQuery("inversionRecovery", ti, tr);
}

function displayDataSet() {
    display3DImage(document.getElementById("imgPD"), array_pd);
    display3DImage(document.getElementById("imgT1"), array_t1);
    display3DImage(document.getElementById("imgT2"), array_t2);
}