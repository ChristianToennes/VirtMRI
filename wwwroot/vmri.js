const jetmap = [[0.0, 0.0, 0.5], [0.0, 0.0, 0.5178], [0.0, 0.0, 0.5357], [0.0, 0.0, 0.5535], [0.0, 0.0, 0.5713], [0.0, 0.0, 0.5891], [0.0, 0.0, 0.607], [0.0, 0.0, 0.6248], [0.0, 0.0, 0.6426], [0.0, 0.0, 0.6604], [0.0, 0.0, 0.6783], [0.0, 0.0, 0.6961], [0.0, 0.0, 0.7139], [0.0, 0.0, 0.7317], [0.0, 0.0, 0.7496], [0.0, 0.0, 0.7674], [0.0, 0.0, 0.7852], [0.0, 0.0, 0.803], [0.0, 0.0, 0.8209], [0.0, 0.0, 0.8387], [0.0, 0.0, 0.8565], [0.0, 0.0, 0.8743], [0.0, 0.0, 0.8922], [0.0, 0.0, 0.91], [0.0, 0.0, 0.9278], [0.0, 0.0, 0.9456], [0.0, 0.0, 0.9635], [0.0, 0.0, 0.9813], [0.0, 0.0, 0.9991], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.002, 1.0], [0.0, 0.0176, 1.0], [0.0, 0.0333, 1.0], [0.0, 0.049, 1.0], [0.0, 0.0647, 1.0], [0.0, 0.0804, 1.0], [0.0, 0.0961, 1.0], [0.0, 0.1118, 1.0], [0.0, 0.1275, 1.0], [0.0, 0.1431, 1.0], [0.0, 0.1588, 1.0], [0.0, 0.1745, 1.0], [0.0, 0.1902, 1.0], [0.0, 0.2059, 1.0], [0.0, 0.2216, 1.0], [0.0, 0.2373, 1.0], [0.0, 0.2529, 1.0], [0.0, 0.2686, 1.0], [0.0, 0.2843, 1.0], [0.0, 0.3, 1.0], [0.0, 0.3157, 1.0], [0.0, 0.3314, 1.0], [0.0, 0.3471, 1.0], [0.0, 0.3627, 1.0], [0.0, 0.3784, 1.0], [0.0, 0.3941, 1.0], [0.0, 0.4098, 1.0], [0.0, 0.4255, 1.0], [0.0, 0.4412, 1.0], [0.0, 0.4569, 1.0], [0.0, 0.4725, 1.0], [0.0, 0.4882, 1.0], [0.0, 0.5039, 1.0], [0.0, 0.5196, 1.0], [0.0, 0.5353, 1.0], [0.0, 0.551, 1.0], [0.0, 0.5667, 1.0], [0.0, 0.5824, 1.0], [0.0, 0.598, 1.0], [0.0, 0.6137, 1.0], [0.0, 0.6294, 1.0], [0.0, 0.6451, 1.0], [0.0, 0.6608, 1.0], [0.0, 0.6765, 1.0], [0.0, 0.6922, 1.0], [0.0, 0.7078, 1.0], [0.0, 0.7235, 1.0], [0.0, 0.7392, 1.0], [0.0, 0.7549, 1.0], [0.0, 0.7706, 1.0], [0.0, 0.7863, 1.0], [0.0, 0.802, 1.0], [0.0, 0.8176, 1.0], [0.0, 0.8333, 1.0], [0.0, 0.849, 1.0], [0.0, 0.8647, 0.9962], [0.0, 0.8804, 0.9836], [0.0, 0.8961, 0.9709], [0.0095, 0.9118, 0.9583], [0.0221, 0.9275, 0.9456], [0.0348, 0.9431, 0.933], [0.0474, 0.9588, 0.9203], [0.0601, 0.9745, 0.9077], [0.0727, 0.9902, 0.895], [0.0854, 1.0, 0.8824], [0.098, 1.0, 0.8697], [0.1107, 1.0, 0.8571], [0.1233, 1.0, 0.8444], [0.136, 1.0, 0.8318], [0.1486, 1.0, 0.8191], [0.1613, 1.0, 0.8065], [0.1739, 1.0, 0.7938], [0.1866, 1.0, 0.7812], [0.1992, 1.0, 0.7685], [0.2119, 1.0, 0.7559], [0.2245, 1.0, 0.7432], [0.2372, 1.0, 0.7306], [0.2498, 1.0, 0.7179], [0.2625, 1.0, 0.7052], [0.2751, 1.0, 0.6926], [0.2878, 1.0, 0.6799], [0.3004, 1.0, 0.6673], [0.3131, 1.0, 0.6546], [0.3257, 1.0, 0.642], [0.3384, 1.0, 0.6293], [0.351, 1.0, 0.6167], [0.3637, 1.0, 0.604], [0.3763, 1.0, 0.5914], [0.389, 1.0, 0.5787], [0.4016, 1.0, 0.5661], [0.4143, 1.0, 0.5534], [0.4269, 1.0, 0.5408], [0.4396, 1.0, 0.5281], [0.4522, 1.0, 0.5155], [0.4649, 1.0, 0.5028], [0.4775, 1.0, 0.4902], [0.4902, 1.0, 0.4775], [0.5028, 1.0, 0.4649], [0.5155, 1.0, 0.4522], [0.5281, 1.0, 0.4396], [0.5408, 1.0, 0.4269], [0.5534, 1.0, 0.4143], [0.5661, 1.0, 0.4016], [0.5787, 1.0, 0.389], [0.5914, 1.0, 0.3763], [0.604, 1.0, 0.3637], [0.6167, 1.0, 0.351], [0.6293, 1.0, 0.3384], [0.642, 1.0, 0.3257], [0.6546, 1.0, 0.3131], [0.6673, 1.0, 0.3004], [0.6799, 1.0, 0.2878], [0.6926, 1.0, 0.2751], [0.7052, 1.0, 0.2625], [0.7179, 1.0, 0.2498], [0.7306, 1.0, 0.2372], [0.7432, 1.0, 0.2245], [0.7559, 1.0, 0.2119], [0.7685, 1.0, 0.1992], [0.7812, 1.0, 0.1866], [0.7938, 1.0, 0.1739], [0.8065, 1.0, 0.1613], [0.8191, 1.0, 0.1486], [0.8318, 1.0, 0.136], [0.8444, 1.0, 0.1233], [0.8571, 1.0, 0.1107], [0.8697, 1.0, 0.098], [0.8824, 1.0, 0.0854], [0.895, 1.0, 0.0727], [0.9077, 1.0, 0.0601], [0.9203, 1.0, 0.0474], [0.933, 1.0, 0.0348], [0.9456, 0.9884, 0.0221], [0.9583, 0.9739, 0.0095], [0.9709, 0.9593, 0.0], [0.9836, 0.9448, 0.0], [0.9962, 0.9303, 0.0], [1.0, 0.9158, 0.0], [1.0, 0.9012, 0.0], [1.0, 0.8867, 0.0], [1.0, 0.8722, 0.0], [1.0, 0.8577, 0.0], [1.0, 0.8431, 0.0], [1.0, 0.8286, 0.0], [1.0, 0.8141, 0.0], [1.0, 0.7996, 0.0], [1.0, 0.785, 0.0], [1.0, 0.7705, 0.0], [1.0, 0.756, 0.0], [1.0, 0.7415, 0.0], [1.0, 0.7269, 0.0], [1.0, 0.7124, 0.0], [1.0, 0.6979, 0.0], [1.0, 0.6834, 0.0], [1.0, 0.6688, 0.0], [1.0, 0.6543, 0.0], [1.0, 0.6398, 0.0], [1.0, 0.6253, 0.0], [1.0, 0.6107, 0.0], [1.0, 0.5962, 0.0], [1.0, 0.5817, 0.0], [1.0, 0.5672, 0.0], [1.0, 0.5527, 0.0], [1.0, 0.5381, 0.0], [1.0, 0.5236, 0.0], [1.0, 0.5091, 0.0], [1.0, 0.4946, 0.0], [1.0, 0.48, 0.0], [1.0, 0.4655, 0.0], [1.0, 0.451, 0.0], [1.0, 0.4365, 0.0], [1.0, 0.4219, 0.0], [1.0, 0.4074, 0.0], [1.0, 0.3929, 0.0], [1.0, 0.3784, 0.0], [1.0, 0.3638, 0.0], [1.0, 0.3493, 0.0], [1.0, 0.3348, 0.0], [1.0, 0.3203, 0.0], [1.0, 0.3057, 0.0], [1.0, 0.2912, 0.0], [1.0, 0.2767, 0.0], [1.0, 0.2622, 0.0], [1.0, 0.2476, 0.0], [1.0, 0.2331, 0.0], [1.0, 0.2186, 0.0], [1.0, 0.2041, 0.0], [1.0, 0.1895, 0.0], [1.0, 0.175, 0.0], [1.0, 0.1605, 0.0], [1.0, 0.146, 0.0], [1.0, 0.1314, 0.0], [1.0, 0.1169, 0.0], [1.0, 0.1024, 0.0], [1.0, 0.0879, 0.0], [0.9991, 0.0733, 0.0], [0.9813, 0.0588, 0.0], [0.9635, 0.0443, 0.0], [0.9456, 0.0298, 0.0], [0.9278, 0.0153, 0.0], [0.91, 0.0007, 0.0], [0.8922, 0.0, 0.0], [0.8743, 0.0, 0.0], [0.8565, 0.0, 0.0], [0.8387, 0.0, 0.0], [0.8209, 0.0, 0.0], [0.803, 0.0, 0.0], [0.7852, 0.0, 0.0], [0.7674, 0.0, 0.0], [0.7496, 0.0, 0.0], [0.7317, 0.0, 0.0], [0.7139, 0.0, 0.0], [0.6961, 0.0, 0.0], [0.6783, 0.0, 0.0], [0.6604, 0.0, 0.0], [0.6426, 0.0, 0.0], [0.6248, 0.0, 0.0], [0.607, 0.0, 0.0], [0.5891, 0.0, 0.0], [0.5713, 0.0, 0.0], [0.5535, 0.0, 0.0], [0.5357, 0.0, 0.0], [0.5178, 0.0, 0.0], [0.5, 0.0, 0.0]];

const datasets = {
    //0: ["BrainWeb 05 - 3T", "./3t/05", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    //1: ["BrainWeb 05 - 1.5T", "./1.5T/05", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    2: ["BrainWeb 54 - 3T + Na", "./3t/54", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    3: ["BrainWeb 54 - 1.5T", "./1.5T/54", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    4: ["BrainWeb colin27 - 3T + Na", "./3t/bw", "https://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27Highres"],
    5: ["BrainWeb colin27 - 1.5T", "./1.5T/bw", "https://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27Highres"],
    6: ["Phantomag - 1.5T", "./1.5T/phantomag", "http://lab.ibb.cnr.it/Phantomag_Desc.htm"],
    7: ["Phantomag - 1T", "./1t/phantomag", "http://lab.ibb.cnr.it/Phantomag_Desc.htm"],
}

function QueryableWorker(url) {
    var instance = this,
        worker = new Worker(url //, {type: "module"}
        ),
        listeners = {};
    
    worker.onerror = function (event) {
        console.log(event.message, event);
    };

    this.defaultListener = function () {};

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
const w = new QueryableWorker("./intensityCalculations.js");

//var xdim = 256;
//var ydim = 256;
//var zdim = 256;
var array_pd, array_t1, array_t2, array_t2s, array_na_mm, array_na_t1, array_na_t2s, array_na_t2f, result, k_data_im_re, slice_data;

const na_tabs = ["params-Na-tab", "params-SQ-tab", "params-TQ-tab", "params-TQF-tab", "params-TQSQR-tab"];
const h_tabs = ["params-IR-tab", "params-SE-tab", "params-bSSFP-tab", "params-FISP-tab", "params-PSIF-tab", "params-FLASH-tab", "params-SGRE-tab"];
var current_tab = "params-IR-tab";
var selected_tab = "params-IR-tab";

function sequenceParametersKeyDown(e) {
    if(e.code == "Enter") {
        startScan();
    }
}


const loadDataMessageHandler = function (data) {
    array_pd = data[0];
    array_t1 = data[1];
    array_t2 = data[2];
    array_t2s = data[3];
    array_na_mm = data[4];
    array_na_t1 = data[5];
    array_na_t2s = data[6];
    array_na_t2f = data[7];
    var zdim = data[8];
    var ydim = data[9];
    var xdim = data[10];

    slice = document.getElementById("slice");
    slice.max = zdim;
    slice.value = Math.round(zdim/2);

    slice = document.getElementById("or_slice");
    slice.max = zdim;
    slice.value = Math.round(zdim/2);

    slice = document.getElementById("xdim")
    slice.max = xdim;
    slice.value = xdim;
    slice = document.getElementById("ydim")
    slice.max = ydim;
    slice.value = ydim;
    slice = document.getElementById("zdim")
    slice.max = zdim;
    slice.value = zdim;

    r = document.getElementById("content");
    r.classList.remove("hidden");
    spin = document.getElementById("datasetLoading");
    spin.classList.add("hidden");
    displayDataSet();
    setSequence("IR");
    if (array_na_mm == null) {
        for (var tab in na_tabs) {
            t = document.getElementById(tab);
            if (t != null) {
                t.classList.add("hidden");
            }
        }
    } else {
        for (var tab in na_tabs) {
            t = document.getElementById(tab);
            if (t!=null) {
                t.classList.remove("hidden");
            }
        }
    }
};
w.addListener('loadData', loadDataMessageHandler);

var imgResult;
var kResult;
const resultMessageHandler = function (data) {
    //console.log("image result");
    if (data[0] != undefined) {
        imgResult = data[0];
        //console.log(imgResult.data.reduce((a,b)=>Math.max(a,b), 0));
        if (data[1] != undefined) {
            kResult = data[1];
            if (data[2] != undefined) {
                var kSpace_filt = data[2];
                setKSpaceFilt(...kSpace_filt);
            } else {
                setKSpaceFilt(256,256,0,256);
            }
        }
    }
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    r.classList.remove("hidden");
    spin.classList.add("hidden");
    document.getElementById("kspace").removeAttribute("disabled");

    var windowing = document.getElementById("windowing");
    var colorBar = document.getElementById("colorBarContainer");
    if(isCurrentTabNa()) {
        windowing.classList.add("hidden");
        plot_colormap("colorBar");
        colorBar.classList.remove("hidden");
    } else {
        windowing.classList.remove("hidden");
        colorBar.classList.add("hidden");
    }
    
    
    var xdim = imgResult.xdim;
    var ydim = imgResult.ydim;
    var zdim = imgResult.zdim;

    var canvas = document.getElementById("imgResult");
    canvas.width = xdim;
    canvas.height = ydim;

    var slice = document.getElementById("r_slice");
    slice.max = zdim-1;
    slice.value = Math.round(zdim/2);

    slice = document.getElementById("or_slice");
    slice.max = zdim-1;
    slice.value = Math.round(zdim/2);

    slice = document.getElementById("xdim");
    slice.value = xdim;
    slice = document.getElementById("ydim");
    slice.value = ydim;
    slice = document.getElementById("zdim");
    slice.value = zdim;

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

function autoWindow2() {
    var min = imgResult[0];
    var max = imgResult[0];
    for (var i = 0; i < imgResult.length; i++) {
        if (imgResult[i] > max) max = imgResult[i];
        if (imgResult[i] < min) min = imgResult[i];
    }
    var ww = 4096 * max - 4096 * min;
    var wc = 4096 * min + ww * 0.5;
    document.getElementById("wc").value = wc;
    document.getElementById("ww").value = ww;
    displayAndWindow3DImage();
}

function autoWindow() {
    var histo = new Int32Array(2048);
    for (var i = 0; i < histo.length; i++) {
        histo[i] = 0;
    }
    for (var i = 0; i < imgResult.length; i++) {
        histo[Math.floor(imgResult[i] * histo.length)] += 1;
    }
    var thresh = 0.9 * histo.reduce((p, c) => p + c, -histo[0]);
    var wc = 1;
    var ww = 0;
    var max = histo[1];
    for (var i = 1; i < histo.length; i++) {
        if (histo[i] > max) {
            max = histo[i]
            wc = i;
        }
    }
    var sum = histo[wc];
    for (ww = 1; ww < histo.length; ww++) {
        if (ww % 2 == 0) {
            if (wc + Math.ceil(ww * 0.5) < histo.length) {
                sum += histo[wc + Math.ceil(ww * 0.5)];
            }
        } else {
            if (wc - Math.ceil(ww * 0.5) > 0) {
                sum += histo[wc - Math.ceil(ww * 0.5)];
            }
        }
        if (sum > thresh) {
            break;
        }
    }
    var start = Math.max(1, Math.ceil(wc - ww * 0.5));
    var end = Math.min(histo.length, Math.ceil(wc + ww * 0.5));
    ww = end - start;
    wc = start + Math.floor(ww * 0.5);
    document.getElementById("wc").value = (4096 / histo.length) * wc;
    document.getElementById("ww").value = (4096 / histo.length) * ww;
    displayAndWindow3DImage();
}

function scrollDataset(event) {
    slice = document.getElementById("slice")
    slice.value = slice.valueAsNumber - Math.sign(event.deltaY) * Math.max(1, Math.abs(Math.round(event.deltaY / 100)))
    if ("createEvent" in document) {
        evt = new Event("change", {"bubbles": false, "canceable": true});
        slice.dispatchEvent(evt);
    } else
        slice.fireEvent("onchange");
    event.preventDefault();
}

function scrollResult(event) {
    if (isCurrentTabNa()) {
        r_slice = document.getElementById("or_slice");
    } else {
        r_slice = document.getElementById("r_slice")
    }
    r_slice.value = r_slice.valueAsNumber - Math.sign(event.deltaY) * Math.max(1, Math.abs(Math.round(event.deltaY / 100)))
    if ("createEvent" in document) {
        evt = new Event("change", {"bubbles": false, "cancelable": true});
        r_slice.dispatchEvent(evt);
    } else
        r_slice.fireEvent("onchange");
    event.preventDefault();
}

var windowing = false;

function startWindowing() {
    windowing = true;
}

function endWindowing() {
    windowing = false;
}

function windowResult(event) {
    if (windowing) {
        ww = document.getElementById("ww")
        wc = document.getElementById("wc")

        ww.value = ww.valueAsNumber - event.movementY * 4
        wc.value = wc.valueAsNumber - event.movementX * 4

        evt = new Event("change", {"bubbles": false, "canceable": true});
        ww.dispatchEvent(evt);
    }
}

function isCurrentTabNa() {
    return na_tabs.includes(current_tab);
}

function displayAndWindow3DImage() {
    //canvas = document.createElement('canvas');
    //res_element = document.getElementById("result");
    var canvas = document.getElementById("imgResult");
    
    var xdim = imgResult.xdim;
    var ydim = imgResult.ydim;
    var zdim = imgResult.zdim;

    ctx = canvas.getContext('2d');
    //canvas.width = xdim;
    //canvas.height = ydim;
    idata = ctx.createImageData(xdim, ydim);

    var image_result = new Uint8ClampedArray(xdim * ydim * 4);
    if(isCurrentTabNa()) {
        in_slice = document.getElementById("or_slice");
    } else {
        in_slice = document.getElementById("r_slice");
    }
    slice = parseInt(in_slice.value);
    
    if (isCurrentTabNa()) {
        for (var x = 0; x < xdim * ydim; x++) {
            var val = Math.round(imgResult.data[x + slice * xdim * ydim] * 255);
            if (val > 255) { val = 255; }
            if (val < 0) {val = 0;}
            var c = jetmap[val];
            image_result[4 * x] = c[0]*255
            image_result[4 * x + 1] = c[1]*255
            image_result[4 * x + 2] = c[2]*255
            image_result[4 * x + 3] = 255
        }
    } else {
        in_ww = document.getElementById("ww");
        in_wc = document.getElementById("wc");
        var ww = parseFloat(in_ww.value) * 0.5
        var wc = parseFloat(in_wc.value)
        for (var x = 0; x < xdim * ydim; x++) {
            val = imgResult.data[x + slice * xdim * ydim] * 4096
            if (val <= (wc - ww)) {
                val = 0
            } else if (val >= (wc + ww)) {
                val = 255;
            } else {
                val = 255 * (val - (wc - ww)) / (ww)
            }
            image_result[4 * x] = val
            image_result[4 * x + 1] = val
            image_result[4 * x + 2] = val
            image_result[4 * x + 3] = 255
        }
    }


    idata.data.set(image_result);
    ctx.putImageData(idata, 0, 0);

    k_canvas = document.getElementById("kResult");
    var xdim = 256;
    var ydim = 256;
    k_canvas.width = xdim;
    k_canvas.height = ydim
    k_ctx = k_canvas.getContext('2d');
    kdata = k_ctx.createImageData(xdim, ydim);

    var mult = document.getElementById("kspacemult").valueAsNumber;
    var xlines = document.getElementById("k_xline").valueAsNumber;
    var ylines = document.getElementById("k_yline").valueAsNumber;
    var fmin = document.getElementById("k_fmin").valueAsNumber;
    var fmax = document.getElementById("k_fmax").valueAsNumber;
    var dx = xdim / xlines;
    var dy = ydim / ylines;
    k_result = new Uint8ClampedArray(xdim * ydim * 4);
    var zend = 10;
    for (var index = 0; index < xdim * ydim; index++) {
        val = kResult[index + slice * xdim * ydim] * mult
        var _fmax = fmax;
        if (slice<zend) {
            _fmax = (slice/zend)*(slice/zend) * fmax;
        }
        if (slice > zdim-zend) {
            _fmax = (zdim-slice)/zend*(zdim-slice)/zend * fmax;
        }
        var y = Math.floor(index / xdim);
        var x = index % xdim;
        var f = Math.sqrt((x - xdim / 2) * (x - xdim / 2) + (y - ydim / 2) * (y - ydim / 2));
        var res1 = f >= fmin && f <= _fmax;
        var res = res1 && ((Math.floor(x / dx) * dx - x) < 1 && (Math.floor(x / dx) * dx - x) > -1) && ((Math.floor(y / dy) * dy - y) < 1 && (Math.floor(y / dy) * dy - y) > -1);

        if (!res) {
            k_result[4 * index] = 255
            k_result[4 * index + 1] = 0
            k_result[4 * index + 2] = 0
            k_result[4 * index + 3] = 255
        } else {
            k_result[4 * index] = val
            k_result[4 * index + 1] = val
            k_result[4 * index + 2] = val
            k_result[4 * index + 3] = 255
        }
    }

    kdata.data.set(k_result);
    k_ctx.putImageData(kdata, 0, 0);
}

function fillDatasets() {
    sel = document.getElementById("datasetPath");
    for(var p in datasets) {
        op = document.createElement("option");
        op.value = p;
        op.innerText = datasets[p][0];
        sel.appendChild(op);
    }
    changedDataset();
}

function changedDataset() {
    url = datasets[document.getElementById("datasetPath").value][2];
    a = document.getElementById("datasetURL")
    a.setAttribute("href", url);
    a.innerText = "Dataset source";
}

function loadFuzzyDataSet() {
    content = document.getElementById("content");
    content.classList.add("hidden");
    r = document.getElementById("result");
    r.classList.add("hidden");

    path = datasets[document.getElementById("datasetPath").value][1];
    return new Promise(async (resolve, reject) => {

        spin = document.getElementById("datasetLoading");
        spin.classList.remove("hidden");

        w.sendQuery("loadData", path);
        resolve();
    });
}

function plot_colormap(canvas_id) {
    let canvas = document.getElementById(canvas_id);
    let ctx = canvas.getContext("2d");
    for (let x = 0; x < 256; x++) {
      let color = jetmap[x];
      let r = color[0]*256;
      let g = color[1]*256;
      let b = color[2]*256;
      ctx.fillStyle = 'rgb(' + r + ',' + g + ',' + b + ')';
      //ctx.fillRect(0, x * canvas.height / 256, canvas.width, canvas.height / 256);
      ctx.fillRect(Math.floor(x*(canvas.width/256)), 0, Math.ceil(canvas.width / 256), canvas.height);
    }
}

function setTabs() {
    var tabId = "params-" + selectedSequence;
    var tabHeadId = "params-" + selectedSequence + "-tab";
    selected_tab = tabHeadId;
    elems = document.getElementById("sequence").getElementsByClassName("nav-link")
    for(var x=0;x<elems.length;x++) { 
        if(elems[x].id == tabHeadId) {
            elems[x].classList.add("active");
        } else {
            elems[x].classList.remove("active");
        }
    }
    elems = document.getElementById("sequence").getElementsByClassName("tab-pane")
    for(var x=0;x<elems.length;x++) { 
        if(elems[x].id == tabId) {
            elems[x].classList.add("active", "show");
        }
        else {
            elems[x].classList.remove("active", "show");
        }
    }
}

var selectedSequence = "IR";

function setSequence(sequence) {
    selectedSequence = sequence;
    setTabs();
    updateTime();
}

function setKSpaceFilt(xlines, ylines, fmin, fmax) {
    document.getElementById("k_xline").value = xlines;
    document.getElementById("k_yline").value = ylines;
    document.getElementById("k_fmin").value = fmin;
    document.getElementById("k_fmax").value = fmax;
    document.getElementById("k_xline_number").value = xlines;
    document.getElementById("k_yline_number").value = ylines;
    document.getElementById("k_fmin_number").value = fmin;
    document.getElementById("k_fmax_number").value = fmax;
}

function reco(update_slider, noIfft = false) {
    document.getElementById("kspace").setAttribute("disabled","disabled");
    var xlines, ylines, fmin, fmax;
    if (update_slider) {
        xlines = document.getElementById("k_xline_number").valueAsNumber;
        ylines = document.getElementById("k_yline_number").valueAsNumber;
        fmin = document.getElementById("k_fmin_number").valueAsNumber;
        fmax = document.getElementById("k_fmax_number").valueAsNumber;
        document.getElementById("k_xline").value = xlines;
        document.getElementById("k_yline").value = ylines;
        document.getElementById("k_fmin").value = fmin;
        document.getElementById("k_fmax").value = fmax;
    } else {
        xlines = document.getElementById("k_xline").valueAsNumber;
        ylines = document.getElementById("k_yline").valueAsNumber;
        fmin = document.getElementById("k_fmin").valueAsNumber;
        fmax = document.getElementById("k_fmax").valueAsNumber;
        document.getElementById("k_xline_number").value = xlines;
        document.getElementById("k_yline_number").value = ylines;
        document.getElementById("k_fmin_number").value = fmin;
        document.getElementById("k_fmax_number").value = fmax;
    }
    if (!noIfft) {
        var param_div = document.getElementById("params-general");
        var params = read_params(param_div, {});
        var xdim = Math.round(params["xdim"]);
        xdim = xdim > 0 ? xdim : k_xdim;
        xdim = xdim > k_xdim ? k_xdim : xdim;
        var ydim = Math.round(params["ydim"]);
        ydim = ydim > 0 ? ydim : k_ydim;
        ydim = ydim > k_ydim ? k_ydim : ydim;
        var zdim = Math.round(params["zdim"]);
        zdim = zdim > 0 ? zdim : k_zdim;
        zdim = zdim > k_zdim ? k_zdim : zdim;
        var params = [xdim, ydim, zdim, xlines, ylines, fmin, fmax, noIfft];
        //console.log(params);
        w.sendQuery("reco", params);
    } else {
        document.getElementById("kspace").removeAttribute("disabled");
        displayAndWindow3DImage();
    }
}

function read_params(param_div, params) {
    for (var child_id in param_div.children) {
        var child = param_div.children[child_id];
        if (child.children == undefined) { continue;}
        for(var input_id in child.children) {
            var input = child.children[input_id];
            switch (input.type) {
                case "number":
                case "range":
                    params[input.name] = parseFloat(input.value);
                    break;
                case "select-one":
                    params[input.name] = input.value;
                    break;
            }
        }
    }
    return params;
}

function startScan() {
    current_tab = selected_tab;
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var params = {sequence: selectedSequence};

    var param_div = document.getElementById("params-general");
    params = read_params(param_div, params);

    param_div = document.getElementById("params-" + selectedSequence);
    params = read_params(param_div, params);
    
    if(selectedSequence == "TQSQR") {
        var param_div = document.getElementById("params-SQ");
        var tq_params = {};
        tq_params = read_params(param_div, tq_params);
        var param_div = document.getElementById("params-TQ");
        var sq_params = {};
        sq_params = read_params(param_div, sq_params);
        params["tq_params"] = tq_params;
        params["sq_params"] = sq_params;
    }
    //console.log(params);
    w.sendQuery("simulateImage", params);
}

function displayDataSet() {
    //display3DImage(document.getElementById("imgPD"), array_pd);
    //display3DImage(document.getElementById("imgT1"), array_t1);
    //display3DImage(document.getElementById("imgT2"), array_t2);
    //display3DImage(document.getElementById("imgT2s"), array_t2s);
}

function updateTime() {
    var functionName = "update" + selectedSequence + "Time";
    if (this.hasOwnProperty(functionName)) {
        this[functionName]();
    }
}

function updateIRTime() {
    var ti = parseFloat(document.getElementById("ir_ti").value)
    var te = parseFloat(document.getElementById("ir_te").value)
    var tr = parseFloat(document.getElementById("ir_tr").value)
    var time = document.getElementById("ir_time");

    if(ti+te >= tr) {
        time.innerText = "TE+TI has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr*ydim*zdim);
    }
}

function updateSETime() {
    var te = parseFloat(document.getElementById("se_te").value)
    var tr = parseFloat(document.getElementById("se_tr").value)
    var time = document.getElementById("se_time");
    
    if(te >= tr) {
        time.innerText = "TE has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr*ydim*zdim);
    }
}

function updatePSIFTime() {
    var te = parseFloat(document.getElementById("psif_te").value)
    var tr = parseFloat(document.getElementById("psif_tr").value)
    var fa = parseFloat(document.getElementById("psif_fa").value)
    var time = document.getElementById("psif_time");
    
    if(fa>180) { document.getElementById("psif_fa").value = 180; }
    if(fa<-180) { document.getElementById("psif_fa").value = -180; }

    time.innerText = formatTime(tr*ydim*zdim);
}

function updateFISPTime() {
    var te = parseFloat(document.getElementById("fisp_te").value)
    var tr = parseFloat(document.getElementById("fisp_tr").value)
    var fa = parseFloat(document.getElementById("fisp_fa").value)
    var time = document.getElementById("fisp_time");
        
    if(fa>180) { document.getElementById("fisp_fa").value = 180; }
    if(fa<-180) { document.getElementById("fisp_fa").value = -180; }

    time.innerText = formatTime(te*2*ydim*zdim);
}

function updateSGRETime() {
    var te = parseFloat(document.getElementById("sgre_te").value)
    var tr = parseFloat(document.getElementById("sgre_tr").value)
    var fa = parseFloat(document.getElementById("sgre_fa").value)
    var time = document.getElementById("sgre_time");
        
    if(fa>180) { document.getElementById("sgre_fa").value = 180; }
    if(fa<-180) { document.getElementById("sgre_fa").value = -180; }

    if(te >= tr) {
        time.innerText = "TE has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr*ydim*zdim);
    }
}

function updateBalancedSSFPTime() {
    var te = parseFloat(document.getElementById("bssfp_te").value)
    var fa = parseFloat(document.getElementById("bssfp_fa").value)
    var time = document.getElementById("bssfp_time");
        
    if(fa>180) { document.getElementById("bssfp_fa").value = 180; }
    if(fa<-180) { document.getElementById("bssfp_fa").value = -180; }

    document.getElementById("bssfp_tr").value = te*2;

    time.innerText = formatTime(te*2*ydim*zdim);
}

function updateNaTime() {
    var te = parseFloat(document.getElementById("na_te").value)
    var tr = parseFloat(document.getElementById("na_tr").value)
    var time = document.getElementById("na_time");
    
    if(te >= tr) {
        time.innerText = "TE has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr*ydim*zdim);
    }
}

function updateSQTime() {
    var te_start = parseFloat(document.getElementById("sq_te_start").value)
    var te_end = parseFloat(document.getElementById("sq_te_end").value)
    var te_step = parseFloat(document.getElementById("sq_te_step").value)
    var te = document.getElementById("sq_te");

    var tes = ""+te_start;
    for (var t=te_start+te_step;t<=te_end;t+=te_step) {
        tes += ", " + t;
    }
    te.innerText = tes;
}

function updateTQTime() {
    var te_start = parseFloat(document.getElementById("tq_te_start").value)
    var te_end = parseFloat(document.getElementById("tq_te_end").value)
    var te_step = parseFloat(document.getElementById("tq_te_step").value)
    var te = document.getElementById("tq_te");

    var tes = ""+te_start;
    for (var t=te_start+te_step;t<=te_end;t+=te_step) {
        tes += ", " + t;
    }
    te.innerText = tes;
}

function updateTQSQRTime() {
}

function updateTQFTime() {
}

function formatTime(time) {
    d = new Date(time);
    timeString = "";
    if(d.getFullYear() > 1970) {
        timeString += (d.getFullYear()-1970) + " years ";
    }
    if(d.getDate() > 1) {
        timeString += (d.getDate()-1) + " days ";
    }
    timeString += (d.getHours()<11?"0":"") + (d.getHours()-1) + ":" + (d.getMinutes()<10?"0":"") + d.getMinutes() + ":" + (d.getSeconds()<10?"0":"") + d.getSeconds();
    return timeString;
}
