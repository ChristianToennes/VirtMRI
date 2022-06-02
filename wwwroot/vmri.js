const datasets = {
    //0: ["BrainWeb 05 - 3T", "./3t/05", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    //1: ["BrainWeb 05 - 1.5T", "./1.5T/05", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    2: ["BrainWeb 54 - 3T + Na", "./3t/54", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    3: ["BrainWeb 54 - 1.5T + Na", "./1.5T/54", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    4: ["BrainWeb colin27 - 3T + Na", "./3t/bw", "https://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27Highres"],
    5: ["BrainWeb colin27 - 1.5T + Na", "./1.5T/bw", "https://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27Highres"],
    6: ["Phantomag - 1.5T", "./1.5T/phantomag", "http://lab.ibb.cnr.it/Phantomag_Desc.htm"],
    7: ["Phantomag - 1T", "./1t/phantomag", "http://lab.ibb.cnr.it/Phantomag_Desc.htm"],
}

function QueryableWorker(url) {
    var instance = this,
        worker = new Worker(url //, {type: "module"}
        ),
        listeners = {};

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

var xdim = 256;
var ydim = 256;
var zdim = 256;
var array_pd, array_t1, array_t2, array_t2s, array_na_mm, result, k_data_im_re, slice_data;

const na_tabs = ["params-sq-tab", "params-tq-tab"];
const h_tabs = ["params-ir-tab", "params-se-tab", "params-bssfp-tab", "params-fisp-tab", "params-psif-tab", "params-flash-tab", "params-sgre-tab", "params-sq-tab", "params-tq-tab"];

const loadDataMessageHandler = function (data) {
    array_pd = data[0];
    array_t1 = data[1];
    array_t2 = data[2];
    array_t2s = data[3];
    //array_na_t2s = data[4];
    //array_na_t2f = data[5];
    array_na_mm = data[4];
    zdim = data[5];
    ydim = data[6];
    xdim = data[7];

    slice = document.getElementById("slice");
    slice.max = zdim;
    slice.value = Math.round(zdim/2);

    r = document.getElementById("content");
    r.classList.remove("hidden");
    spin = document.getElementById("datasetLoading");
    spin.classList.add("hidden");
    displayDataSet();
    setInversionRecovery();
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
    }
    if (data[1] != undefined) {
        kResult = data[1];
    }
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    r.classList.remove("hidden");
    spin.classList.add("hidden");
    document.getElementById("kspace").removeAttribute("disabled");
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
    r_slice = document.getElementById("r_slice")
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

        if ("createEvent" in document) {
            var evt = document.createEvent("HTMLEvents");
            evt.initEvent("change", false, true);
            ww.dispatchEvent(evt);
        } else {
            ww.fireEvent("onchange");
        }
    }
}

function displayAndWindow3DImage() {
    //canvas = document.createElement('canvas');
    //res_element = document.getElementById("result");
    canvas = document.getElementById("imgResult");

    ctx = canvas.getContext('2d');
    canvas.width = xdim;
    canvas.height = ydim;
    idata = ctx.createImageData(xdim, ydim);

    in_ww = document.getElementById("ww");
    in_wc = document.getElementById("wc");
    in_slice = document.getElementById("r_slice");
    slice = parseInt(in_slice.value);

    image_result = new Uint8ClampedArray(xdim * ydim * 4);
    var ww = parseFloat(in_ww.value) * 0.5
    var wc = parseFloat(in_wc.value)

    for (var x = 0; x < xdim * ydim; x++) {
        val = imgResult[x + slice * xdim * ydim] * 4096
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

    idata.data.set(image_result);
    ctx.putImageData(idata, 0, 0);

    k_canvas = document.getElementById("kResult");
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
    for (var index = 0; index < xdim * ydim; index++) {
        val = kResult[index + slice * xdim * ydim] * mult

        var y = Math.floor(index / xdim);
        var x = index % xdim;
        var f = Math.sqrt((x - xdim / 2) * (x - xdim / 2) + (y - ydim / 2) * (y - ydim / 2));
        var res1 = f >= fmin && f <= fmax;
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

function setTabs(tabId, tabHeadId) {
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

var selectedSequence = inversionRecovery;

function setSpinEcho() {
    setTabs("params-se", "params-se-tab");
    updateSETime();
    selectedSequence = spinEcho;
}

function spinEcho() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var te = parseFloat(document.getElementById("se_te").value)
    var tr = parseFloat(document.getElementById("se_tr").value)

    w.sendQuery("spinEcho", te, tr);
}

function setInversionRecovery() {
    setTabs("params-ir", "params-ir-tab");
    updateIRTime();
    selectedSequence = inversionRecovery;
}

function inversionRecovery() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var ti = parseFloat(document.getElementById("ir_ti").value)
    var tr = parseFloat(document.getElementById("ir_tr").value)
    var te = parseFloat(document.getElementById("ir_te").value)

    w.sendQuery("inversionRecovery", te, tr, ti);
}

function setFlash() {
    setTabs("params-flash", "params-flash-tab");
    updateFlashTime();
    selectedSequence = flash;
}

function flash() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var te = parseFloat(document.getElementById("flash_te").value)
    var tr = parseFloat(document.getElementById("flash_tr").value)
    var fa = parseFloat(document.getElementById("flash_fa").value)

    w.sendQuery("flash", te, tr, fa);
}

function setBalancedSSFP() {
    setTabs("params-bssfp", "params-bssfp-tab");
    updateBalancedSSFPTime();
    selectedSequence = balancedSSFP;
}

function balancedSSFP() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var te = parseFloat(document.getElementById("fisp_te").value)
    var fa = parseFloat(document.getElementById("fisp_fa").value)
    var tr = te*2;

    w.sendQuery("fisp", te, tr, fa);
}

function setFISP() {
    setTabs("params-fisp", "params-fisp-tab");
    updateFISPTime();
    selectedSequence = FISP;
}

function FISP() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var te = parseFloat(document.getElementById("fisp_te").value)
    var fa = parseFloat(document.getElementById("fisp_fa").value)
    var tr = te*2;

    w.sendQuery("fisp", te, tr, fa);
}

function setPSIF() {
    setTabs("params-psif", "params-psif-tab");
    updatePSIFTime();
    selectedSequence = PSIF;
}

function PSIF() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var te = parseFloat(document.getElementById("psif_te").value)
    var fa = parseFloat(document.getElementById("psif_fa").value)
    var tr = te*2;

    w.sendQuery("psif", te, tr, fa);
}

function setSGRE() {
    setTabs("params-sgre", "params-sgre-tab");
    updateSGRETime();
    selectedSequence = SGRE;
}

function SGRE() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var te = parseFloat(document.getElementById("sgre_te").value)
    var fa = parseFloat(document.getElementById("sgre_fa").value)
    var tr = parseFloat(document.getElementById("sgre_tr").value)

    w.sendQuery("sgre", te, tr, fa);
}

function setSQ() {
    setTabs("params-sq", "params-sq-tab");
    updateSQTime();
    selectedSequence = SQ;
}

function SQ() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var te_start = parseFloat(document.getElementById("sq_te_start").value)
    var te_end = parseFloat(document.getElementById("sq_te_end").value)
    var te_step = parseFloat(document.getElementById("sq_te_step").value)

    w.sendQuery("sq", te_start, te_end, te_step);
}

function setTQ() {
    setTabs("params-tq", "params-tq-tab");
    updateTQTime();
    selectedSequence = TQ;
}

function TQ() {
    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    slice = document.getElementById("r_slice");
    slice.max = zdim;
    r.classList.add("hidden");
    spin.classList.remove("hidden");

    var te_start = parseFloat(document.getElementById("tq_te_start").value)
    var te_end = parseFloat(document.getElementById("tq_te_end").value)
    var te_step = parseFloat(document.getElementById("tq_te_step").value)

    w.sendQuery("tq", te_start, te_end, te_step);
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
    w.sendQuery("reco", xlines, ylines, fmin, fmax, noIfft);
}

function startScan() {
    selectedSequence();
}

function displayDataSet() {
    display3DImage(document.getElementById("imgPD"), array_pd);
    display3DImage(document.getElementById("imgT1"), array_t1);
    display3DImage(document.getElementById("imgT2"), array_t2);
    display3DImage(document.getElementById("imgT2s"), array_t2s);
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

    if(te >= tr) {
        time.innerText = "TE has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr*ydim*zdim);
    }
}

function updateFlashTime() {
    var te = parseFloat(document.getElementById("flash_te").value)
    var tr = parseFloat(document.getElementById("flash_tr").value)
    var fa = parseFloat(document.getElementById("flash_fa").value)
    var time = document.getElementById("flash_time");
        
    if(fa>180) { document.getElementById("flash_fa").value = 180; }
    if(fa<-180) { document.getElementById("flash_fa").value = -180; }

    if(te >= tr) {
        time.innerText = "TE has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr*ydim*zdim);
    }
}

function updateFISPTime() {
    var te = parseFloat(document.getElementById("fisp_te").value)
    var tr = parseFloat(document.getElementById("fisp_tr").value)
    var fa = parseFloat(document.getElementById("fisp_fa").value)
    var time = document.getElementById("fisp_time");
        
    if(fa>180) { document.getElementById("fisp_fa").value = 180; }
    if(fa<-180) { document.getElementById("fisp_fa").value = -180; }

    if(te >= tr) {
        time.innerText = "TE has to be smaller than TR";
    } else {
        time.innerText = formatTime(te*2*ydim*zdim);
    }
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
