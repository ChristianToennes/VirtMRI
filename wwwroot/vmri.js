const jetmap = [ [0.0, 0.0, 0.5], [0.0, 0.0, 0.5178], [0.0, 0.0, 0.5357], [0.0, 0.0, 0.5535], [0.0, 0.0, 0.5713], [0.0, 0.0, 0.5891], [0.0, 0.0, 0.607], [0.0, 0.0, 0.6248], [0.0, 0.0, 0.6426], [0.0, 0.0, 0.6604], [0.0, 0.0, 0.6783], [0.0, 0.0, 0.6961], [0.0, 0.0, 0.7139], [0.0, 0.0, 0.7317], [0.0, 0.0, 0.7496], [0.0, 0.0, 0.7674], [0.0, 0.0, 0.7852], [0.0, 0.0, 0.803], [0.0, 0.0, 0.8209], [0.0, 0.0, 0.8387], [0.0, 0.0, 0.8565], [0.0, 0.0, 0.8743], [0.0, 0.0, 0.8922], [0.0, 0.0, 0.91], [0.0, 0.0, 0.9278], [0.0, 0.0, 0.9456], [0.0, 0.0, 0.9635], [0.0, 0.0, 0.9813], [0.0, 0.0, 0.9991], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.002, 1.0], [0.0, 0.0176, 1.0], [0.0, 0.0333, 1.0], [0.0, 0.049, 1.0], [0.0, 0.0647, 1.0], [0.0, 0.0804, 1.0], [0.0, 0.0961, 1.0], [0.0, 0.1118, 1.0], [0.0, 0.1275, 1.0], [0.0, 0.1431, 1.0], [0.0, 0.1588, 1.0], [0.0, 0.1745, 1.0], [0.0, 0.1902, 1.0], [0.0, 0.2059, 1.0], [0.0, 0.2216, 1.0], [0.0, 0.2373, 1.0], [0.0, 0.2529, 1.0], [0.0, 0.2686, 1.0], [0.0, 0.2843, 1.0], [0.0, 0.3, 1.0], [0.0, 0.3157, 1.0], [0.0, 0.3314, 1.0], [0.0, 0.3471, 1.0], [0.0, 0.3627, 1.0], [0.0, 0.3784, 1.0], [0.0, 0.3941, 1.0], [0.0, 0.4098, 1.0], [0.0, 0.4255, 1.0], [0.0, 0.4412, 1.0], [0.0, 0.4569, 1.0], [0.0, 0.4725, 1.0], [0.0, 0.4882, 1.0], [0.0, 0.5039, 1.0], [0.0, 0.5196, 1.0], [0.0, 0.5353, 1.0], [0.0, 0.551, 1.0], [0.0, 0.5667, 1.0], [0.0, 0.5824, 1.0], [0.0, 0.598, 1.0], [0.0, 0.6137, 1.0], [0.0, 0.6294, 1.0], [0.0, 0.6451, 1.0], [0.0, 0.6608, 1.0], [0.0, 0.6765, 1.0], [0.0, 0.6922, 1.0], [0.0, 0.7078, 1.0], [0.0, 0.7235, 1.0], [0.0, 0.7392, 1.0], [0.0, 0.7549, 1.0], [0.0, 0.7706, 1.0], [0.0, 0.7863, 1.0], [0.0, 0.802, 1.0], [0.0, 0.8176, 1.0], [0.0, 0.8333, 1.0], [0.0, 0.849, 1.0], [0.0, 0.8647, 0.9962], [0.0, 0.8804, 0.9836], [0.0, 0.8961, 0.9709], [0.0095, 0.9118, 0.9583], [0.0221, 0.9275, 0.9456], [0.0348, 0.9431, 0.933], [0.0474, 0.9588, 0.9203], [0.0601, 0.9745, 0.9077], [0.0727, 0.9902, 0.895], [0.0854, 1.0, 0.8824], [0.098, 1.0, 0.8697], [0.1107, 1.0, 0.8571], [0.1233, 1.0, 0.8444], [0.136, 1.0, 0.8318], [0.1486, 1.0, 0.8191], [0.1613, 1.0, 0.8065], [0.1739, 1.0, 0.7938], [0.1866, 1.0, 0.7812], [0.1992, 1.0, 0.7685], [0.2119, 1.0, 0.7559], [0.2245, 1.0, 0.7432], [0.2372, 1.0, 0.7306], [0.2498, 1.0, 0.7179], [0.2625, 1.0, 0.7052], [0.2751, 1.0, 0.6926], [0.2878, 1.0, 0.6799], [0.3004, 1.0, 0.6673], [0.3131, 1.0, 0.6546], [0.3257, 1.0, 0.642], [0.3384, 1.0, 0.6293], [0.351, 1.0, 0.6167], [0.3637, 1.0, 0.604], [0.3763, 1.0, 0.5914], [0.389, 1.0, 0.5787], [0.4016, 1.0, 0.5661], [0.4143, 1.0, 0.5534], [0.4269, 1.0, 0.5408], [0.4396, 1.0, 0.5281], [0.4522, 1.0, 0.5155], [0.4649, 1.0, 0.5028], [0.4775, 1.0, 0.4902], [0.4902, 1.0, 0.4775], [0.5028, 1.0, 0.4649], [0.5155, 1.0, 0.4522], [0.5281, 1.0, 0.4396], [0.5408, 1.0, 0.4269], [0.5534, 1.0, 0.4143], [0.5661, 1.0, 0.4016], [0.5787, 1.0, 0.389], [0.5914, 1.0, 0.3763], [0.604, 1.0, 0.3637], [0.6167, 1.0, 0.351], [0.6293, 1.0, 0.3384], [0.642, 1.0, 0.3257], [0.6546, 1.0, 0.3131], [0.6673, 1.0, 0.3004], [0.6799, 1.0, 0.2878], [0.6926, 1.0, 0.2751], [0.7052, 1.0, 0.2625], [0.7179, 1.0, 0.2498], [0.7306, 1.0, 0.2372], [0.7432, 1.0, 0.2245], [0.7559, 1.0, 0.2119], [0.7685, 1.0, 0.1992], [0.7812, 1.0, 0.1866], [0.7938, 1.0, 0.1739], [0.8065, 1.0, 0.1613], [0.8191, 1.0, 0.1486], [0.8318, 1.0, 0.136], [0.8444, 1.0, 0.1233], [0.8571, 1.0, 0.1107], [0.8697, 1.0, 0.098], [0.8824, 1.0, 0.0854], [0.895, 1.0, 0.0727], [0.9077, 1.0, 0.0601], [0.9203, 1.0, 0.0474], [0.933, 1.0, 0.0348], [0.9456, 0.9884, 0.0221], [0.9583, 0.9739, 0.0095], [0.9709, 0.9593, 0.0], [0.9836, 0.9448, 0.0], [0.9962, 0.9303, 0.0], [1.0, 0.9158, 0.0], [1.0, 0.9012, 0.0], [1.0, 0.8867, 0.0], [1.0, 0.8722, 0.0], [1.0, 0.8577, 0.0], [1.0, 0.8431, 0.0], [1.0, 0.8286, 0.0], [1.0, 0.8141, 0.0], [1.0, 0.7996, 0.0], [1.0, 0.785, 0.0], [1.0, 0.7705, 0.0], [1.0, 0.756, 0.0], [1.0, 0.7415, 0.0], [1.0, 0.7269, 0.0], [1.0, 0.7124, 0.0], [1.0, 0.6979, 0.0], [1.0, 0.6834, 0.0], [1.0, 0.6688, 0.0], [1.0, 0.6543, 0.0], [1.0, 0.6398, 0.0], [1.0, 0.6253, 0.0], [1.0, 0.6107, 0.0], [1.0, 0.5962, 0.0], [1.0, 0.5817, 0.0], [1.0, 0.5672, 0.0], [1.0, 0.5527, 0.0], [1.0, 0.5381, 0.0], [1.0, 0.5236, 0.0], [1.0, 0.5091, 0.0], [1.0, 0.4946, 0.0], [1.0, 0.48, 0.0], [1.0, 0.4655, 0.0], [1.0, 0.451, 0.0], [1.0, 0.4365, 0.0], [1.0, 0.4219, 0.0], [1.0, 0.4074, 0.0], [1.0, 0.3929, 0.0], [1.0, 0.3784, 0.0], [1.0, 0.3638, 0.0], [1.0, 0.3493, 0.0], [1.0, 0.3348, 0.0], [1.0, 0.3203, 0.0], [1.0, 0.3057, 0.0], [1.0, 0.2912, 0.0], [1.0, 0.2767, 0.0], [1.0, 0.2622, 0.0], [1.0, 0.2476, 0.0], [1.0, 0.2331, 0.0], [1.0, 0.2186, 0.0], [1.0, 0.2041, 0.0], [1.0, 0.1895, 0.0], [1.0, 0.175, 0.0], [1.0, 0.1605, 0.0], [1.0, 0.146, 0.0], [1.0, 0.1314, 0.0], [1.0, 0.1169, 0.0], [1.0, 0.1024, 0.0], [1.0, 0.0879, 0.0], [0.9991, 0.0733, 0.0], [0.9813, 0.0588, 0.0], [0.9635, 0.0443, 0.0], [0.9456, 0.0298, 0.0], [0.9278, 0.0153, 0.0], [0.91, 0.0007, 0.0], [0.8922, 0.0, 0.0], [0.8743, 0.0, 0.0], [0.8565, 0.0, 0.0], [0.8387, 0.0, 0.0], [0.8209, 0.0, 0.0], [0.803, 0.0, 0.0], [0.7852, 0.0, 0.0], [0.7674, 0.0, 0.0], [0.7496, 0.0, 0.0], [0.7317, 0.0, 0.0], [0.7139, 0.0, 0.0], [0.6961, 0.0, 0.0], [0.6783, 0.0, 0.0], [0.6604, 0.0, 0.0], [0.6426, 0.0, 0.0], [0.6248, 0.0, 0.0], [0.607, 0.0, 0.0], [0.5891, 0.0, 0.0], [0.5713, 0.0, 0.0], [0.5535, 0.0, 0.0], [0.5357, 0.0, 0.0], [0.5178, 0.0, 0.0], [0.5, 0.0, 0.0]];

const datasets = {
    //0: ["BrainWeb 05 - 3T", "./3t/05", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    //1: ["BrainWeb 05 - 1.5T", "./1.5T/05", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    2: ["BrainWeb colin27 - 3T + Na", "./3t/bw", "https://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27Highres"],
    3: ["BrainWeb colin27 - 1.5T", "./1.5T/bw", "https://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27Highres"],
    4: ["BrainWeb 54 - 3T + Na", "./3t/54", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    5: ["BrainWeb 54 - 1.5T", "./1.5T/54", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    6: ["HighRes BrainWeb colin27 - 3T + Na", "./3t/bw_2", "https://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27Highres"],
    7: ["HighRes BrainWeb colin27 - 1.5T", "./1.5T/bw_2", "https://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27Highres"],
    8: ["HighRes BrainWeb 54 - 3T + Na", "./3t/54_2", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    9: ["HighRes BrainWeb 54 - 1.5T", "./1.5T/54_2", "https://brainweb.bic.mni.mcgill.ca/anatomic_normal_20.html"],
    10: ["Phantomag - 1.5T", "./1.5T/phantomag", "http://lab.ibb.cnr.it/Phantomag_Desc.htm"],
    11: ["Phantomag - 1T", "./1t/phantomag", "http://lab.ibb.cnr.it/Phantomag_Desc.htm"],
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

const na_tabs = ["params-Na-tab", "params-SQ-tab", "params-TQ-tab", "params-TQSQR-tab"];
const h_tabs = ["params-IR-tab", "params-SE-tab", "params-bSSFP-tab", "params-FISP-tab", "params-PSIF-tab", "params-SGRE-tab"];
var current_tab = "params-IR-tab";
var selected_tab = "params-IR-tab";

var imgResultCache = {
    "cur": undefined,
    "pre": undefined,
    "cs": undefined
};

function sequenceParametersKeyDown(e) {
    if (e.code == "Enter") {
        startScan();
    }
}

const debugMessageHandler = function(data) {
    console.log(data);
};
w.addListener('debug', debugMessageHandler);

function debug_send_message(msg) {
    w.sendQuery("debug", msg);
}

const loadDataMessageHandler = function (data) {
    var array_na_mm = data[4];
    var zdim = data[8];
    var ydim = data[9];
    var xdim = data[10];

    for (var node in document.querySelectorAll('input[type="range"]')) {
        if (node != undefined && node.id != undefined && node.id.endsWith("slice")) {
            node.max = zdim;
            node.value = Math.round(zdim / 2);
        }
    }

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
    setSequence("IR");
    var na = document.getElementById("sodium_sequences");
    if (array_na_mm == null) {
        na.classList.add("hidden");
        /*for (var tab in na_tabs) {
            t = document.getElementById(na_tabs[tab]);
            if (t != null) {
                t.classList.add("hidden");
            }
        }*/
    } else {
        na.classList.remove("hidden");
        /*for (var tab in na_tabs) {
            t = document.getElementById(na_tabs[tab]);
            if (t != null) {
                t.classList.remove("hidden");
            }
        }*/
    }
};
w.addListener('loadData', loadDataMessageHandler);

var resultTimer;
const resultMessageHandler = function (data) {
    if (Array.isArray(data)) {
        imgResultCache["pre"] = data[0];
        imgResultCache["cur"] = data[1];
        imgResultCache["cs"] = data[2];
        document.getElementById("pre_label").textContent = "Full Image";
        document.getElementById("cur_label").textContent = "Undersampled Image";
        document.getElementById("pre_visible_label").textContent = "Show Full Image";
        document.getElementById("cur_visible_label").textContent = "Show Undersampled Image";
        document.getElementById("cs_visible").parentElement.style.display = "block";
    } else if (typeof (data) == "string") {
        console.log("Error occured: ", data);
    } else {
        document.getElementById("pre_label").textContent = "Previous Image Result";
        document.getElementById("cur_label").textContent = "Current Image Result";
        document.getElementById("pre_visible_label").textContent = "Show Previous Image";
        document.getElementById("cur_visible_label").textContent = "Show Current Image";
        document.getElementById("cs_visible").parentElement.style.display = "none";
        if ("cs" in imgResultCache) {
            imgResultCache["pre"] = imgResultCache["cs"];
            delete imgResultCache["cs"];
        } else if ("cur" in imgResultCache) {
            imgResultCache["pre"] = imgResultCache["cur"];
        }
        imgResultCache["cur"] = data;
    }

    console.log((performance.now() - resultTimer) / 1000);

    if (imgResultCache["cur"] != undefined) {
        document.getElementById("cur_traslice").max = imgResultCache["cur"].zdim;
        document.getElementById("cur_traslice").value = Math.round(imgResultCache["cur"].zdim / 2);
        document.getElementById("cur_sagslice").max = imgResultCache["cur"].xdim;
        document.getElementById("cur_sagslice").value = Math.round(imgResultCache["cur"].xdim / 2);
        document.getElementById("cur_corslice").max = imgResultCache["cur"].ydim;
        document.getElementById("cur_corslice").value = Math.round(imgResultCache["cur"].ydim / 2);
    }
    if (imgResultCache["pre"] != undefined) {
        document.getElementById("pre_traslice").max = imgResultCache["pre"].zdim;
        document.getElementById("pre_traslice").value = Math.round(imgResultCache["pre"].zdim / 2);
        document.getElementById("pre_sagslice").max = imgResultCache["pre"].xdim;
        document.getElementById("pre_sagslice").value = Math.round(imgResultCache["pre"].xdim / 2);
        document.getElementById("pre_corslice").max = imgResultCache["pre"].ydim;
        document.getElementById("pre_corslice").value = Math.round(imgResultCache["pre"].ydim / 2);
    }
    if (imgResultCache["cs"] != undefined) {
        document.getElementById("cs_traslice").max = imgResultCache["cs"].zdim;
        document.getElementById("cs_traslice").value = Math.round(imgResultCache["cs"].zdim / 2);
        document.getElementById("cs_sagslice").max = imgResultCache["cs"].xdim;
        document.getElementById("cs_sagslice").value = Math.round(imgResultCache["cs"].xdim / 2);
        document.getElementById("cs_corslice").max = imgResultCache["cs"].ydim;
        document.getElementById("cs_corslice").value = Math.round(imgResultCache["cs"].ydim / 2);
    }

    r = document.getElementById("result");
    spin = document.getElementById("scanningSpinner");
    var progress = document.getElementById("cs_progress").parentElement;
    if (!progress.classList.contains("hidden")) {
        progress.classList.add("hidden");
    }
    r.classList.remove("hidden");
    if (!spin.classList.contains("hidden")) {
        spin.classList.add("hidden");
    }

    var windowing = document.getElementById("cur_windowing");
    var colorBar = document.getElementById("cur_colorBarContainer");
    if (isCurrentTabNa("cur")) {
        //windowing.classList.add("hidden");
        plot_colormap("cur_colorBar");
        colorBar.classList.remove("hidden");
        if (imgResultCache["cur"] != undefined && (imgResultCache["cur"].params["sequence"] == "SQ" || imgResultCache["cur"].params["sequence"] == "Na")) {
            colorBar.classList.remove("colorBarContainerTQ");
            colorBar.classList.add("colorBarContainerSQ");
        } else {
            colorBar.classList.remove("colorBarContainerSQ");
            colorBar.classList.add("colorBarContainerTQ");
        }
    } else {
        windowing.classList.remove("hidden");
        colorBar.classList.add("hidden");
    }

    var windowing = document.getElementById("pre_windowing");
    var colorBar = document.getElementById("pre_colorBarContainer");
    if (isCurrentTabNa("pre")) {
        //windowing.classList.add("hidden");
        plot_colormap("pre_colorBar");
        colorBar.classList.remove("hidden");
        if (imgResultCache["pre"] != undefined && (imgResultCache["pre"].params["sequence"] == "SQ" || imgResultCache["pre"].params["sequence"] == "Na")) {
            colorBar.classList.remove("colorBarContainerTQ");
            colorBar.classList.add("colorBarContainerSQ");
        } else {
            colorBar.classList.remove("colorBarContainerSQ");
            colorBar.classList.add("colorBarContainerTQ");
        }
    } else {
        windowing.classList.remove("hidden");
        colorBar.classList.add("hidden");
    }

    var windowing = document.getElementById("cs_windowing");
    var colorBar = document.getElementById("cs_colorBarContainer");
    if (isCurrentTabNa("cs")) {
        //windowing.classList.add("hidden");
        plot_colormap("cs_colorBar");
        colorBar.classList.remove("hidden");
        if (imgResultCache["cs"] != undefined && (imgResultCache["cs"].params["sequence"] == "SQ" || imgResultCache["cs"].params["sequence"] == "Na")) {
            colorBar.classList.remove("colorBarContainerTQ");
            colorBar.classList.add("colorBarContainerSQ");
        } else {
            colorBar.classList.remove("colorBarContainerSQ");
            colorBar.classList.add("colorBarContainerTQ");
        }
    } else {
        windowing.classList.remove("hidden");
        colorBar.classList.add("hidden");
    }

    displayAndWindow3DImage();
    toggleImg();
};
w.addListener('result', resultMessageHandler);

const kspaceMessageHandler = function (kspace) {

};
w.addListener('kspace', kspaceMessageHandler);

const progressMessageHandler = function (progress) {
    var bar = document.getElementById("cs_progress");
    if (bar.parentElement.classList.contains("hidden")) {
        bar.parentElement.classList.remove("hidden");
    }
    bar.style.cssText = "width:" + progress + "%";
};
w.addListener('progress', progressMessageHandler);

const profileMessageHandler = function (data) {
    var [params, times] = data;

    var n = times.length
    var mean = times.reduce((a, b) => a + b) / n
    var std = Math.sqrt(times.map(x => Math.pow(x - mean, 2)).reduce((a, b) => a + b) / n)

    /*var spin = document.getElementById("scanningSpinner");
    spin.classList.add("hidden");*/

    console.log(`${params['sequence']}: ${params['xdim']}x${params['ydim']}x${params['zdim']} & ${params['nearest']==0?'Nearest':'Average'} & ${params['fft'].toUpperCase()} & ${params['img_noise_type']==0?'No':'Yes'}${params['cs']==0?'':'+CS'} & ${params['compute']} & $${mean.toFixed(2)} \\pm ${std.toFixed(2)}$\\\\`, n);
};
w.addListener('profile', profileMessageHandler);

const endProfileMessageHandler = function () {
    var spin = document.getElementById("scanningSpinner");
    spin.classList.add("hidden");
};
w.addListener('endprofile', endProfileMessageHandler);

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

function autoWindow(which) {
    var histo = new Int32Array(2048);
    for (var i = 0; i < histo.length; i++) {
        histo[i] = 0;
    }
    var imgResult = imgResultCache[which].data;
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
    document.getElementById(which + "_wc").value = (4096 / histo.length) * wc;
    document.getElementById(which + "_ww").value = (4096 / histo.length) * ww;
    displayAndWindow3DImage(which);
}

function scrollDataset(event) {
    slice = document.getElementById("slice")
    slice.value = slice.valueAsNumber - Math.sign(event.deltaY) * Math.max(1, Math.abs(Math.round(event.deltaY / 100)))
    if ("createEvent" in document) {
        evt = new Event("change", {
            "bubbles": false,
            "canceable": true
        });
        slice.dispatchEvent(evt);
    } else
        slice.fireEvent("onchange");
    event.preventDefault();
}

function scrollResult(event, plane) {
    var keys = Object.getOwnPropertyNames(imgResultCache)
    for (var key in keys) {
        which = keys[key]
        r_slice = document.getElementById(which + "_" + plane + "slice")
        if(r_slice == undefined) {continue;}
        r_slice.value = r_slice.valueAsNumber - Math.sign(event.deltaY) * Math.max(1, Math.abs(Math.round(event.deltaY / 100)))
        if ("createEvent" in document) {
            evt = new Event("change", {
                "bubbles": false,
                "cancelable": true
            });
            r_slice.dispatchEvent(evt);
        } else {
            r_slice.fireEvent("onchange");
        }
    }
    event.preventDefault();
}

var windowing = {
    "cur": false,
    "pre": false,
    "cs": false
};

function startWindowing(which) {
    windowing[which] = true;
}

function endWindowing(which) {
    windowing[which] = false;
}

function windowResult(event, which) {
    if (windowing[which]) {
        ww = document.getElementById(which + "_ww")
        wc = document.getElementById(which + "_wc")

        ww.value = ww.valueAsNumber - event.movementY * 4
        if(ww.value < 0){ww.value=0};
        if(ww.value > 4096){ww.value=4096};
        wc.value = wc.valueAsNumber - event.movementX * 4
        if(wc.value < 0){wc.value=0};
        if(wc.value > 4096){wc.value=4096};

        evt = new Event("change", {
            "bubbles": false,
            "canceable": true
        });
        ww.dispatchEvent(evt);
    }
}

var rotation = {
    "cur": false,
    "pre": false,
    "cs": false
};
var viewPort = {
    "cur": [4.22, -4.27, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
    "pre": [4.22, -4.27, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
    "cs": [4.22, -4.27, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
}

function rot(x, y, z, a, t) {
    var pos = [a[0]*x+a[1]*y+a[2]*z+a[3], 
                       a[4]*x+a[5]*y+a[6]*z+a[7], 
                       a[8]*x+a[9]*y+a[10]*z+a[11],
                       a[12]*x+a[13]*y+a[14]*z+a[15]];
    if(pos[3]!=1) {
        pos[0] = pos[0]/pos[3];
        pos[1] = pos[1]/pos[3];
        pos[2] = pos[2]/pos[3];
        pos[3] = 1;
    }
    pos[0] += t[0];
    pos[1] += t[1];
    pos[2] += t[2];
    return pos;
}

function startRotation(event, which) {
    rotation[which] = {x: event.x, y: event.y};
}

function endRotation(which) {
    rotation[which] = false;
}

function rotateResult(event, which) {
    if (rotation[which]) {
        var dx = (rotation[which].x-event.x) / 100;
        var dy = (rotation[which].y-event.y) / 100;
        /*var dz = 0;
        var h = event.target.clientHeight;
        var w = event.target.clientWidth;
        var y = event.pageY-event.target.offsetTop;
        var x = event.pageX-event.target.offsetLeft;
        var v = [(x-0.5*w)/w, (y-0.5*h)/w, 1];
        v[2] = 1-Math.sqrt(v[1]*v[1]+v[0]*v[0]);
        
        var rz = [v[1], -v[0], 0];
        var ry = [-v[2], 0, v[0]];
        var rx = [0, -v[2], v[1]];
        var l = Math.sqrt(rz[0]*rz[0]+rz[1]*rz[1]);
        rz[0] = rz[0]/l;
        rz[1] = rz[1]/l;
        l = Math.sqrt(ry[0]*ry[0]+ry[2]*ry[2]);
        ry[0] = ry[0]/l;
        ry[2] = ry[2]/l;
        l = Math.sqrt(rx[1]*rx[1]+rx[2]*rx[2]);
        rx[1] = rx[1]/l;
        rx[2] = rx[2]/l;
        */

        rotation[which] = {x: event.x, y: event.y};
        viewPort[which][0] += dy;
        while(viewPort[which][0] < 0) {
            viewPort[which][0] = 2*Math.PI + viewPort[which][0];
        }
        while(viewPort[which][0] >= 2*Math.PI) {
            viewPort[which][0] -= 2*Math.PI;
        }
        viewPort[which][1] -= dx;
        while(viewPort[which][1] < 0) {
            viewPort[which][1] = 2*Math.PI + viewPort[which][1];
        }
        while(viewPort[which][1] >= 2*Math.PI) {
            viewPort[which][1] -= 2*Math.PI;
        }
        //viewPort[which][0] += dy*ry[2];
        /*
        a0 a1 a2       1   0   0     a0  ca*a1+sa*a2  -sa*a1+ca*a2     cb  0   sb     cb*a0+sa*sb*a1-ca*sb*a2  ca*a1+sa*a2  sb*a0-sa*cb*a1+ca*cb*a2       cc  -sc 0   cb*cc*a0+sa*sb*cc*a1-ca*sb*cc*a2+ca*cc*a1+sa*cc*a2   -cb*sc*a0-sa*sb*sc*a1+ca*sb*sc*a2+ca*cc*a1+sa*cc*a2   sb*a0-sa*cb*a1+ca*cb*a2
        a4 a5 a6    *  0   ca  -sa = a4  ca*a5+sa*a6  -sa*a5+ca*a6  *  0   1   0   =  cb*a4+sa*sb*a5-ca*sb*a6  ca*a5+sa*a6  sb*a4-sa*cb*a5+ca*cb*a6    *  sc  cc  0 = cb*cc*a4+sa*sb*cc*a5-ca*sb*cc*a6+ca*cc*a5+sa*cc*a6   -cb*sc*a4-sa*sb*sc*a5+ca*sb*sc*a6+ca*cc*a5+sa*cc*a6   sb*a4-sa*cb*a5+ca*cb*a6
        a8 a9 a10      0   sa  ca    a8  ca*a9+sa*a10 -sa*a9+sa*a10    -sb 0   cb     cb*a8+sa*sb*a9-ca*sb*a10  ca*a9+sa*a10  sb*a8-sa*cb*a9+ca*cb*a10    0   0   1   cb*cc*a8+sa*sb*cc*a9-ca*sb*cc*a10+ca*cc*a9+sa*cc*a10   -cb*sc*a8-sa*sb*sc*a9+ca*sb*sc*a10+ca*cc*a9+sa*cc*a10   sb*a8-sa*cb*a9+ca*cb*a10
        */
        /*var sa = Math.sin(dx);
        var ca = Math.cos(dx);
        var sb = Math.sin(dy);
        var cb = Math.cos(dy);
        var sc = Math.sin(dz);
        var cc = Math.cos(dz);
        //console.log(dx, dy, dz, sa, ca, sb, cb, sc, cc, viewPort[which], which);
        var a = [...viewPort[which]];
        viewPort[which][0] = cb*cc*a[0]+sa*sb*cc*a[1]-ca*sb*cc*a[2]+ca*cc*a[1]+sa*cc*a[2];
        viewPort[which][1] = -cb*sc*a[0]-sa*sb*sc*a[1]+ca*sb*sc*a[2]+ca*cc*a[1]+sa*cc*a[2];
        viewPort[which][2] = sb*a[0]-sa*cb*a[1]+ca*cb*a[2];
        viewPort[which][4] = cb*cc*a[4]+sa*sb*cc*a[5]-ca*sb*cc*a[6]+ca*cc*a[5]+sa*cc*a[6];
        viewPort[which][5] = -cb*sc*a[4]-sa*sb*sc*a[5]+ca*sb*sc*a[6]+ca*cc*a[5]+sa*cc*a[6];
        viewPort[which][6] = sb*a[4]-sa*cb*a[5]+ca*cb*a[6];
        viewPort[which][8] = cb*cc*a[8]+sa*sb*cc*a[9]-ca*sb*cc*a[10]+ca*cc*a[9]+sa*cc*a[10];
        viewPort[which][9] = -cb*sc*a[8]-sa*sb*sc*a[9]+ca*sb*sc*a[10]+ca*cc*a[9]+sa*cc*a[10];
        viewPort[which][10] = sb*a[8]-sa*cb*a[9]+ca*cb*a[10];*/
        /*evt = new Event("change", {
            "bubbles": false,
            "canceable": true
        });
        event.target.dispatchEvent(evt);*/
        displayAndWindow3DImage(which);
    }
}

function isCurrentTabNa(which) {
    if (imgResultCache[which] != undefined) {
        return na_tabs.includes("params-" + imgResultCache[which].params["sequence"] + "-tab");
    }
    return na_tabs.includes(current_tab);
}

function toggleImg() {
    var keys = ["pre", "cur", "cs"];
    var kspace = document.getElementById("kspace_visible");
    for (var key in keys) {
        var which = keys[key];
        var div = document.getElementById(which + "_image");
        var checkbox = document.getElementById(which + "_visible");
        if (imgResultCache[which] == undefined) {
            div.classList.add("hidden");
        } else {
            if (checkbox.checked) {
                div.classList.remove("hidden");
            } else {
                div.classList.add("hidden");
            }
            var div = document.getElementById(which + "_kResult").parentElement.parentElement;
            var div2 = document.getElementById(which + "_3dResult").parentElement.parentElement;
            if (kspace != undefined && kspace.checked) {
                div.classList.remove("hidden");
                div2.classList.add("hidden");
            } else {
                div.classList.add("hidden");
                div2.classList.remove("hidden");
            }
        }
        var keys2 = ["cor", "tra", "sag"];
        for (var key in keys2) {
            var plane = keys2[key];
            var checkbox = document.getElementById(plane + "_visible");
            if(checkbox == undefined) continue;
            var div = document.getElementById(which + "_" + plane + "Result");
            if (checkbox.checked) {
                div.classList.remove("hidden");
            } else {
                div.classList.add("hidden");
            }
            var div = document.getElementById(which + "_" + plane + "slice").parentElement;
            if (checkbox.checked) {
                div.classList.remove("hidden");
            } else {
                div.classList.add("hidden");
            }
        }
    }
}

function toggleAccordion(which) {
    var content = document.getElementById('acc-' + which + '-content');
    var head = document.getElementById('acc-' + which + '-head').children[0];
    if (content.classList.contains("show")) {
        content.classList.remove("show");
        head.classList.add("collapsed");
    } else {
        content.classList.add("show");
        head.classList.remove("collapsed");
    }
}

function displayAndWindow3DImage(which) {
    if (which == undefined) {
        var keys = Object.getOwnPropertyNames(imgResultCache)
        for (var key in keys) {
            displayAndWindow3DImage(keys[key]);
        }
        return;
    }
    var imgResult = imgResultCache[which];
    if (imgResult == undefined) {
        return;
    }

    var xdim = imgResult.xdim;
    var ydim = imgResult.ydim;
    var zdim = imgResult.zdim;

    var size = Math.max(xdim,ydim,zdim);
    var xoff = Math.round((size-xdim)/2);
    var yoff = Math.round((size-ydim)/2);
    var zoff = Math.round((size-zdim)/2);

    var in_slice = document.getElementById(which + "_corslice");
    var cor_slice = parseInt(in_slice.value);
    in_slice = document.getElementById(which + "_sagslice");
    var sag_slice = parseInt(in_slice.value);
    in_slice = document.getElementById(which + "_traslice");
    var tra_slice = parseInt(in_slice.value);

    var show_crosshair = document.getElementById("crosshair_visible").checked;

    var canvas = document.getElementById(which + "_corResult");
    ctx = canvas.getContext('2d');
    canvas.width = size;
    canvas.height = size;
    idata = ctx.createImageData(size, size);

    var cor_result = new Uint8ClampedArray(size * size * 4);
    var tra_result = new Uint8ClampedArray(size * size * 4);
    var sag_result = new Uint8ClampedArray(size * size * 4);
    in_ww = document.getElementById(which + "_ww");
    in_wc = document.getElementById(which + "_wc");
    var ww = parseFloat(in_ww.value) * 0.5
    var wc = parseFloat(in_wc.value)
    for (var x = 0; x < xdim; x++) {
        for (var z = 0; z < zdim; z++) {
            var val = imgResult.data[2*(x + (ydim-cor_slice) * xdim + (zdim - z) * xdim * ydim)] * 4096;
            if (val <= (wc - ww)) {
                val = 0;
            } else if (val >= (wc + ww)) {
                val = 1;
            } else {
                val = (val - (wc - ww)) / (2*ww);
            }
            if (isCurrentTabNa(which)) {
                val = jetmap[Math.round(255*val)];
                if (val == undefined) {
                    val = [0, 0, 0];
                }
            } else {
                val = [val, val, val];
            }
            cor_result[4 * (x+xoff + (z+zoff) * size)] = val[0]*255;
            cor_result[4 * (x+xoff + (z+zoff) * size) + 1] = val[1]*255;
            cor_result[4 * (x+xoff + (z+zoff) * size) + 2] = val[2]*255;
            cor_result[4 * (x+xoff + (z+zoff) * size) + 3] = 255;
        }
    }
    if (isCurrentTabNa(which)) {
        updateCB(which, 2*ww, wc);
    }
    if (show_crosshair) {
        for (var x = 0; x < xdim; x++) {
            if (x > sag_slice - 10 && x < sag_slice + 10) {
                continue;
            }
            cor_result[4 * (x+xoff + (zdim - tra_slice+zoff) * size)] = 255
            cor_result[4 * (x+xoff + (zdim - tra_slice+zoff) * size) + 1] = 0
            cor_result[4 * (x+xoff + (zdim - tra_slice+zoff) * size) + 2] = 0
            cor_result[4 * (x+xoff + (zdim - tra_slice+zoff) * size) + 3] = 255
        }
        for (var z = 0; z < zdim; z++) {
            if (z > (zdim - tra_slice - 10) && z < (zdim - tra_slice + 10)) {
                continue;
            }
            cor_result[4 * (sag_slice+xoff + (z+zoff) * size)] = 0
            cor_result[4 * (sag_slice+xoff + (z+zoff) * size) + 1] = 255
            cor_result[4 * (sag_slice+xoff + (z+zoff) * size) + 2] = 0
            cor_result[4 * (sag_slice+xoff + (z+zoff) * size) + 3] = 255
        }
    }

    idata.data.set(cor_result);
    ctx.putImageData(idata, 0, 0);

    canvas = document.getElementById(which + "_sagResult");
    ctx = canvas.getContext('2d');
    canvas.width = size;
    canvas.height = size;
    idata = ctx.createImageData(size, size);

    in_ww = document.getElementById(which + "_ww");
    in_wc = document.getElementById(which + "_wc");
    var ww = parseFloat(in_ww.value) * 0.5
    var wc = parseFloat(in_wc.value)
    for (var y = 0; y < ydim; y++) {
        for (var z = 0; z < zdim; z++) {
            var val = imgResult.data[2*(sag_slice + y * xdim + (zdim - z) * xdim * ydim)] * 4096;
            if (val <= (wc - ww)) {
                val = 0
            } else if (val >= (wc + ww)) {
                val = 1;
            } else {
                val = (val - (wc - ww)) / (2.0*ww)
            }
            if(isCurrentTabNa(which)){
                val = jetmap[Math.round(255*val)];
                if (val == undefined) {
                    val = [0, 0, 0];
                }
            } else {
                val = [val, val, val];
            }
            sag_result[4 * (y+yoff + (z+zoff) * size)] = val[0]*255;
            sag_result[4 * (y+yoff + (z+zoff) * size) + 1] = val[1]*255;
            sag_result[4 * (y+yoff + (z+zoff) * size) + 2] = val[2]*255;
            sag_result[4 * (y+yoff + (z+zoff) * size) + 3] = 255;
        }
    }
    if (show_crosshair) {
        for (var y = 0; y < ydim; y++) {
            if (y > cor_slice - 10 && y < cor_slice + 10) {
                continue;
            }
            sag_result[4 * (y+yoff + (zdim - tra_slice+zoff) * size)] = 255
            sag_result[4 * (y+yoff + (zdim - tra_slice+zoff) * size) + 1] = 0
            sag_result[4 * (y+yoff + (zdim - tra_slice+zoff) * size) + 2] = 0
            sag_result[4 * (y+yoff + (zdim - tra_slice+zoff) * size) + 3] = 255
        }
        for (var z = 0; z < zdim; z++) {
            if (z > (zdim - tra_slice - 10) && z < (zdim - tra_slice + 10)) {
                continue;
            }
            sag_result[4 * (cor_slice+yoff + (z+zoff) * size)] = 0
            sag_result[4 * (cor_slice+yoff + (z+zoff) * size) + 1] = 0
            sag_result[4 * (cor_slice+yoff + (z+zoff) * size) + 2] = 255
            sag_result[4 * (cor_slice+yoff + (z+zoff) * size) + 3] = 255
        }
    }

    idata.data.set(sag_result);
    ctx.putImageData(idata, 0, 0);

    canvas = document.getElementById(which + "_traResult");
    ctx = canvas.getContext('2d');
    canvas.width = size;
    canvas.height = size;
    idata = ctx.createImageData(size, size);

    in_ww = document.getElementById(which + "_ww");
    in_wc = document.getElementById(which + "_wc");
    var ww = parseFloat(in_ww.value) * 0.5
    var wc = parseFloat(in_wc.value)
    for (var x = 0; x < xdim; x++) {
        for (var y = 0; y < ydim; y++) {
            var val = imgResult.data[2*(x + y * xdim + tra_slice * xdim * ydim)] * 4096;
            if (val <= (wc - ww)) {
                val = 0
            } else if (val >= (wc + ww)) {
                val = 1;
            } else {
                val = (val - (wc - ww)) / (2.0*ww)
            }
            if(isCurrentTabNa(which)){
                val = jetmap[Math.round(255*val)];
                if (val == undefined) {
                    val = [0, 0, 0];
                }
            } else {
                val = [val, val, val];
            }
            tra_result[4 * (x+xoff + (ydim-y+yoff) * size)] = val[0]*255;
            tra_result[4 * (x+xoff + (ydim-y+yoff) * size) + 1] = val[1]*255;
            tra_result[4 * (x+xoff + (ydim-y+yoff) * size) + 2] = val[2]*255;
            tra_result[4 * (x+xoff + (ydim-y+yoff) * size) + 3] = 255;
        }
    }
    if (show_crosshair) {
        for (var x = 0; x < xdim; x++) {
            if (x > sag_slice - 10 && x < sag_slice + 10) {
                continue;
            }
            tra_result[4 * (x+xoff + (cor_slice+yoff) * size)] = 0
            tra_result[4 * (x+xoff + (cor_slice+yoff) * size) + 1] = 0
            tra_result[4 * (x+xoff + (cor_slice+yoff) * size) + 2] = 255
            tra_result[4 * (x+xoff + (cor_slice+yoff) * size) + 3] = 255
        }
        for (var y = 0; y < ydim; y++) {
            if (y > cor_slice - 10 && y < cor_slice + 10) {
                continue;
            }
            tra_result[4 * (sag_slice+xoff + (y+yoff) * size)] = 0
            tra_result[4 * (sag_slice+xoff + (y+yoff) * size) + 1] = 255
            tra_result[4 * (sag_slice+xoff + (y+yoff) * size) + 2] = 0
            tra_result[4 * (sag_slice+xoff + (y+yoff) * size) + 3] = 255
        }
    }

    idata.data.set(tra_result);
    ctx.putImageData(idata, 0, 0);

    canvas = document.getElementById(which + "_3dResult");
    ctx = canvas.getContext('2d');

    var result = new Uint8ClampedArray(size * size * 4);
    var zindex = new Array(size*size);
    canvas.width = size;
    canvas.height = size;
    idata = ctx.createImageData(size, size);


    var angles = viewPort[which];
    var sa = Math.sin(angles[0]);
    var ca = Math.cos(angles[0]);
    var sb = Math.sin(angles[1]);
    var cb = Math.cos(angles[1]);
    var sc = Math.sin(angles[2]);
    var cc = Math.cos(angles[2]);
    //console.log(dx, dy, dz, sa, ca, sb, cb, sc, cc, viewPort[which], which);
    var a = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
    /*a[0] = ca*cb;
    a[1] = ca*sb*sc-sa*cc;
    a[2] = ca*sb*cc+sa*sc;
    a[4] = sa*cb;
    a[5] = sa*sb*sc+ca*cc;
    a[6] = sa*sb*cc-ca*sc;
    a[8] = -sb;
    a[9] = cb*sc;
    a[10] = cb*cc;*/
    a[0] = cb;
    a[1] = -cc*sb;
    a[2] = sb*sc;
    a[4] = ca*sb;
    a[5] = ca*cb*cc-sa*sc;
    a[6] = -cc*sa-ca*cb*sc;
    a[8] = sa*sb;
    a[9] = ca*sc+cb*cc*sa;
    a[10] = ca*cc-cb*sa*sc;

    var t = [xoff, 0, 0];

    for(var x=0;x<xdim;x++) {
        for(var y=0;y<ydim;y++) {
            var z = zdim-tra_slice;
            var pos = rot(x-0.5*xdim, y-0.5*ydim, z-0.5*zdim, a, t);
            var nx = Math.round(pos[0]+0.5*xdim);
            var ny = Math.round(pos[1]+0.5*ydim);
            if(nx<0 || nx>=size || ny<0 || ny>=size) {continue;}
            if (zindex[nx+ny*size]==undefined || pos[2] > zindex[nx+ny*size]) {
                zindex[nx+ny*size] = pos[2];
                if(x==0 || y==0 || x==xdim-1 || y==ydim-1) {
                    result[4*(nx+ny*size)] = 255;
                    result[4*(nx+ny*size)+1] = 0;
                    result[4*(nx+ny*size)+2] = 0;
                    result[4*(nx+ny*size)+3] = 255;
                } else {
                    result[4*(nx+ny*size)] = tra_result[4*(x+xoff+(y+yoff)*size)];
                    result[4*(nx+ny*size)+1] = tra_result[4*(x+xoff+(y+yoff)*size)+1];
                    result[4*(nx+ny*size)+2] = tra_result[4*(x+xoff+(y+yoff)*size)+2];
                    result[4*(nx+ny*size)+3] = tra_result[4*(x+xoff+(y+yoff)*size)+3];
                }
            }
        }
    }
    for(var x=0;x<xdim;x++) {
        for(var z=0;z<zdim;z++) {
            var y = cor_slice;
            var pos = rot(x-0.5*xdim, y-0.5*ydim, z-0.5*zdim, a, t);
            var nx = Math.round(pos[0]+0.5*xdim);
            var ny = Math.round(pos[1]+0.5*ydim);
            if(nx<0 || nx>=size || ny<0 || ny>=size) {continue;}
            if (zindex[nx+ny*size]==undefined || pos[2] > zindex[nx+ny*size]) {
                zindex[nx+ny*size] = pos[2];
                if(x==0 || z==0 || x==xdim-1 || z==zdim-1) {
                    result[4*(nx+ny*size)] = 0;
                    result[4*(nx+ny*size)+1] = 0;
                    result[4*(nx+ny*size)+2] = 255;
                    result[4*(nx+ny*size)+3] = 255;
                } else {
                    result[4*(nx+ny*size)] = cor_result[4*(x+xoff+(z+zoff)*size)];
                    result[4*(nx+ny*size)+1] = cor_result[4*(x+xoff+(z+zoff)*size)+1];
                    result[4*(nx+ny*size)+2] = cor_result[4*(x+xoff+(z+zoff)*size)+2];
                    result[4*(nx+ny*size)+3] = cor_result[4*(x+xoff+(z+zoff)*size)+3];
                }
            }
        }
    }
    for(var z=0;z<zdim;z++) {
        for(var y=0;y<ydim;y++) {
            var x = sag_slice;
            var pos = rot(x-0.5*xdim, y-0.5*ydim, z-0.5*zdim, a, t);
            var nx = Math.round(pos[0]+0.5*xdim);
            var ny = Math.round(pos[1]+0.5*ydim);
            if(nx<0 || nx>=size || ny<0 || ny>=size) {continue;}
            if (zindex[nx+ny*size]==undefined || pos[2] > zindex[nx+ny*size]) {
                zindex[nx+ny*size] = pos[2];
                if(z==0 || y==0 || z==zdim-1 || y==ydim-1) {
                    result[4*(nx+ny*size)] = 0;
                    result[4*(nx+ny*size)+1] = 255;
                    result[4*(nx+ny*size)+2] = 0;
                    result[4*(nx+ny*size)+3] = 255;
                }
                else {
                    result[4*(nx+ny*size)] = sag_result[4*(y+yoff+(z+zoff)*size)];
                    result[4*(nx+ny*size)+1] = sag_result[4*(y+yoff+(z+zoff)*size)+1];
                    result[4*(nx+ny*size)+2] = sag_result[4*(y+yoff+(z+zoff)*size)+2];
                    result[4*(nx+ny*size)+3] = sag_result[4*(y+yoff+(z+zoff)*size)+3];
                }
            }
        }
    }

    idata.data.set(result);
    ctx.putImageData(idata, 0, 0);


    k_canvas = document.getElementById(which + "_kResult");
    k_ctx = k_canvas.getContext('2d');
    size = size<128?128:size;
    k_canvas.width = size;
    k_canvas.height = size;
    kdata = k_ctx.createImageData(size, size);

    var mult = document.getElementById(which + "_kspacemult").valueAsNumber;
    k_result = new Uint8ClampedArray(size * size * 4);

    xoff = Math.round((size-xdim)/2);
    yoff = Math.round((size-ydim)/2);
    for (var x = 0; x < xdim; x++) {
        for (var y = 0; y < ydim; y++) {
            val_re = imgResult.kSpace[2*(x + y*xdim + tra_slice * xdim * ydim)];
            val_im = imgResult.kSpace[2*(x + y*xdim + tra_slice * xdim * ydim)+1];
            val = mult*(Math.sqrt(val_re*val_re+val_im*val_im));
            val = mult*Math.abs(val_re);
            if (val < 0) val = 0;
            if (val > 255) val = 255;
            k_result[4 * (x+xoff+(y+yoff)*size)] = val
            k_result[4 * (x+xoff+(y+yoff)*size) + 1] = val
            k_result[4 * (x+xoff+(y+yoff)*size) + 2] = val
            k_result[4 * (x+xoff+(y+yoff)*size) + 3] = 255
        }
    }
    kdata.data.set(k_result);
    k_ctx.putImageData(kdata, 0, 0);
}

function downloadKSpace(which) {
    downloadArray(imgResultCache[which].raw_data, "kspace_" + which + ".csv");
}

function downloadImage(which) {
    downloadArray(imgResultCache[which].data, "image_" + which + ".csv");
}

function downloadArray(array, name) {
    const csvContent = array.join(",");
    const blob = new Blob([csvContent], {
        type: 'text/csv;charset=utf-8;'
    });
    const url = URL.createObjectURL(blob);

    var link = document.getElementById("downloadLink");
    link.setAttribute("href", url);
    link.setAttribute("download", name);
    link.click()

    URL.revokeObjectURL(link.href)
}

function fillDatasets() {
    sel = document.getElementById("datasetPath");
    for (var p in datasets) {
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
        let r = color[0] * 256;
        let g = color[1] * 256;
        let b = color[2] * 256;
        ctx.fillStyle = 'rgb(' + r + ',' + g + ',' + b + ')';
        //ctx.fillRect(0, x * canvas.height / 256, canvas.width, canvas.height / 256);
        ctx.fillRect(Math.floor(x * (canvas.width / 256)), 0, Math.ceil(canvas.width / 256), canvas.height);
    }
}

function setTabs() {
    var tabId = "params-" + selectedSequence;
    var tabHeadId = "params-" + selectedSequence + "-tab";
    selected_tab = tabHeadId;
    elems = document.getElementById("sequence").getElementsByClassName("nav-link")
    for (var x = 0; x < elems.length; x++) {
        if (elems[x].id == tabHeadId) {
            elems[x].classList.add("active");
        } else {
            elems[x].classList.remove("active");
        }
    }
    elems = document.getElementById("sequence").getElementsByClassName("tab-pane")
    for (var x = 0; x < elems.length; x++) {
        if (elems[x].id == tabId) {
            elems[x].classList.add("active", "show");
        } else {
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
    document.getElementById("kspace").setAttribute("disabled", "disabled");
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
        var ydim = Math.round(params["ydim"]);
        var zdim = Math.round(params["zdim"]);
        var params = [xdim, ydim, zdim, xlines, ylines, fmin, fmax, noIfft];
        w.sendQuery("reco", params);
    } else {
        document.getElementById("kspace").removeAttribute("disabled");
        displayAndWindow3DImage();
    }
}

function read_params(param_div, params) {
    for (var child_id in param_div.children) {
        var child = param_div.children[child_id];
        if (child.children == undefined) {
            continue;
        }
        for (var input_id in child.children) {
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

function check_size_fft(x) {
    if(x>=128) return x;
    if(x>64) return 128;
    if(x>32) return 64;
    if(x>16) return 32;
    if(x>8) return 16;
    if(x>4) return 8;
    if(x>2) return 4;
    return 2;
}

function startScan() {
    resultTimer = performance.now();

    current_tab = selected_tab;
    r = document.getElementById("result");
    r.classList.add("hidden");


    var params = {
        sequence: selectedSequence
    };

    var param_div = document.getElementById("params-general");
    params = read_params(param_div, params);
    params["xdim"] = Math.ceil(params["xdim"] / 2) * 2;
    params["ydim"] = Math.ceil(params["ydim"] / 2) * 2;
    params["zdim"] = Math.ceil(params["zdim"] / 2) * 2;
    param_div = document.getElementById("params-noise");
    params = read_params(param_div, params);
    param_div = document.getElementById("params-kspace");
    params = read_params(param_div, params);
    param_div = document.getElementById("params-cs");
    params = read_params(param_div, params);

    param_div = document.getElementById("params-" + selectedSequence);
    params = read_params(param_div, params);

    if (selectedSequence == "TQSQR") {
        var param_div = document.getElementById("params-SQ");
        var tq_params = {};
        tq_params = read_params(param_div, tq_params);
        var param_div = document.getElementById("params-TQ");
        var sq_params = {};
        sq_params = read_params(param_div, sq_params);
        params["tq_params"] = tq_params;
        params["sq_params"] = sq_params;
    }

    var spin = document.getElementById("scanningSpinner");
    spin.classList.remove("hidden");

    if (params["compute"] == "WebASM") {
        w.sendQuery("simulateImageFast", params);
    } else {
        w.sendQuery("simulateImage", params);
    }
}

function profile() {
    current_tab = selected_tab;
    r = document.getElementById("result");
    r.classList.add("hidden");

    var params = {
        sequence: "SE"
    };

    var param_div = document.getElementById("params-general");
    params = read_params(param_div, params);
    params["xdim"] = Math.ceil(params["xdim"] / 2) * 2;
    params["ydim"] = Math.ceil(params["ydim"] / 2) * 2;
    params["zdim"] = Math.ceil(params["zdim"] / 2) * 2;
    param_div = document.getElementById("params-noise");
    params = read_params(param_div, params);
    param_div = document.getElementById("params-kspace");
    params = read_params(param_div, params);
    param_div = document.getElementById("params-cs");
    params = read_params(param_div, params);

    param_div = document.getElementById("params-SE");
    params = read_params(param_div, params);

    var spin = document.getElementById("scanningSpinner");
    spin.classList.remove("hidden");

    params_list = [10];
    params["count"] = 1;
    params["xdim"] = "256";
    params["ydim"] = "256";
    params["zdim"] = "256";
    params["nearest"] = "0";
    params["fft"] = "3D";
    params["sequence"] = "SE";
    params["img_noise_type"] = "0";
    params["cs"] = "0";

    params["compute"] = "JS";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["fft"] = "2D";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["fft"] = "3D";
    params["nearest"] = "2";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["fft"] = "2D";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["fft"] = "3D";
    params["nearest"] = "0";

    params["zdim"] = "64";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["fft"] = "2D";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["fft"] = "3D";
    params["nearest"] = "2";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["fft"] = "2D";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";

    params["zdim"] = "256";
    params["img_noise_type"] = "2";
    params["nearest"] = "0";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["nearest"] = "2";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["zdim"] = "64";
    params["nearest"] = "0";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";
    params["nearest"] = "2"
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";

    params["cs"] = "2";
    params["nearest"] = "0";
    params_list.push({
        ...params
    });
    params["compute"] = "WebASM";
    params_list.push({
        ...params
    });
    params["compute"] = "JS";

    w.sendQuery("profile", params_list);
}

function displayDataSet() {
    //display3DImage(document.getElementById("imgPD"), array_pd);
    //display3DImage(document.getElementById("imgT1"), array_t1);
    //display3DImage(document.getElementById("imgT2"), array_t2);
    //display3DImage(document.getElementById("imgT2s"), array_t2s);
}

function updateScale() {
    var scale = parseFloat(document.getElementById("scale").value) * 0.01;
    var xdim = document.getElementById("xdim");
    xdim.value = Math.round(scale * parseFloat(xdim.max));
    var ydim = document.getElementById("ydim");
    ydim.value = Math.round(scale * parseFloat(ydim.max));
    var zdim = document.getElementById("zdim");
    zdim.value = Math.round(scale * parseFloat(zdim.max));
    updateRes();
}

function updateRes() {
    var xdim = document.getElementById("xdim");
    var ydim = document.getElementById("ydim");
    var zdim = document.getElementById("zdim");

    scale = xdim.value/362.0;
    var xsize = document.getElementById("xsize");
    xsize.value = (0.5 / scale).toFixed(2);
    scale = ydim.value/434.0;
    var ysize = document.getElementById("ysize");
    ysize.value = (0.5 / scale).toFixed(2);
    scale = zdim.value/362.0;
    var zsize = document.getElementById("zsize");
    zsize.value = (0.5 / scale).toFixed(2);
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
    var ydim = parseInt(document.getElementById("ydim").value)
    var zdim = parseInt(document.getElementById("zdim").value)
    var time = document.getElementById("ir_time");

    if (ti + te >= tr) {
        time.innerText = "TE+TI has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr * ydim * zdim);
    }
}

function updateSETime() {
    var te = parseFloat(document.getElementById("se_te").value)
    var tr = parseFloat(document.getElementById("se_tr").value)
    var ydim = parseInt(document.getElementById("ydim").value)
    var zdim = parseInt(document.getElementById("zdim").value)
    var time = document.getElementById("se_time");

    if (te >= tr) {
        time.innerText = "TE has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr * ydim * zdim);
    }
}

function updatePSIFTime() {
    var te = parseFloat(document.getElementById("psif_te").value)
    var tr = parseFloat(document.getElementById("psif_tr").value)
    var fa = parseFloat(document.getElementById("psif_fa").value)
    var ydim = parseInt(document.getElementById("ydim").value)
    var zdim = parseInt(document.getElementById("zdim").value)
    var time = document.getElementById("psif_time");

    if (fa > 180) {
        document.getElementById("psif_fa").value = 180;
    }
    if (fa < -180) {
        document.getElementById("psif_fa").value = -180;
    }

    time.innerText = formatTime(tr * ydim * zdim);
}

function updateFISPTime() {
    var te = parseFloat(document.getElementById("fisp_te").value)
    var tr = parseFloat(document.getElementById("fisp_tr").value)
    var fa = parseFloat(document.getElementById("fisp_fa").value)
    var ydim = parseInt(document.getElementById("ydim").value)
    var zdim = parseInt(document.getElementById("zdim").value)
    var time = document.getElementById("fisp_time");

    if (fa > 180) {
        document.getElementById("fisp_fa").value = 180;
    }
    if (fa < -180) {
        document.getElementById("fisp_fa").value = -180;
    }

    time.innerText = formatTime(te * 2 * ydim * zdim);
}

function updateSGRETime() {
    var te = parseFloat(document.getElementById("sgre_te").value)
    var tr = parseFloat(document.getElementById("sgre_tr").value)
    var fa = parseFloat(document.getElementById("sgre_fa").value)
    var ydim = parseInt(document.getElementById("ydim").value)
    var zdim = parseInt(document.getElementById("zdim").value)
    var time = document.getElementById("sgre_time");

    if (fa > 180) {
        document.getElementById("sgre_fa").value = 180;
    }
    if (fa < -180) {
        document.getElementById("sgre_fa").value = -180;
    }

    if (te >= tr) {
        time.innerText = "TE has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr * ydim * zdim);
    }
}

function updateBalancedSSFPTime() {
    var te = parseFloat(document.getElementById("bssfp_te").value)
    var fa = parseFloat(document.getElementById("bssfp_fa").value)
    var ydim = parseInt(document.getElementById("ydim").value)
    var zdim = parseInt(document.getElementById("zdim").value)
    var time = document.getElementById("bssfp_time");

    if (fa > 180) {
        document.getElementById("bssfp_fa").value = 180;
    }
    if (fa < -180) {
        document.getElementById("bssfp_fa").value = -180;
    }

    document.getElementById("bssfp_tr").value = te * 2;

    time.innerText = formatTime(te * 2 * ydim * zdim);
}

function updateNaTime() {
    var te = parseFloat(document.getElementById("na_te").value)
    var tr = parseFloat(document.getElementById("na_tr").value)
    var ydim = parseInt(document.getElementById("ydim").value)
    var zdim = parseInt(document.getElementById("zdim").value)
    var time = document.getElementById("na_time");

    if (te >= tr) {
        time.innerText = "TE has to be smaller than TR";
    } else {
        time.innerText = formatTime(tr * ydim * zdim);
    }
}

function updateSQTime() {
    var te_start = parseFloat(document.getElementById("sq_te_start").value)
    var te_end = parseFloat(document.getElementById("sq_te_end").value)
    var te_step = parseFloat(document.getElementById("sq_te_step").value)
    var ydim = parseInt(document.getElementById("ydim").value)
    var zdim = parseInt(document.getElementById("zdim").value)
    var te = document.getElementById("sq_te");

    var tes = "" + te_start;
    for (var t = te_start + te_step; t <= te_end; t += te_step) {
        tes += ", " + t;
    }
    te.innerText = tes;
}

function updateTQTime() {
    var te_start = parseFloat(document.getElementById("tq_te_start").value)
    var te_end = parseFloat(document.getElementById("tq_te_end").value)
    var te_step = parseFloat(document.getElementById("tq_te_step").value)
    var ydim = parseInt(document.getElementById("ydim").value)
    var zdim = parseInt(document.getElementById("zdim").value)
    var te = document.getElementById("tq_te");

    var tes = "" + te_start;
    for (var t = te_start + te_step; t <= te_end; t += te_step) {
        tes += ", " + t;
    }
    te.innerText = tes;
}

function updateTQSQRTime() {}

function updateCB(which, ww, wc) {
    var cb = document.getElementById(which+"_colorBar");
    if (imgResultCache[which] != undefined && (imgResultCache[which].params["sequence"] == "SQ" || imgResultCache[which].params["sequence"] == "Na")) {
        cb.previousElementSibling.textContent = Math.round((wc-0.5*ww)*140/4096) + "mmol";
        cb.nextElementSibling.textContent = Math.round((wc+0.5*ww)*140/4096) + "mmol";
    } else {
        cb.previousElementSibling.textContent = ((wc-0.5*ww)/4096) + "aU";
        cb.nextElementSibling.textContent = ((wc+0.5*ww)/4096) + "aU";
    }
}

function formatTime(time) {
    d = new Date(time);
    timeString = "";
    if (d.getFullYear() > 1970) {
        timeString += (d.getFullYear() - 1970) + " years ";
    }
    if (d.getDate() > 1) {
        timeString += (d.getDate() - 1) + " days ";
    }
    timeString += (d.getHours() < 11 ? "0" : "") + (d.getHours() - 1) + ":" + (d.getMinutes() < 10 ? "0" : "") + d.getMinutes() + ":" + (d.getSeconds() < 10 ? "0" : "") + d.getSeconds();
    return timeString;
}