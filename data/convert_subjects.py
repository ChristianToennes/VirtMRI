import numpy as np
import scipy.ndimage
import gzip
import os
import os.path
import netCDF4
import nibabel

# Name, Label, T1, T2, T2*, PD, filename
params_15 = {
    0: ["Background", 0, 0, 0, 0, 0, "subject{}_bck_v{}"],
    1: ["CSF", 1, 2569.0, 329, 58, 1, "subject{}_csf_v{}"],
    2: ["Grey Matter", 2.0, 833, 83, 69, 0.86, "subject{}_gm_v{}"], # 577
    3: ["White Matter", 3.0, 500, 70, 61, 0.77, "subject{}_wm_v{}"], # 346
    4: ["Fat", 4, 350.0, 70.0, 58, 1, "subject{}_fat_v{}"],
    5: ["Muscle", 5, 900.0, 47, 30, 1, "subject{}_muscles_v{}"],
    6: ["Muscle / Skin", 6, 569.0, 329, 58, 1, "subject{}_muscles_skin_v{}"],
    7: ["Skull", 7, 0, 0, 0, 0, "subject{}_skull_v{}"],
    8: ["Vessels", 8, 0, 0, 0, 0, "subject{}_vessels_v{}"],
    9: ["Around fat", 9, 500.0, 70, 61, 0.77, "subject{}_fat2_v{}"],
    10: ["Dura Matter", 10, 2569.0, 329, 58, 1, "subject{}_dura_v{}"],
    11: ["Bone Marrow", 11, 500.0, 70, 61, 0.77, "subject{}_marrow_v{}"]
}

params_3 = { # 3T params
    0: ["Background", 0, 0, 0, 0, 0, "subject{}_bck_v{}"],
    1: ["CSF", 1, 2000.0, 329, 58, 1, "subject{}_csf_v{}"],
    2: ["Grey Matter", 2.0, 1433.2, 92.6, 69, 0.86, "subject{}_gm_v{}"], # 993
    3: ["White Matter", 3.0, 866.9, 60.8, 61, 0.77, "subject{}_wm_v{}"], # 600
    4: ["Fat", 4, 346, 70.0, 58, 1, "subject{}_fat_v{}"],
    5: ["Muscle", 5, 1232.9, 37.2, 30, 1, "subject{}_muscles_v{}"],
    6: ["Muscle / Skin", 6, 377.0, 97.5, 58, 1, "subject{}_muscles_skin_v{}"],
    7: ["Skull", 7, 0, 0, 0, 0, "subject{}_skull_v{}"],
    8: ["Vessels", 8, 1984.4, 275.0, 0, 0, "subject{}_vessels_v{}"],
    9: ["Around fat", 9, 346, 70, 58, 0.77, "subject{}_fat2_v{}"],
    10: ["Dura Matter", 10, 2569.0, 329, 58, 1, "subject{}_dura_v{}"],
    11: ["Bone Marrow", 11, 365.0, 133.0, 61, 0.77, "subject{}_marrow_v{}"]
}

nii_names = {0: "BCK", 1:"CSF", 2:"GM", 3:"WM", 4:"FAT", 5:"MUSCLES", 6:"SKIN-MUSCLES", 7:"SKULL", 8:"VESSELS", 9:"FAT2", 10:"DURA", 11:"MARROW"}
brainWeb_names = {0: "subject{}_bck_v", 1:"subject{}_csf_v", 2:"subject{}_gm_v", 3:"subject{}_wm_v", 4:"subject{}_fat_v", 5:"subject{}_muscles_v",
 6:"subject{}_muscles_skin_v", 7:"subject{}_skull_v", 8:"subject{}_vessels", 9:"subject{}_fat2_v", 10:"subject{}_dura_v", 11:"subject{}_marrow_v"}

xdim = 362
ydim = 434
zdim = 362
shape = (zdim*ydim*xdim)

def read_bin(params, in_dir, sub):
    print(in_dir, sub)
    t1 = np.zeros(shape, dtype=np.float32)
    t2 = np.zeros(shape, dtype=np.float32)
    pd = np.zeros(shape, dtype=np.float32)
    s = np.zeros(shape, dtype=np.float32)

    for p in params:
        with open(os.path.join(in_dir, params[p][-1].format(*sub)), "rb") as f:
            a = np.array(np.frombuffer(f.read(), dtype=np.int8), dtype=np.float32)
            a = (a+128) / 256
            t1 = t1 + a*params[p][2]
            t2 = t2 + a*params[p][3]
            pd = pd + a*params[p][-2]
            s = s + a

    print(np.min(s), np.max(s), np.mean(s), np.median(s))
    #t1 = t1/s
    #t2 = t2/s
    #pd = pd/s
    zoom = (2**8/xdim, 2**8/ydim, 2**8/zdim)
    t1 = scipy.ndimage.zoom(t1.reshape((xdim,ydim,zdim)), zoom, order=0).flatten()
    pd = scipy.ndimage.zoom(pd.reshape((xdim,ydim,zdim)), zoom, order=0).flatten()
    t2 = scipy.ndimage.zoom(t2.reshape((xdim,ydim,zdim)), zoom, order=0).flatten()

    return t1, t2, pd

def read_minc(params, names, in_dir, sub=(), trans=None, nib=False):
    print(in_dir)

    t1 = np.zeros(shape, dtype=np.float32)
    t2 = np.zeros(shape, dtype=np.float32)
    pd = np.zeros(shape, dtype=np.float32)
    s = np.zeros(shape, dtype=np.float32)

    for p in params:
        print(p, os.path.join(in_dir, names[p].format(*sub)+".mnc"))
        if nib:
            d = nibabel.load(os.path.join(in_dir, names[p].format(*sub)+".mnc"))
            img = d.get_fdata()
            a = np.array(img, dtype=np.float32).flatten()
        else:
            with gzip.open(os.path.join(in_dir, names[p].format(*sub)+".mnc.gz"), "rb") as f:
                d = netCDF4.Dataset("in_memory.mnc", memory=np.array(f.read()))
                img = np.array(d.variables["image"])
                a = np.array(img, dtype=np.float32).flatten()
        if trans is not None:
            print(np.min(a), np.max(a), np.mean(a), np.median(a))
            a = trans(a)
            print(np.min(a), np.max(a), np.mean(a), np.median(a))
        t1 = t1 + a*params[p][2]
        t2 = t2 + a*params[p][3]
        pd = pd + a*params[p][-2]
        s = s + a

    print(np.min(s), np.max(s), np.mean(s), np.median(s), np.count_nonzero(s>1.1), np.count_nonzero(s<0.9))
    s[s==0] = 1
    t1 = t1/s
    t2 = t2/s
    pd = pd/s
    
    zoom = (2**8/xdim, 2**8/ydim, 2**8/zdim)
    t1 = scipy.ndimage.zoom(t1.reshape((xdim,ydim,zdim)), zoom, order=0).flatten()
    pd = scipy.ndimage.zoom(pd.reshape((xdim,ydim,zdim)), zoom, order=0).flatten()
    t2 = scipy.ndimage.zoom(t2.reshape((xdim,ydim,zdim)), zoom, order=0).flatten()

    return t1, t2, pd

def write_files(t1, t2, pd, sub, out_dir):
    if not os.path.exists(os.path.join(out_dir, sub[0])):
        os.makedirs(os.path.join(out_dir, sub[0]))
    with gzip.open(os.path.join(out_dir, sub[0], "t1.bin.gz"), "wb") as f:
        mi = np.min(t1)
        ma = np.max(t1)
        ex = np.array(255 * (t1-mi) / (ma-mi), dtype=np.uint8)
        print("t1", mi, ma, np.mean(ex), np.max(ex), np.max(255 * (t1-mi) / (ma-mi)))
        er = (255 * (t1-mi) / (ma-mi)) - ex
        print("er", np.min(er), np.max(er), np.mean(er), np.median(er))
        f.write(np.array([mi,ma], dtype=np.float32).tobytes())
        f.write(ex.tobytes())

    with gzip.open(os.path.join(out_dir, sub[0], "t2.bin.gz"), "wb") as f:
        mi = np.min(t2)
        ma = np.max(t2)
        ex = np.array(255 * (t2-mi) / (ma-mi), dtype=np.uint8)
        print("t2", mi, ma, np.mean(ex), np.max(ex), np.max(255 * (t2-mi) / (ma-mi)))
        er = (255 * (t2-mi) / (ma-mi)) - ex
        print("er", np.min(er), np.max(er), np.mean(er), np.median(er))
        f.write(np.array([mi,ma], dtype=np.float32).tobytes())
        f.write(ex.tobytes())

    with gzip.open(os.path.join(out_dir, sub[0], "pd.bin.gz"), "wb") as f:
        mi = np.min(pd)
        ma = np.max(pd)
        ex = np.array(255 * (pd-mi) / (ma-mi), dtype=np.uint8)
        print("pd", mi, ma, np.mean(ex), np.max(ex), np.max(255 * (pd-mi) / (ma-mi)))
        er = (255 * (pd-mi) / (ma-mi)) - ex
        print("er", np.min(er), np.max(er), np.mean(er), np.median(er))
        f.write(np.array([mi,ma], dtype=np.float32).tobytes())
        f.write(ex.tobytes())

#t1,t2,pd = read_bin(params_15, "", ("04",".bin"))
#write_files(t1, t2, pd, ("04",), "1t")
#t1,t2,pd = read_bin(params_15, "", ("05",".rawb"))
#write_files(t1, t2, pd, ("05",), "1t")
#t1,t2,pd = read_bin(params_3, "", ("04",".bin"))
#write_files(t1, t2, pd, ("04",), "3t")
#t1,t2,pd = read_bin(params_3, "", ("05",".rawb"))
#write_files(t1, t2, pd, ("05",), "3t")

#t1,t2,pd = read_minc(params_15, nii_names, "mni_colin27_2008_fuzzy_minc2", nib=True)
#write_files(t1, t2, pd, ("bw",), "1t")
#t1,t2,pd = read_minc(params_3, nii_names, "mni_colin27_2008_fuzzy_minc2", nib=True)
#write_files(t1, t2, pd, ("bw",), "3t")

t1,t2,pd = read_minc(params_15, brainWeb_names, "", ("05",), trans=lambda a: (a+128)/255)
write_files(t1, t2, pd, ("05",), "1t")
t1,t2,pd = read_minc(params_3, brainWeb_names, "", ("05",), trans=lambda a: (a+128)/255)
write_files(t1, t2, pd, ("05",), "3t")

t1,t2,pd = read_minc(params_15, brainWeb_names, "", ("54",), trans=lambda a: (a+128)/255)
write_files(t1, t2, pd, ("54",), "1t")
t1,t2,pd = read_minc(params_3, brainWeb_names, "", ("54",), trans=lambda a: (a+128)/255)
write_files(t1, t2, pd, ("54",), "3t")