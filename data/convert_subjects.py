import numpy as np
import scipy.ndimage

# Name, Label, T1, T2, T2*, PD, filename
params = {
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

xdim = 362
ydim = 434
zdim = 362
shape = (zdim*ydim*xdim)

t1 = np.zeros(shape, dtype=np.float32)
t2 = np.zeros(shape, dtype=np.float32)
pd = np.zeros(shape, dtype=np.float32)
s = np.zeros(shape, dtype=np.float32)

for p in params:
    with open(params[p][-1], "rb") as f:
        a = np.array(np.frombuffer(f.read(), dtype=np.int8), dtype=np.float32)
        a = (a+128) / 256
        t1 = t1 + a*params[p][2]
        t2 = t2 + a*params[p][3]
        pd = pd + a*params[p][5]
        s = s + a

t1 = t1/s
t2 = t2/s
pd = pd/s
zoom = (2**8/xdim, 2**8/ydim, 2**8/zdim)
t1 = scipy.ndimage.zoom(t1.reshape((xdim,ydim,zdim)), zoom, order=0).flatten()
pd = scipy.ndimage.zoom(pd.reshape((xdim,ydim,zdim)), zoom, order=0).flatten()
t2 = scipy.ndimage.zoom(t2.reshape((xdim,ydim,zdim)), zoom, order=0).flatten()
with open("t1.bin", "wb") as f:
    mi = np.min(t1)
    ma = np.max(t1)
    ex = np.array(256 * (t1-mi) / (ma-mi), dtype=np.uint8)
    f.write(np.array([mi,ma], dtype=np.float32))
    f.write(ex.tobytes())

with open("t2.bin", "wb") as f:
    mi = np.min(t2)
    ma = np.max(t2)
    ex = np.array(256 * (t2-mi) / (ma-mi), dtype=np.uint8)
    f.write(np.array([mi,ma], dtype=np.float32))
    f.write(ex.tobytes())

with open("pd.bin", "wb") as f:
    mi = np.min(pd)
    ma = np.max(pd)
    ex = np.array(256 * (pd-mi) / (ma-mi), dtype=np.uint8)
    f.write(np.array([mi,ma], dtype=np.float32))
    f.write(ex.tobytes())