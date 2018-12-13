from symsplash.collapse import SSHalo, NFWHalo, symmetron
import numpy as np
from matplotlib import pyplot as plt


dir = "output/SS_smooth_z/"
halo_class = "SS"
cache_dir = "cache_smooth/"

f_u = np.linspace(0.1, 1, 4)
z_ssb_u = np.linspace(1, 3, 20)
size_u = np.linspace(10, 50, 5)

f, z_ssb, size = np.meshgrid(f_u, z_ssb_u, size_u)
f, z_ssb, size = f.flatten(), z_ssb.flatten(), size.flatten()
beta = np.sqrt(size**2.*f/(1+z_ssb)**3.)

res = np.zeros(len(z_ssb))

for i in xrange(len(z_ssb)):
    print z_ssb[i], beta[i], size[i]

    if(halo_class=="SS"):
        halo = SSHalo(z_ssb=z_ssb[i], beta=beta[i], size=size[i], load_dir=cache_dir)

    sym = symmetron(halo, ts=np.linspace(0.19, 1, 20))
    pos2 = halo.pos(np.linspace(0.2, 1, 100), verbose=False, symmetron=sym)

    x = pos2[:, 0]
    x_p = pos2[:, 1]
    res[i] = x[x_p>0].max()

    np.save(dir+"/raw/x"+str(i)+".npy", x)
    np.save(dir+"/raw/x_p"+str(i)+".npy", x_p)

    y = np.linspace(1e-2, 2, 1e2)
    chichipout = np.zeros((20, len(y)))
    for j in xrange(20):
        chichipout[j, :] = sym.chichip[j](y)
    np.save(dir+"/raw/chichip"+str(i)+".npy", chichipout)

np.save(dir+'/rsp.npy', res)
np.save(dir+'/beta.npy', beta)
np.save(dir+'/z_ssb.npy', z_ssb)
np.save(dir+'/size.npy', size)
