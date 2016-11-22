from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import types

def pressure_in_abaqus(coord_list_pressure_data):
    x = []
    y = []
    z = []
    with open(coord_list_pressure_data) as coords_file:
        for line in coords_file:
            line = line.strip().split(';')
            if line != []:
                x.append(float(line[0]))
                y.append(float(line[1]))
                z.append(float(line[3]))
            else:
                continue

    #with open(pore_pressure_data) as pore_pr:
    #    for line in pore_pr:
    #        line = line.strip().split()
    #        if line != []:
    #            z.append(float(line[1]))
    #        else:
    #            continue

    xi = np.linspace(min(x),max(x), 500)
    yi = np.linspace(min(y), max(y), 500)
    xig, yig = np.meshgrid(xi, yi)
    zi = interpolate.griddata((x,y), z, (xig, yig), method='cubic')


    return xig, yig, zi

def pressure_in_exp(pore_pressure_mat):

    f = h5py.File(pore_pressure_mat, 'r')

    pressure_dict = {}
    keys_list = list(f.keys())

    for elem in keys_list:
        k_value = np.array(f.get(elem))
        pressure_dict[elem] = k_value

    x = pressure_dict['xp']/1000
    y = pressure_dict['yp']/1000
    p_in_time = pressure_dict['P'].transpose()*10**6
    print(np.shape(p_in_time))
    Pinj = np.ones((20621,1))*1500000
    Pprod = np.ones((20621,1))*0

    #x = np.delete(x, [13])
    #y = np.delete(y, [13])
    p_in_time = np.delete(p_in_time,0,1)
    print(np.shape(p_in_time))

    x = np.append(x, 0.121)
    x = np.append(x, -0.121)
    y = np.append(y, 0.121)
    y = np.append(y, -0.121)
    p_in_time = np.append(p_in_time, Pinj, 1)
    p_in_time = np.append(p_in_time, Pprod, 1)

    x = x.reshape(1, 15)
    y = y.reshape(1, 15)
    p_in_time = p_in_time.reshape(20621, 15)

    xi = np.linspace(x.min(), x.max(), 100)
    yi = np.linspace(y.min(), y.max(), 100)
    xig_2, yig_2 = np.meshgrid(xi, yi)

    return x, y, p_in_time, xig_2, yig_2

fig2 = plt.figure()
xig,yig,zi = pressure_in_abaqus('por_vs_coord_2_exp_abaqus_stat.txt')
levels = list(range(0,1600000,50000))
surf = plt.contourf(xig, yig, zi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(zi), vmax=np.nanmax(zi),linewidth=0.2, levels = levels)
fig2.colorbar(surf, shrink=0.5, aspect=5)

x, y, p_in_time, xig_2, yig_2 = pressure_in_exp('P Exp3_2_1.mat')

def animate(i):

    zi_2 = interpolate.griddata((x[0], y[0]), p_in_time[100*i], (xig_2, yig_2), method='cubic')
    im = plt.contourf(xig_2, yig_2, zi_2, cmap=cm.jet, antialiased=True, vmin=np.nanmin(zi), vmax=np.nanmax(zi), linewidth=0.2, levels = levels)
    plt.scatter(x,y)
    print(i)
    return im,

im_ani = animation.FuncAnimation(fig2, animate, 206, interval=1000, blit=False, repeat=False)
im_ani.save('animation_cubic_interpol_St3Exp2.mp4', writer='ffmpeg')

#plt.show()