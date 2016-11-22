from scipy import interpolate
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt

from from_work.piezo_for_wells_and_frac_new_BC import PorePressure_in_Time



alpha = 0.8 * 10 ** -12
beta = 0.17 * 10 ** -9
hx = 0.005
hy = 0.005
hz = 0.07

t_step = 0.05
T_exp = 200
Lx = 0.5
Ly = 0.5

N = int(Lx / hx)  # количество ячеек вдоль оси х
M = int(Ly / hy)
print(N, M)

wells_with_Q = {}

# wells_with_Q = {(int((Lx / 2) / hx), int((Ly / 2 - 0.121) / hy)): -0.00000003}

frac_with_P = {}

Pres = 1 * 10 ** 5  # давление в пласте

V = hx * hy * hz
coeff_1 = hx * hz / hy
coeff_2 = hy * hz / hx
Pres_distrib = np.ones((N, M)) * Pres
Pprod = 100000
Pinj = 100000

f = h5py.File('data3-7.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value

x = pressure_dict['xp']/1000
y = pressure_dict['yp']/1000
pt = pressure_dict['p']*10**6
xhole = pressure_dict['xhole']/1000
yhole = pressure_dict['yhole']/1000

t = pressure_dict['t']
print(np.shape(y))
print(np.shape(pt))

x_list = []
y_list = []

for elem in x[0]:
    x_list.append(elem)

for elem in y[0]:
    y_list.append(elem)

set_coord = list(zip(x_list, y_list))
print(set_coord)


def pressure_in_exp(Pprod, Pinj, Lx, Ly, set_coord, t_current):

    set_with_P = {}
    i = 0
    for couple in set_coord:
        set_with_P[couple] = pt[i][t_current]
        i += 1


    set_with_P[(xhole[0][0], yhole[0][0])] = Pprod
    set_with_P[(xhole[2][0], yhole[2][0])] = Pinj

    set_with_P_mesh = {}
    for key in set_with_P:
        set_with_P_mesh[(int((key[0] + Lx / 2) / hx), int((key[1] + Ly / 2) / hy))] = set_with_P[key]

    wells_with_P = {}
    wells_with_P[(xhole[0][0], yhole[0][0])] = Pprod
    wells_with_P[(xhole[2][0], yhole[2][0])] = Pinj
    wells_with_P_mesh = {}
    for key in wells_with_P:
        wells_with_P_mesh[(int((key[0]+Lx/2)/hx), int((key[1]+Ly/2)/hy))] = wells_with_P[key]


    return set_with_P_mesh, set_with_P, wells_with_P_mesh


X = np.zeros((N, M))
Y = np.zeros((N, M))
for m in range(M):
    for n in range(N):
        X[n][m] = n * hx
        Y[n][m] = m * hy

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]
xi = np.linspace(min(X_list), max(X_list), 1000)
yi = np.linspace(min(Y_list), max(Y_list), 1000)
xig, yig = np.meshgrid(xi, yi)
levels = list(range(0, 1600000, 50000))

def animate(t_current):
    set_with_P_mesh, set_with_P, wells_with_P_mesh = pressure_in_exp(Pprod, Pinj, Lx, Ly, set_coord, t_current*10)
    print(t_current)

    global Pres_distrib
    print(np.min(Pres_distrib), np.max(Pres_distrib))
    P_total = PorePressure_in_Time(alpha, beta, t_step, N, M, wells_with_Q, set_with_P_mesh, frac_with_P, Pres, V, coeff_1, coeff_2, Pres_distrib)
    Pres_distrib = P_total


    P_list = [k for k in P_total.flat]
    CP_list = zip(X_list, Y_list, P_list)
    Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')
    #print(np.min(P_list), np.max(P_list))


    im = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=0, vmax=1600000,
                        linewidth=0.2, levels=levels)

    for key in set_with_P:
        plt.scatter(key[0] + Lx / 2, key[1] + Ly / 2, c='r')


    return im,


   ## surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
    #fig.colorbar(surf, shrink=0.5, aspect=5)

fig = plt.figure()

circle = np.arange(0, 2 * np.pi, 0.01)
r = 0.215
plt.plot(r * np.sin(circle) + Lx / 2, r * np.cos(circle) + Ly / 2, 'r')
plt.axis('equal')

im_ani = animation.FuncAnimation(fig, animate, 804, interval=100, blit=False, repeat=False)
#im_ani.save('animation_st2_exp7_at1.mp4', writer='ffmpeg')

plt.show()



