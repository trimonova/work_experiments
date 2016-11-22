#  Можно задавать много точек с источниками и давлениями. Можно задавать распределение давления в трещине. Если задавать только давления (граничные условия) то задача устойчива при любых шагах времени и координаты,
# если задавать еще источники, то задача устойчива при каком-то соотношении t_step и hx, граничное условие-новое, градиента давления на границе равен 0.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# если в скважинах задавать расход, а не давление, то задача НЕ всегда устойчива. Нужно подбирать шаги повремени и по координатам.
if __name__ == '__main__':
    alpha = 0.8*10**-12
    beta = 0.17*10**-9
    hx = 0.01
    hy = 0.01
    hz = 0.07

    t_step = 0.01
    T_exp = 100
    Lx = 0.5
    Ly = 0.5

    N = int(Lx/hx) # количество ячеек вдоль оси х
    M = int(Ly/hy)
    print(N,M)

    #wells_with_Q = {}
    wells_with_P = {(int((Lx/2+0.121)/hx),int((Ly/2+0.121)/hy)): 20*10**5, (int((Lx/2-0.121)/hx), int((Ly/2-0.121)/hy)): 1*10**5}
    frac_with_P = {(int(N/2), int(M/2)):25*10**5, (int(N/2)-1, int(M/2)):24*10**5, (int(N/2)+1, int(M/2)):24*10**5, (int(N/2)-2, int(M/2)):23*10**5, (int(N/2)+2, int(M/2)):23*10**5, (int(N/2)-3, int(M/2)):22*10**5, (int(N/2)+3, int(M/2)):22*10**5, (int(N/2)+4, int(M/2)):21*10**5, (int(N/2)-4, int(M/2)):21*10**5}
    wells_with_Q = {(int((Lx/2)/hx),int((Ly/2-0.121)/hy)): -0.00000003}
    #wells_with_P = {(int((Lx/2+0.121)/hx),int((Ly/2+0.121)/hy)): 20*10**5, (int((Lx/2-0.121)/hx), int((Ly/2-0.121)/hy)): 1*10**5, (int((Lx/2-0.121)/hx), int((Ly/2)/hy)): 5*10**5, (int((Lx/2+0.121)/hx), int((Ly/2)/hy)): 5*10**5, (int((Lx/2)/hx), int((Ly/2-0.121)/hy)): 5*10**5, (int((Lx/2)/hx), int((Ly/2+0.121)/hy)): 5*10**5}
    #frac_with_P = {}

    Pres = 1*10**5 # давление в пласте

    V = hx*hy*hz
    coeff_1 = hx*hz/hy
    coeff_2 = hy*hz/hx
    Pres_distrib = np.ones((N, M)) * Pres

def PorePressure_in_Time(alpha, beta, t_step, N, M, wells_with_Q, wells_with_P, frac_with_P, Pres, V, coeff_1, coeff_2, Pres_distrib):
      # пластовое давление во всей области на нулевом временном шаге
    indic = []
    P_total = np.ones((N, 1)) * Pres
    for m in range(0,M):
        A = np.zeros((N,N))
        B = np.zeros((N,1))

        for n in range(1, N-1):
            A[n][n-1] = alpha*coeff_2
            A[n][n] = (-2*coeff_2*alpha - V*beta/t_step)
            A[n][n+1] = alpha*coeff_2

        A[0][0] = -2*coeff_2*alpha - V*beta/t_step
        A[0][1] = 2*alpha*coeff_2
        A[N-1][N-1] = A[0][0]
        A[N-1][N-2] = A[0][1]

        for n in range(0,N):
            if m == 0:
                B[n][0] = -V*beta/t_step*Pres_distrib[n][m]- alpha*coeff_1*(- 2*Pres_distrib[n][m] + 2*Pres_distrib[n][m+1])

            elif m == M-1:
                B[n][0] = -V*beta/t_step*Pres_distrib[n][m]- alpha*coeff_1*(2*Pres_distrib[n][m-1] - 2*Pres_distrib[n][m])

            else:
                B[n][0] = -V*beta/t_step*Pres_distrib[n][m]- alpha*coeff_1*(Pres_distrib[n][m-1] - 2*Pres_distrib[n][m] + Pres_distrib[n][m+1])

            for coord_key in wells_with_Q:
                if (n,m) == coord_key:
                    B[n][0] = -V*beta/t_step*Pres_distrib[n][m]- alpha*coeff_1*(Pres_distrib[n][m-1] - 2*Pres_distrib[n][m] + Pres_distrib[n][m+1]) + wells_with_Q[coord_key]

        for n in range(0, N):

            for coord_key in wells_with_P:
                if (n,m) == coord_key:
                    indic.append(coord_key)
                elif (n-1,m) == coord_key:
                    A[n][n - 1] = 0
                    B[n][0] = B[n][0] - alpha * coeff_2 * wells_with_P[coord_key]
                elif (n+1,m) == coord_key:
                    A[n][n+1] = 0
                    B[n][0] = B[n][0] - alpha * coeff_2 * wells_with_P[coord_key]

            for coord_key_fr in frac_with_P:
                if (n,m) == coord_key_fr:
                    indic.append(coord_key_fr)
                elif (n-1,m) == coord_key_fr:
                    A[n][n - 1] = 0
                    B[n][0] = B[n][0] - alpha * coeff_2 * frac_with_P[coord_key_fr]
                elif (n+1,m) == coord_key_fr:
                    A[n][n+1] = 0
                    B[n][0] = B[n][0] - alpha * coeff_2 * frac_with_P[coord_key_fr]

        #print(type(indic))
        if indic != []:
            counter = 0
            indic.sort()
            for element in indic:
                A = np.delete(A, element[0]-counter, axis=0)
                A = np.delete(A, element[0]-counter, axis=1)
                B = np.delete(B, element[0]-counter)
                counter += 1


        P_new = np.linalg.solve(A,B)
#
        if indic != []:

            for element in indic:
                if element in wells_with_P:
                    P_new = np.insert(P_new,element[0],wells_with_P[element])
                else:
                    P_new = np.insert(P_new, element[0], frac_with_P[element])

            indic = []
            P_new = P_new.reshape(N, 1)

        P_total = np.hstack((P_total,P_new))

    P_total = np.delete(P_total, 0, axis=1)

    Pres_distrib = np.array(P_total.copy())


#---------------------------------------------------------------------------
    indic = []

    P_total = np.ones((1, M)) * Pres
    for n in range(0, N):
        A = np.zeros((M, M))
        B = np.zeros((M, 1))
        for m in range(1, M - 1):
            A[m][m - 1] = alpha * coeff_1
            A[m][m] = (-2 * coeff_1 * alpha - V * beta / t_step)
            A[m][m + 1] = alpha * coeff_1
        A[0][0] = -2 * coeff_1 * alpha - V * beta / t_step
        A[0][1] = 2*alpha * coeff_1
        A[M - 1][M - 1] = A[0][0]
        A[M - 1][M - 2] = A[0][1]

        for m in range(0, M):
            if n == 0:
                B[m][0] = -V * beta / t_step * Pres_distrib[n][m] - alpha * coeff_2 * (-2*Pres_distrib[n][m] + 2*Pres_distrib[n+1][m])

            elif n == N-1:
                B[m][0] = -V * beta / t_step * Pres_distrib[n][m] - alpha * coeff_2 * (-2*Pres_distrib[n][m] + 2 * Pres_distrib[n - 1][m])

            else:
                B[m][0] = -V * beta / t_step * Pres_distrib[n][m] - alpha * coeff_2 *(Pres_distrib[n-1][m] - 2 * Pres_distrib[n][m] + Pres_distrib[n+1][m])

            for coord_key in wells_with_Q:
                if (n, m) == coord_key:
                    B[m][0] = -V * beta / t_step * Pres_distrib[n][m] - alpha * coeff_2 * (Pres_distrib[n][m - 1] - 2 * Pres_distrib[n][m] + Pres_distrib[n][m + 1]) + wells_with_Q[coord_key]


        for m in range(0, M):
            for coord_key in wells_with_P:
                if (n,m) == coord_key:
                    indic.append(coord_key)

                elif (n,m-1) == coord_key:
                    A[m][m - 1] = 0
                    B[m][0] = B[m][0] - alpha * coeff_1 * wells_with_P[coord_key]

                elif (n,m+1) == coord_key:
                    A[m][m + 1] = 0
                    B[m][0] = B[m][0] - alpha * coeff_1 * wells_with_P[coord_key]

            for coord_key_fr in frac_with_P:
                if (n,m) == coord_key_fr:
                    indic.append(coord_key_fr)
                elif (n,m-1) == coord_key_fr:
                    A[m][m - 1] = 0
                    B[m][0] = B[m][0] - alpha * coeff_1 * frac_with_P[coord_key_fr]
                elif (n,m+1) == coord_key_fr:
                    A[m][m+1] = 0
                    B[m][0] = B[m][0] - alpha * coeff_1 * frac_with_P[coord_key_fr]


        if indic != []:
            counter = 0
            indic.sort()
            for element in indic:
                A = np.delete(A, element[1]-counter, axis=0)
                A = np.delete(A, element[1]-counter, axis=1)
                B = np.delete(B, element[1]-counter)
                counter += 1

        P_new = np.linalg.solve(A,B)

        if indic != []:

            for element in indic:
                if element in wells_with_P:
                    P_new = np.insert(P_new, element[1], wells_with_P[element])
                else:
                    P_new = np.insert(P_new, element[1], frac_with_P[element])

            indic = []
            P_new = P_new.reshape(M, 1)


        P_total = np.vstack((P_total, P_new.T))

    P_total = np.delete(P_total, 0, axis=0)
    Pres_distrib = np.array(P_total.copy())


    return P_total

#----------------------------------------------------------------------------

if __name__ == '__main__':
    for t in range(T_exp):
        P_total = PorePressure_in_Time(alpha, beta, t_step, N, M, wells_with_Q, wells_with_P, frac_with_P, Pres, V, coeff_1, coeff_2, Pres_distrib)
        Pres_distrib = P_total

    X = np.zeros((N,M))
    Y = np.zeros((N, M))
    for m in range(M):
        for n in range(N):
            X[n][m] = n*hx
            Y[n][m] = m*hy

    X_list = [i for i in X.flat]
    Y_list = [j for j in Y.flat]
    P_list = [k for k in P_total.flat]


    CP_list = zip(X_list, Y_list, P_list)

if __name__ == '__main__':
    print(min(P_list), max(P_list))

    xi = np.linspace(min(X_list),max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')

    levels = list(range(0,2600000,50000))
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)

    t = np.arange(0, 2 * np.pi, 0.01)
    r = 0.215
    plt.plot(r * np.sin(t) + Lx/2, r * np.cos(t) + Ly/2)
    #ax = fig.gca(projection='3d')

    #surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()



















#