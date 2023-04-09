# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


import numpy as np
from numpy.linalg import inv
from numpy.linalg import det
from numpy.linalg import solve
import pandas as pd
from openpyxl.workbook import Workbook

E1=55e6
E2=35e6
E3=23e9

v1=0.22
v2=0.22
v3=0.25

d1=2300
d2=2300
dw=1000

water_level=290




def get_data():
    x=np.matrix
    y=np.matrix
    assoc=np.matrix
    input = np.loadtxt("xydata.txt", dtype='f', delimiter=',')
    x=input[:,1]
    y=input[:,2]
    input1 = np.loadtxt("assoc.txt", dtype='f', delimiter=',')
    assoc=input1[:,1:4]
    nds = np.prod(x.shape)
    els = int(np.prod(assoc.shape)/3)
    return x, y, assoc, nds, els


def get_nodes_of_materials(xc, yc, nds):
    ige1 = np.array
    ige2 = np.array
    base = np.array
    kol = 0
    kol1 = 0
    kol2 = 0
    for i in range(0, nds):
        xi = xc[i]
        yi = yc[i]
        if ((xi < 9 / 5 * (yi + 71 / 3) or abs(yi - 5 / 9 * xi + 71 / 3) < pow(10, -4)) and yi >= 0 and xi >= 0):
            ige1 = np.append(ige1, i)
            kol = kol + 1
        if ((xi > 9 / 5 * (yi + 71 / 3) or abs(yi - 5 / 9 * xi + 71 / 3) < pow(10, -4)) and yi >= 0 and xi <= 1029):
            ige2 = np.append(ige2, i)
            kol1 = kol1 + 1
        if (yi <= 0):
            base = np.append(base, i)
            kol2 = kol2 + 1
    return ige1,ige2,base,kol,kol1,kol2


def get_elements_of_materials(xc, yc, assoc, els):
    ige1mt = np.array
    ige2mt = np.array
    basemt = np.array
    kol = 0
    kol1 = 0
    kol2 = 0
    for i in range(0, els):
        xi = xc[int(assoc[i, 0]-1)]
        yi = yc[int(assoc[i, 0]-1)]

        xj = xc[int(assoc[i, 1] - 1)]
        yj = yc[int(assoc[i, 1] - 1)]

        xk = xc[int(assoc[i, 2] - 1)]
        yk = yc[int(assoc[i, 2] - 1)]

        xcm = 1/3*(xi + xj + xk)
        ycm = 1/3*(yi + yj + yk)

        if ((xcm < 9 / 5 * (ycm + 71 / 3) ) and ycm >= 0 and xcm >= 0):
            ige1mt = np.append(ige1mt, i)
            kol = kol + 1
        if ((xcm > 9 / 5 * (ycm + 71 / 3)) and ycm >= 0 and xcm <= 1029):
            ige2mt = np.append(ige2mt, i)
            kol1 = kol1 + 1
        if (ycm <= 0):
            basemt = np.append(basemt, i)
            kol2 = kol2 + 1

    return ige1mt, ige2mt, basemt, kol, kol1, kol2


def get_square(xi,yi,xj,yj,xk,yk):
    Sqm = np.zeros((3, 3))
    Sqm[:, 2] = 1
    Sqm[0, 0] = xi
    Sqm[0, 1] = yi
    Sqm[1, 0] = xj
    Sqm[1, 1] = yj
    Sqm[2, 0] = xk
    Sqm[2, 1] = yk
    S = 1 / 2 * abs(det(Sqm))
    return S



def k_matrix1(ige1m, assoc, nds):
    #для основание
    kol_els_1 = int(np.prod(ige1m.shape))
    J = np.matrix
    D = np.matrix
    K = np.zeros((2*nds, 2*nds))
    #print(kol_els_1)
    for i in range(1, kol_els_1):
        xi = xc[int(assoc[ige1m[i], 0] - 1)]
        yi = yc[int(assoc[ige1m[i], 0] - 1)]

        xj = xc[int(assoc[ige1m[i], 1] - 1)]
        yj = yc[int(assoc[ige1m[i], 1] - 1)]

        xk = xc[int(assoc[ige1m[i], 2] - 1)]
        yk = yc[int(assoc[ige1m[i], 2] - 1)]

        D = E1 * (1 - v1) / ((1 + v1) * (1 - 2 * v1)) * np.matrix([[1, v1 / (1 - v1), 0], [v1 / (1 - v1), 1, 0],[0, 0, (1 - 2 * v1) / (2 * (1 - v1))]])

        der = np.matrix([[-1, 0, 1], [-1, 1, 0]])
        J = der*np.matrix([[xi, yi], [xj, yj], [xk, yk]])

        B=np.zeros((3,6))
        invJ = inv(J)
        tmp = invJ*der[:, 0]
        B[0, 0] = tmp[0]
        B[2, 0] = tmp[1]
        B[1, 1] = tmp[1]
        B[2, 1] = tmp[0]
        tmp = invJ * der[:, 1]
        B[0, 2] = tmp[0]
        B[2, 2] = tmp[1]
        B[1, 3] = tmp[1]
        B[2, 3] = tmp[0]
        tmp = invJ * der[:,2]
        B[0, 4] = tmp[0]
        B[2, 4] = tmp[1]
        B[1, 5] = tmp[1]
        B[2, 5] = tmp[0]

        S = get_square(xi, yi, xj, yj, xk, yk)

        ke = np.transpose(B)*D*B*S

        A = np.zeros((6, 2*nds))

        it = int(assoc[i, 0] - 1)
        jt = int(assoc[i, 1] - 1)
        kt = int(assoc[i, 2] - 1)

        A[0, 2*it] = 1
        A[1, 2*it+1] = 1
        A[2, 2 * jt] = 1
        A[3, 2 * jt + 1] = 1
        A[4, 2 * kt] = 1
        A[5, 2 * kt + 1] = 1

        K = K + np.transpose(A)*ke*A

    return K

def k_matrix2(ige1m, assoc, nds):
    #для основание
    kol_els_1 = int(np.prod(ige1m.shape))
    J = np.matrix
    D = np.matrix
    K = np.zeros((2*nds, 2*nds))
    #print(kol_els_1)
    for i in range(1, kol_els_1):
        xi = xc[int(assoc[ige1m[i], 0] - 1)]
        yi = yc[int(assoc[ige1m[i], 0] - 1)]

        xj = xc[int(assoc[ige1m[i], 1] - 1)]
        yj = yc[int(assoc[ige1m[i], 1] - 1)]

        xk = xc[int(assoc[ige1m[i], 2] - 1)]
        yk = yc[int(assoc[ige1m[i], 2] - 1)]

        D = E2 * (1 - v2) / ((1 + v2) * (1 - 2 * v2)) * np.matrix([[1, v2 / (1 - v2), 0], [v2 / (1 - v2), 1, 0],[0, 0, (1 - 2 * v2) / (2 * (1 - v2))]])

        der = np.matrix([[-1, 0, 1], [-1, 1, 0]])
        J = der*np.matrix([[xi, yi], [xj, yj], [xk, yk]])

        B=np.zeros((3,6))
        invJ = inv(J)
        tmp = invJ*der[:, 0]
        B[0, 0] = tmp[0]
        B[2, 0] = tmp[1]
        B[1, 1] = tmp[1]
        B[2, 1] = tmp[0]
        tmp = invJ * der[:, 1]
        B[0, 2] = tmp[0]
        B[2, 2] = tmp[1]
        B[1, 3] = tmp[1]
        B[2, 3] = tmp[0]
        tmp = invJ * der[:,2]
        B[0, 4] = tmp[0]
        B[2, 4] = tmp[1]
        B[1, 5] = tmp[1]
        B[2, 5] = tmp[0]

        S = get_square(xi, yi, xj, yj, xk, yk)

        ke = np.transpose(B)*D*B*S

        A = np.zeros((6, 2*nds))

        it = int(assoc[i, 0] - 1)
        jt = int(assoc[i, 1] - 1)
        kt = int(assoc[i, 2] - 1)

        A[0, 2*it] = 1
        A[1, 2*it+1] = 1
        A[2, 2 * jt] = 1
        A[3, 2 * jt + 1] = 1
        A[4, 2 * kt] = 1
        A[5, 2 * kt + 1] = 1

        K = K + np.transpose(A)*ke*A

    return K


def k_matrix3(ige1m, assoc, nds):
    #для основание
    kol_els_1 = int(np.prod(ige1m.shape))
    J = np.matrix
    D = np.matrix
    K = np.zeros((2*nds, 2*nds))
    #print(kol_els_1)
    for i in range(1, kol_els_1):
        xi = xc[int(assoc[ige1m[i], 0] - 1)]
        yi = yc[int(assoc[ige1m[i], 0] - 1)]

        xj = xc[int(assoc[ige1m[i], 1] - 1)]
        yj = yc[int(assoc[ige1m[i], 1] - 1)]

        xk = xc[int(assoc[ige1m[i], 2] - 1)]
        yk = yc[int(assoc[ige1m[i], 2] - 1)]

        D = E3 * (1 - v3) / ((1 + v3) * (1 - 2 * v3)) * np.matrix([[1, v3 / (1 - v3), 0], [v3 / (1 - v3), 1, 0],[0, 0, (1 - 2 * v3) / (2 * (1 - v3))]])

        der = np.matrix([[-1, 0, 1], [-1, 1, 0]])
        J = der*np.matrix([[xi, yi], [xj, yj], [xk, yk]])

        B=np.zeros((3,6))
        invJ = inv(J)
        tmp = invJ*der[:, 0]
        B[0, 0] = tmp[0]
        B[2, 0] = tmp[1]
        B[1, 1] = tmp[1]
        B[2, 1] = tmp[0]
        tmp = invJ * der[:, 1]
        B[0, 2] = tmp[0]
        B[2, 2] = tmp[1]
        B[1, 3] = tmp[1]
        B[2, 3] = tmp[0]
        tmp = invJ * der[:,2]
        B[0, 4] = tmp[0]
        B[2, 4] = tmp[1]
        B[1, 5] = tmp[1]
        B[2, 5] = tmp[0]

        S = get_square(xi, yi, xj, yj, xk, yk)

        ke = np.transpose(B)*D*B*S

        A = np.zeros((6, 2*nds))

        it = int(assoc[i, 0] - 1)
        jt = int(assoc[i, 1] - 1)
        kt = int(assoc[i, 2] - 1)

        A[0, 2*it] = 1
        A[1, 2*it+1] = 1
        A[2, 2 * jt] = 1
        A[3, 2 * jt + 1] = 1
        A[4, 2 * kt] = 1
        A[5, 2 * kt + 1] = 1

        K = K + np.transpose(A)*ke*A

    return K


def gravity_force(ige1, ige2, k1, k2,xc, yc, assoc, els, fe):
    for i in range(1, k1+1):
        tmp = np.array([])
        a = ige1[i]
        for j in range(0, els):
            it = int(assoc[j, 0] - 1)
            jt = int(assoc[j, 1] - 1)
            kt = int(assoc[j, 2] - 1)
            if(a == it or a==jt or a==kt):
                tmp = np.append(tmp, j)
        V = 0
        for k in range(0 , int(np.prod(tmp.shape))-1):
            xi = xc[int(assoc[k, 0] - 1)]
            yi = yc[int(assoc[k, 0] - 1)]

            xj = xc[int(assoc[k, 1] - 1)]
            yj = yc[int(assoc[k, 1] - 1)]

            xk = xc[int(assoc[k, 2] - 1)]
            yk = yc[int(assoc[k, 2] - 1)]
            V = V + get_square(xi, yi, xj, yj, xk, yk)
        fe[2*a+1] = fe[2*a+1] - d1*9.8*V
    for i in range(1, k2):
        tmp = np.array([])
        a = ige2[i]
        for j in range(0,els):
            it = int(assoc[j, 0] - 1)
            jt = int(assoc[j, 1] - 1)
            kt = int(assoc[j, 2] - 1)
            if(a == it or a==jt or a==kt):
                tmp = np.append(tmp, j)
        V = 0
        for k in range(0 , int(np.prod(tmp.shape))-1):
            xi = xc[int(assoc[k, 0] - 1)]
            yi = yc[int(assoc[k, 0] - 1)]

            xj = xc[int(assoc[k, 1] - 1)]
            yj = yc[int(assoc[k, 1] - 1)]

            xk = xc[int(assoc[k, 2] - 1)]
            yk = yc[int(assoc[k, 2] - 1)]
            V = V + get_square(xi, yi, xj, yj, xk, yk)
        fe[2*a+1] = fe[2*a+1] - d1*9.8*V

    return fe


def get_bot_bc(yc, nds):
    kol=0
    bbot=np.array([])
    for i in range(0, nds):
        if(yc[i] == -293):
            kol = kol+1
            bbot = np.append(bbot, i)
    return bbot, kol
def get_right_left_bc(xc, nds):
    kol1 = 0
    kol2 = 0
    bleft = np.array([])
    bright = np.array([])
    for i in range(0, nds):
        if (abs(xc[i]+1028.4)<pow(10,-4)):
            kol1 = kol1 + 1
            bleft = np.append(bleft, i)
        elif (abs(xc[i]-2056.8)<pow(10,-4)):
            kol2 = kol2 + 1
            bright = np.append(bright, i)
    return bleft, bright, kol1, kol2


def set_bc(bleft, k1, bright, k2, bbot, k3, m):
    k = m
    for i in range(0, k1):
        k[2 * int(bleft[i]), :] = 0
        k[:, 2 * int(bleft[i])] = 0
        k[2 * int(bleft[i]), 2 * int(bleft[i])] = 1
    for i in range(0, k2):
        k[2 * int(bright[i]), :] = 0
        k[:, 2 * int(bright[i])] = 0
        k[2 * int(bright[i]), 2 * int(bright[i])] = 1
    for i in range(0, k3):
        k[2 * int(bbot[i])+1, :] = 0
        k[:, 2 * int(bbot[i])+1] = 0
        k[2 * int(bbot[i])+1 , 2 * int(bbot[i])+1] = 1

    return k

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    xc = np.matrix
    yc = np.matrix
    assoc = np.matrix
    xc, yc, assoc, nds, els = get_data()
    print("nodes "+str(nds) + " els "+str(els))
#считаем матрицу K

#получвем узлы, принадлежащие материалам
    ige1, ige2, base, kolig1n, kolig2n, kolbasen = get_nodes_of_materials(xc, yc, nds)
    ige1m,ige2m,basem, kolig1e, kolig2e, kolbasee = get_elements_of_materials(xc, yc, assoc, els)
    #print(ige1, ige2)
    #print(kolig1n, kolig2n, kolbasen)
    #print(kolig1e, kolig2e, kolbasee)


    Kige1 = k_matrix1(ige1m, assoc, nds)

    Kige2 = k_matrix2(ige2m, assoc, nds)

    Kbase = k_matrix3(basem, assoc, nds)
    #print(nds)
    K = Kige1 + Kige2 + Kbase
    fe = np.zeros((2*nds, 1))
    print(base)
    # составляем нагрузки
    fe = gravity_force(ige1, ige2, kolig1n, kolig2n, xc, yc, assoc, els, fe)
    # получаем узла под гу
    bcbot,kolbot=get_bot_bc(yc, nds)
    bcleft, bcright, kolleft,kolright = get_right_left_bc(xc, nds)
    # уставнавливаем ГУ
    K = set_bc(bcleft,  kolleft, bcright, kolright, bcbot, kolbot, K)

    print(K[22, 22])
    df = pd.DataFrame(K)

    # save to xlsx file

    #filepath = 'my_excel_file.xlsx'

    #df.to_excel(filepath, index=False)
    u = solve(K, fe)







# See PyCharm help at https://www.jetbrains.com/help/pycharm/
