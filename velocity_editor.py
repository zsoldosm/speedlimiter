import subprocess
import os
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import math

#-----------HASZNÁLATHOZ SZÜKSÉGES VÁLTOZÓK MEGADÁSA ITT!-----------
filename = "egyetem_parkolo_teteje.csv"
limit = 0.9     # [m/s^2]
#-------------------------------------------------------------------

print("Futtatáshoz használt fájl: %s" % (filename))
print("Használt gyorsulási limit: %f m/s^2" % (limit))
print("A diagramok felugró ablakokban jelennek meg!")

#----------------ALPROGRAMOK-----------------
def data(filename):
    names = ["x", "y", "z", "yaw", "velocity", "change_flag"]
    df = pd.read_csv(filename, sep=",", names=names, skiprows=1)
    x = len(df) - 1
    df._set_value(0, "velocity", 0.0)
    df._set_value(x, "velocity", 0.0)
    return df

def delta_speed(filename):
    ds = []
    for i in range(1, len(data(filename))):
        ds.append((data(filename)["velocity"][i])-(data(filename)["velocity"][i-1]))
        ds[i-1] = round(ds[i-1], 4)
    return ds


def delta_dist(filename):
    dd = []
    dat = data(filename)
    for j in range(1, len(dat)):
        s1 = abs(dat["x"][j] - (dat["x"][j-1]))
        s0 = abs(dat["y"][j] - (dat["y"][j-1]))
        del_s = round(math.sqrt(s1**2 + s0**2), 4)
        dd.append(del_s)
    return dd

def distance(filename):
    d = []
    del_d = delta_dist(filename)
    for item in del_d:
        d.append(item)
    for i in range(0, len(del_d)):
        for j in range(0,i):
            d[i] += del_d[j] 
        d[i] = round(d[i], 4)
        a = [0]
    return a + d

def delta_time(filename):
    v = (data(filename)["velocity"])
    v.pop(0)
    s = delta_dist(filename)
    t = []
    for i in range(1, len(v)):
        x = round((s[i-1] / v[i])*3.6, 4)   #km/h -> m/s
        t.append(x)
    cor = (s[len(s)-1] / (v[len(s)-1]))*3.6
    cor = round(cor, 4)
    a = [cor]
    return t + a

def time(filename):
    ti = []
    del_t = delta_time(filename)
    for item in del_t:
        ti.append(item)
    for i in range(0, len(del_t)):
        for j in range(0,i):
            ti[i] += del_t[j] 
        ti[i] = round(ti[i], 4)
    a = [0]
    return a + ti

def del_accelerate(filename):
    v = delta_speed(filename)
    t = delta_time(filename)
    acc = []
    for i in range(0,len(v)):
        x = (v[i]/3.6) / t[i]
        x = round(x, 4)
        acc.append(x)
    return acc

def plot_vt(filename, limit):

    x1 = time(filename)
    x2 = time(filename)

    y1 = data(filename)["velocity"]
    y2 = limited_velo(filename, limit)

    fig, (ax1, ax2) = plt.subplots(2, 1)
    fig.suptitle('Original and limited velocity')

    ax1.plot(x1, y1, 'x-b')
    ax1.set_ylabel('Original [km/h]')

    ax2.plot(x2, y2, 'x-g')
    ax2.set_xlabel('time (s)')
    ax2.set_ylabel('Limited [km/h]')
    print("")

def limit_acc(filename, limit):
    acc = del_accelerate(filename)
    for i in range(0, len(acc)):
        if acc[i] >= limit:
            acc[i] = limit
        elif acc[i] <= (limit * -1):
            acc[i] = (limit * -1)
    return acc

def limited_delta_velo(filename, limit):
    a = limit_acc(filename, limit)
    t = delta_time(filename)
    velo_d = []
    for i in range(0, len(a)):
        v = round(((a[i] * t[i])*3.6), 4)
        velo_d.append(v)
    return velo_d

def limited_velo(filename, limit):
    velo = []
    del_v = limited_delta_velo(filename, limit)
    for item in del_v:
        velo.append(item)
    for i in range(0, len(del_v)):
        for j in range(0,i):
            velo[i] += del_v[j] 
        velo[i] = round(velo[i], 4)
    a = [0]
    velo[len(del_v)-1] = 0
    return a + velo

def plot_at(filename, limit):

    x1 = time(filename)[:-1]
    x2 = time(filename)[:-1]

    y1 = del_accelerate(filename)
    y2 = limit_acc(filename, limit)

    fig, (ax1, ax2) = plt.subplots(2, 1)
    fig.suptitle('Original and limited accelerate')

    ax1.plot(x1, y1, 'ob')
    ax1.set_ylabel('Original [m/s^2]')

    ax2.plot(x2, y2, 'og')
    ax2.set_xlabel('time (s)')
    ax2.set_ylabel('Limited [m/s^2]')
    return plt.show()

#-----------------SZÜKSÉGES ADATOK FÁJLBA KIÍRÁSA------------------
lim = np.asarray(limited_velo(filename, limit), order='F')
df1 = pd.DataFrame({
    'Original': data(filename)["velocity"],
    'Limited': lim
})
fn = filename[:-4] + "_compare_velo.csv"
df1.to_csv(fn, index=False)
print("Eredeti és limitált sebesség-értékek fájlba kiírva!")
ora = np.asarray(del_accelerate(filename), order='F')
lima = np.asarray(limit_acc(filename, limit), order='F')
df2 = pd.DataFrame({
    'Original': ora,
    'Limited': lima
})
fn_a = filename[:-4] + "_compare_acc.csv"
df2.to_csv(fn_a, index=False)
print("Eredeti és limitált gyorsulás-értékek fájlba kiírva!")
plot_vt(filename, limit)
plot_at(filename, limit)

#------------SZÁMOLT ADATOK KIIRATÁSA, HA SZÜKSÉGES---------------
"""
print("\nEredeti sebesség értékek [km/h]:\n")
print(data(filename)["velocity"])
print("\nKét pont közötti sebesség különbség [km/h]:\n")
print(delta_speed(filename))
print("\nKét pont között megtett út [m]:\n")
print(delta_dist(filename))
print("\nAdott pontig összesen megtett út [m]:\n")
print(distance(filename))
print("\n2 pont között eltelt idő [sec]:\n")
print(delta_time(filename))
print("\nAdott pontig összesen eltelt idő [sec]:\n")
print(time(filename))
print("\nAdott pontban számolt gyorsulás értékek [m/s^2]:\n")
print(del_accelerate(filename))
#print(limit_acc(filename, 0.98))
#print(limited_delta_velo(filename, 0.98))
#print(data(filename)["velocity"])
print("\nLimitált gyorsulásból sebesség értékek [km/h]:\n")
print(limited_velo(filename, limit))
"""