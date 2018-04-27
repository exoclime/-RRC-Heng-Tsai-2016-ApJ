# Python script to compute atmospheric chemistry (C-H-O-N system with 9 molecules)
# by Kevin Heng (8th March 2016)
# modified by Shang-Min Tsai (7th June 2016)
# This script is made to reproduce the three plots shown in Figure 1. in Heng & Tsai (2016)
# note: K, K2, K3 have Gibbs free energies expressed in J/K/mol.
#       For K4, K5, K6, K7, I have used kJ/K/mol.

# choose to play 1a, 1b or 1c (corresponding to the top, middle, and bottom panels in Figure 2.)
plot_fig = '2a'

from numpy import mean,arange,zeros,polynomial,array,interp,exp,sqrt
from matplotlib import pyplot as plt
import numpy as np

# function to compute first equilibrium constant (K)
def kk(my_temperature,pbar):
    runiv = 8.3144621   # J/K/mol
    temperatures = arange(500.0, 3100.0, 100.0)
    dg = [96378.0, 72408.0, 47937.0, 23114.0, -1949.0, -27177.0, -52514.0, -77918.0, -103361.0, -128821.0, -154282.0, -179733.0, -205166.0, -230576.0, -255957.0, -281308.0, -306626.0, -331911.0, -357162.0, -382380.0, -407564.0, -432713.0, -457830.0, -482916.0, -507970.0, -532995.0]
    my_dg = interp(my_temperature,temperatures,dg)
    result = exp(-my_dg/runiv/my_temperature)/pbar/pbar
    return result

# function to compute second equilibrium constant (K2)
def kk2(my_temperature):
    runiv = 8.3144621   # J/K/mol
    temperatures = arange(500.0, 3100.0, 100.0)
    dg2 = [20474.0, 16689.0, 13068.0, 9593.0, 6249.0, 3021.0, -107.0, -3146.0, -6106.0, -8998.0, -11828.0, -14600.0, -17323.0, -20000.0, -22634.0, -25229.0, -27789.0, -30315.0, -32809.0, -35275.0, -37712.0, -40123.0, -42509.0, -44872.0, -47211.0, -49528.0]
    my_dg = interp(my_temperature,temperatures,dg2)
    result = exp(-my_dg/runiv/my_temperature)
    return result

# function to compute second equilibrium constant (K3)
def kk3(my_temperature,pbar):
    runiv = 8.3144621   # J/K/mol
    temperatures = arange(500.0, 3100.0, 100.0)
    dg3 = [262934.0, 237509.0, 211383.0, 184764.0, 157809.0, 130623.0, 103282.0, 75840.0, 48336.0, 20797.0, -6758.0, -34315.0, -61865.0, -89403.0, -116921.0, -144422.0, -171898.0, -199353.0, -226786.0, -254196.0, -281586.0, -308953.0, -336302.0, -363633.0, -390945.0, -418243.0]
    my_dg = interp(my_temperature,temperatures,dg3)
    result = exp(-my_dg/runiv/my_temperature)/pbar/pbar
    return result

# function to compute second equilibrium constant (K4)
def kk4(my_temperature,pbar):
    runiv = 8.3144621e-3   # kJ/K/mol
    temperatures = arange(500.0, 3100.0, 100.0)
    dg4 = [116.519, 103.718, 90.63, 77.354, 63.959, 50.485, 36.967, 23.421, 9.864, -3.697, -17.253, -30.802, -44.342, -57.87, -71.385, -84.888, -98.377, -111.855, -125.322, -138.777, -152.222, -165.657, -179.085, -192.504, -205.916, -219.322]
    my_dg = interp(my_temperature,temperatures,dg4)
    result = exp(-my_dg/runiv/my_temperature)/pbar
    return result

# function to compute second equilibrium constant (K5)
def kk5(my_temperature,pbar):
    runiv = 8.3144621e-3   # kJ/K/mol
    temperatures = arange(500.0, 3100.0, 100.0)
    dg5 = [-9.6, -31.758, -54.38, -77.324, -100.494, -123.82, -147.25, -170.746, -194.282, -217.836, -241.392, -264.938, -288.468, -311.972, -335.45, -358.894, -382.304, -405.68, -429.018, -452.32, -475.584, -498.812, -522.006, -545.162, -568.286, -591.378]
    my_dg = interp(my_temperature,temperatures,dg5)
    result = exp(-my_dg/runiv/my_temperature)/pbar/pbar
    return result

# function to compute second equilibrium constant (K6)
def kk6(my_temperature,pbar):
    runiv = 8.3144621e-3   # kJ/K/mol
    temperatures = arange(500.0, 3100.0, 100.0)
    dg6 = [145.71, 121.401, 96.516, 71.228, 45.662, 19.906, -5.977, -31.942, -57.955, -83.992, -110.035, -136.073, -162.096, -188.098, -214.075, -240.023, -265.94, -291.826, -317.679, -343.5, -369.29, -395.047, -420.775, -446.472, -472.141, -497.784]
    my_dg = interp(my_temperature,temperatures,dg6)
    result = exp(-my_dg/runiv/my_temperature)/pbar/pbar
    return result

# "old" function to compute mixing ratio for CO (without nitrogen, from Heng & Lyons 2016)
def n_cmono_old(n_o,n_c,temp,pbar):
    k1 = kk(temp,pbar)
    k2 = kk2(temp)
    k3 = kk3(temp,pbar)
    a0 = 8.0*k1*k3*k3/k2
    a1 = 8.0*k1*k3/k2
    a2 = 2.0*k1/k2*( 1.0 + 8.0*k3*(n_o-n_c) ) + 2.0*k1*k3
    a3 = 8.0*k1/k2*(n_o-n_c) + 2.0*k3 + k1
    a4 = 8.0*k1/k2*(n_o-n_c)*(n_o-n_c) + 1.0 + 2.0*k1*(n_o-n_c)
    a5 = -2.0*n_c
    result = polynomial.polynomial.polyroots([a5,a4,a3,a2,a1,a0])
    n_ch4 = result[4]  # picks out correct root (computer-dependent)
    n_h2o = 2.0*k3*n_ch4*n_ch4 + n_ch4 + 2.0*(n_o-n_c)
    n_co = k1*n_ch4*n_h2o
    return n_co

# "old" function to compute mixing ratio for H2O (without nitrogen, from Heng & Lyons 2016)
def n_water_old(n_o,n_c,temp,pbar):
    k1 = kk(temp,pbar)
    k2 = kk2(temp)
    k3 = kk3(temp,pbar)
    a0 = 8.0*k1*k3*k3/k2
    a1 = 8.0*k1*k3/k2
    a2 = 2.0*k1/k2*( 1.0 + 8.0*k3*(n_o-n_c) ) + 2.0*k1*k3
    a3 = 8.0*k1/k2*(n_o-n_c) + 2.0*k3 + k1
    a4 = 8.0*k1/k2*(n_o-n_c)*(n_o-n_c) + 1.0 + 2.0*k1*(n_o-n_c)
    a5 = -2.0*n_c
    result = polynomial.polynomial.polyroots([a5,a4,a3,a2,a1,a0])
    n_ch4 = result[4]  # picks out correct root (computer-dependent)
    n_h2o = 2.0*k3*n_ch4*n_ch4 + n_ch4 + 2.0*(n_o-n_c)
    return n_h2o

# "old" function to compute mixing ratio for CH4 (without nitrogen, from Heng & Lyons 2016)
def n_methane_old(n_o,n_c,temp,pbar):
    k1 = kk(temp,pbar)
    k2 = kk2(temp)
    k3 = kk3(temp,pbar)
    a0 = 8.0*k1*k3*k3/k2
    a1 = 8.0*k1*k3/k2
    a2 = 2.0*k1/k2*( 1.0 + 8.0*k3*(n_o-n_c) ) + 2.0*k1*k3
    a3 = 8.0*k1/k2*(n_o-n_c) + 2.0*k3 + k1
    a4 = 8.0*k1/k2*(n_o-n_c)*(n_o-n_c) + 1.0 + 2.0*k1*(n_o-n_c)
    a5 = -2.0*n_c
    result = polynomial.polynomial.polyroots([a5,a4,a3,a2,a1,a0])
    n_ch4 = result[4]  # picks out correct root (computer-dependent)
    return n_ch4

# "old" function to compute mixing ratio for C2H2 (without nitrogen, from Heng & Lyons 2016)
def n_acet_old(n_o,n_c,temp,pbar):
    k1 = kk(temp,pbar)
    k2 = kk2(temp)
    k3 = kk3(temp,pbar)
    a0 = 8.0*k1*k3*k3/k2
    a1 = 8.0*k1*k3/k2
    a2 = 2.0*k1/k2*( 1.0 + 8.0*k3*(n_o-n_c) ) + 2.0*k1*k3
    a3 = 8.0*k1/k2*(n_o-n_c) + 2.0*k3 + k1
    a4 = 8.0*k1/k2*(n_o-n_c)*(n_o-n_c) + 1.0 + 2.0*k1*(n_o-n_c)
    a5 = -2.0*n_c
    result = polynomial.polynomial.polyroots([a5,a4,a3,a2,a1,a0])
    n_ch4 = result[4]  # picks out correct root (computer-dependent)
    n_c2h2 = k3*n_ch4*n_ch4
    return n_c2h2

# 9-molecule function to solve 10th order polynomial equation for CO mixing ratio
def n_cmono3(n_o,n_c,n_n,my_temperature,pbar):
    k1 = kk(my_temperature,pbar)
    k2 = kk2(my_temperature)
    k3 = kk3(my_temperature,pbar)
    k4 = kk4(my_temperature,pbar)
    k5 = kk5(my_temperature,pbar)
    k6 = kk6(my_temperature,pbar)
    d2 = 1.0 + 2.0*k1*(n_o+n_c)
    d3 = 1.0 + 1.0/k4
    c2 = 1.0/k2
    xx = c2*n_o
    f0 = 8.0*(k1**2)*(n_o**2)*n_c
    f1 = 2.0*k1*n_o*( -1.0 + 2.0*k1*( 2.0*n_c*(2.0*xx-1.0) - n_o ) - 4.0*k1*c2*(n_o**2) )
    f2 = k1*(1.0-8.0*xx) - 2.0*k3*d3 + 2.0*(k1**2)*( n_c*(1.0-8.0*xx) + 2.0*n_o*(1.0+xx) )
    f3 = 4.0*k1*c2*(1.0-2.0*xx) - 12.0*k3*d3*c2 + (k1**2)*( 2.0*c2*(2.0*n_c+n_o) - 1.0 )
    f4 = k1*c2*(4.0*c2-k1) - 24.0*k3*d3*(c2**2)
    f5 = -16.0*k3*d3*(c2**3)
    j0 = f0**2
    j1 = 2.0*f0*f1
    j2 = f1**2 + 2.0*f0*f2
    j3 = 2.0*f0*f3 + 2.0*f1*f2
    j4 = 2.0*f0*f4 + 2.0*f1*f3 + f2**2
    j5 = 2.0*f0*f5 + 2.0*f1*f4 + 2.0*f2*f3
    j6 = 2.0*f1*f5 + 2.0*f2*f4 + f3**2
    j7 = 2.0*f2*f5 + 2.0*f3*f4
    j8 = 2.0*f3*f5 + f4**2
    j9 = 2.0*f4*f5
    j10 = f5**2
    a0 = 2.0*k5*j0
    a1 = 2.0*k5*j1 + 2.0*k1*k6*f0*n_o
    a2 = 2.0*k5*j2 + k1*k6*( 2.0*n_o*f1 + f0*(8.0*xx-1.0) ) + k6*k6*f0 - 8.0*((k1*k6*n_o)**2)*n_n
    a3 = 2.0*k5*j3 + k1*k6*( 2.0*n_o*f2 + f1*(8.0*xx-1.0) + 4.0*f0*c2*(2.0*xx-1.0) ) + (k6**2)*( f1 + 6.0*f0*c2 ) + 8.0*((k1*k6)**2)*n_n*n_o*(1.0-8.0*xx)
    a4 = 2.0*k5*j4 + k1*k6*( 2.0*n_o*f3 + f2*(8.0*xx-1.0) + 4.0*f1*c2*(2.0*xx-1.0) - 4.0*f0*c2*c2 ) + (k6**2)*( f2 + 6.0*f1*c2 + 12.0*f0*(c2**2) ) - 2.0*((k1*k6)**2)*n_n*(1.0 - 32.0*xx + 96.0*(xx**2) )
    a5 = 2.0*k5*j5 + k1*k6*( 2.0*n_o*f4 + f3*(8.0*xx-1.0) + 4.0*f2*c2*(2.0*xx-1.0) - 4.0*f1*c2*c2 ) + (k6**2)*( f3 + 6.0*f2*c2 + 12.0*f1*(c2**2) + 8.0*f0*(c2**3) ) - 16.0*((k1*k6)**2)*n_n*c2*(1.0 - 12.0*xx + 16.0*(xx**2) )
    a6 = 2.0*k5*j6 + k1*k6*( 2.0*n_o*f5 + f4*(8.0*xx-1.0) + 4.0*f3*c2*(2.0*xx-1.0) - 4.0*f2*c2*c2 ) + (k6**2)*( f4 + 6.0*f3*c2 + 12.0*f2*(c2**2) + 8.0*f1*(c2**3) ) - 16.0*((k1*k6)**2)*n_n*(c2**2)*(3.0 - 16.0*xx + 8.0*(xx**2) )
    a7 = 2.0*k5*j7 + k1*k6*( f5*(8.0*xx-1.0) + 4.0*f4*c2*(2.0*xx-1.0) - 4.0*f3*c2*c2 ) + k6*k6*f5 + 2.0*(k6**2)*c2*( 3.0*f4 + 6.0*f3*c2 + 4.0*f2*(c2**2) ) - 64.0*((k1*k6)**2)*n_n*(c2**3)*(1.0 - 2.0*xx )
    a8 = 2.0*k5*j8 + 4.0*k1*k6*c2*( f5*(2.0*xx-1.0) - f4*c2 ) + 2.0*(k6**2)*c2*( 3.0*f5 + 6.0*f4*c2 + 4.0*f3*(c2**2) ) - 32.0*((k1*k6)**2)*n_n*(c2**4)
    a9 = 2.0*k5*j9 - 4.0*k1*k6*f5*c2*c2 + 4.0*((k6*c2)**2)*( 3.0*f5 + 2.0*f4*c2 )
    a10 = 2.0*k5*j10 + 8.0*k6*k6*f5*(c2**3)
    result = polynomial.polynomial.polyroots([a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10])
    result = result[result.real>0.0]
    result = result[result.real<2.0*n_o]
    if (len(result) > 1):
        real_result = result[0]
    else:
        real_result = result
    return real_result.real

# 6-molecule function to solve 5th order polynomial equation for CO mixing ratio
def n_cmono(n_o,n_c,n_n,my_temperature,pbar):
    k1 = kk(my_temperature,pbar)
    k5 = kk5(my_temperature,pbar)
    k6 = kk6(my_temperature,pbar)
    d2 = 1.0 + 2.0*k1*(n_o+n_c)
    a0 = 256.0*(k1**3)*k5*(n_o**3)*n_c*n_c
    a1 = 32.0*((k1*n_o)**2)*n_c*( k6 - 4.0*k5*( d2 + k1*n_c ) )
    a2 = 16.0*k1*k5*n_o*( 8.0*k1*k1*n_o*n_c + d2*d2 + 4.0*k1*d2*n_c ) + 8.0*k1*k6*n_o*( 2.0*k6*(n_c-n_n) - 2.0*k1*n_c - d2 )
    a3 = -8.0*k1*k5*( 4.0*k1*d2*n_o + 8.0*k1*k1*n_o*n_c + d2*d2 ) + 4.0*k1*k6*( 2.0*k1*n_o + d2 ) + 4.0*k6*k6*( 2.0*k1*n_n - d2 )
    a4 = 16.0*k1*k1*k5*( k1*n_o + d2 ) + 4.0*k1*k6*(k6-k1)
    a5 = -8.0*(k1**3)*k5
    result = polynomial.polynomial.polyroots([a0,a1,a2,a3,a4,a5])
    result = result[result.real>0.0]
    result = result[result.real<2.0*n_o]
    if (len(result) > 1):
        real_result = result[0]
    else:
        real_result = result
    return real_result.real


# 8-molecule function to solve 6th order polynomial equation for CO mixing ratio
def n_cmono2(n_o,n_c,n_n,my_temperature,pbar):
    k1 = kk(my_temperature,pbar)
    k3 = kk3(my_temperature,pbar)*(1.0 + 1.0/kk4(my_temperature,pbar))  # this sets K3 really to K3*D3 to include C2H4
    k5 = kk5(my_temperature,pbar)
    k6 = kk6(my_temperature,pbar)
    f0 = 8.0*k1*k1*n_c*n_o*n_o
    f1 = -2.0*k1*n_o*( 1.0 + 2.0*k1*(n_o + 2.0*n_c) )
    f2 = k1 - 2.0*k3 + 2.0*k1*k1*(n_c + 2.0*n_o)
    f3 = -k1*k1
    j0 = f0**2
    j1 = 2.0*f0*f1
    j2 = f1**2 + 2.0*f0*f2
    j3 = 2.0*f0*f3 + 2.0*f1*f2
    j4 = 2.0*f1*f3 + f2**2
    j5 = 2.0*f2*f3
    j6 = f3**2
    a0 = 2.0*k5*j0
    a1 = 2.0*k5*j1 + 2.0*k1*k6*n_o*f0
    a2 = 2.0*k5*j2 + k1*k6*(2.0*n_o*f1-f0) + k6*k6*f0 - 8.0*((k1*k6*n_o)**2)*n_n 
    a3 = 2.0*k5*j3 + k1*k6*(2.0*n_o*f2-f1) + k6*k6*f1 + 8.0*((k1*k6)**2)*n_n*n_o 
    a4 = 2.0*k5*j4 + k1*k6*(2.0*n_o*f3-f2) + k6*k6*f2 - 2.0*((k1*k6)**2)*n_n
    a5 = 2.0*k5*j5 + k6*f3*(k6-k1)
    a6 = 2.0*k5*j6
    result = polynomial.polynomial.polyroots([a0,a1,a2,a3,a4,a5,a6])
    result = result[result.real>0.0]
    result = result[result.real<2.0*n_o]
    if (len(result) > 1):
        real_result = result[0]
    else:
        real_result = result
    return real_result.real

# function to compute water mixing ratio
def n_water(n_o,n_c,n_n,my_temperature,pbar):
    c2 = 1.0/kk2(my_temperature)
    n_co = n_cmono(n_o,n_c,n_n,my_temperature,pbar)
    n_h2o = (2.0*n_o - n_co)/(1.0 + 2.0*c2*n_co)
    return n_h2o

# function to compute CH4 mixing ratio
def n_methane(n_o,n_c,n_n,my_temperature,pbar):
    k1 = kk(my_temperature,pbar)
    n_co = n_cmono(n_o,n_c,n_n,my_temperature,pbar)
    n_h2o = n_water(n_o,n_c,n_n,my_temperature,pbar)
    n_ch4 = n_co/k1/n_h2o
    return n_ch4

# function to compute NH3 mixing ratio
def n_ammonia(n_o,n_c,n_n,my_temperature,pbar):
    k5 = kk5(my_temperature,pbar)
    k6 = kk6(my_temperature,pbar)
    n_ch4 = n_methane(n_o,n_c,n_n,my_temperature,pbar)
    term1 = 1.0 + k6*n_ch4
    n_nh3 = ( sqrt( term1**2 + 16.0*n_n*k5 ) - term1 )/4.0/k5
    return n_nh3

# function to compute N2 mixing ratio
def n_nitrogen(n_o,n_c,n_n,my_temperature,pbar):
    k5 = kk5(my_temperature,pbar)
    n_nh3 = n_ammonia(n_o,n_c,n_n,my_temperature,pbar)
    result = k5*n_nh3*n_nh3
    return result

# function to compute HCN mixing ratio
def n_hcyan(n_o,n_c,n_n,my_temperature,pbar):
    k6 = kk6(my_temperature,pbar)
    n_nh3 = n_ammonia(n_o,n_c,n_n,my_temperature,pbar)
    n_ch4 = n_methane(n_o,n_c,n_n,my_temperature,pbar)
    result = k6*n_nh3*n_ch4
    return result

# function to compute CO2 mixing ratio
def n_cdio(n_o,n_c,n_n,my_temperature,pbar):
    k2 = kk2(my_temperature)
    n_co = n_cmono(n_o,n_c,n_n,my_temperature,pbar)
    n_h2o = n_water(n_o,n_c,n_n,my_temperature,pbar)
    result = n_co*n_h2o/k2
    return result

# function to compute C2H2 mixing ratio
def n_acet(n_o,n_c,n_n,my_temperature,pbar):
    k3 = kk3(my_temperature,pbar)
    n_ch4 = n_methane(n_o,n_c,n_n,my_temperature,pbar)
    result = k3*n_ch4*n_ch4
    return result

# function to compute C2H4 mixing ratio
def n_ethy(n_o,n_c,n_n,my_temperature,pbar):
    k4 = kk4(my_temperature,pbar)
    n_c2h2 = n_acet(n_o,n_c,n_n,my_temperature,pbar)
    result = n_c2h2/k4
    return result
    
# load TEA ouput
teafile = 'f'+plot_fig+'.tea' 
tea = np.loadtxt(teafile,skiprows=8)
c_o_tea = np.append(arange(0.1, 1., 0.1),arange(1,10.1,1.) )  

# Open the atmospheric file and read
f = open(teafile, 'r')
lines = np.asarray(f.readlines())
f.close()

# locate the molecule name line
imol = np.where(lines == "#SPECIES\n")[0][0] + 1
# creat the list of molecules
mol_list = lines[imol].split()

# read the column for a specific species: (need to use JANAF naming)
h2_col = mol_list.index('H2_ref')+2
ch4_col = mol_list.index('CH4_g')+2
h2o_col = mol_list.index('H2O_g')+2
co_col = mol_list.index('CO_g')+2
co2_col = mol_list.index('CO2_g')+2
c2h2_col = mol_list.index('C2H2_g')+2
c2h4_col = mol_list.index('C2H4_g')+2
hcn_col = mol_list.index('HCN_g')+2
n2_col = mol_list.index('N2_ref')+2
nh3_col = mol_list.index('NH3_g')+2

            
# MAIN #
pbar = 1.0e0 # pressure (bar)
n_o = 5e-4
n_c = np.logspace(-1.,1.,100) *n_o
CtoO = np.logspace(-1.,1.,100) # the C/O array
n_n = 0.2*n_o
if plot_fig == '2a':
    temperature = 800.
elif plot_fig == '2b':
    temperature = 1500.
else:
    print 'Please specify plot_fig ("2a","2b")'    

if plot_fig == '2a':
    title = r'$P=1$ bar, $T$ = 800 K, $\tilde{n}_{\rm O}=5\times 10^{-4}$, $\tilde{n}_{\rm N}=10^{-4}$ '
elif plot_fig == '2b':
    title = r'$P=1$ bar, $T$ = 1500 K, $\tilde{n}_{\rm O}=5\times 10^{-4}$, $\tilde{n}_{\rm N}=10^{-4}$ '

nt = len(n_c)

n_co = zeros(nt)
n_co_old = zeros(nt)
n_h2o = zeros(nt)
n_h2o_old = zeros(nt)
n_ch4 = zeros(nt)
n_ch4_old = zeros(nt)
n_nh3 = zeros(nt)
n_n2 = zeros(nt)
n_hcn = zeros(nt)
n_co2 = zeros(nt)
n_c2h2 = zeros(nt)
n_c2h2_old = zeros(nt)
n_c2h4 = zeros(nt)

for i in range(0,nt):
    n_co[i] = n_cmono(n_o,n_c[i],n_n,temperature,pbar)
    n_co_old[i] = n_cmono_old(n_o,n_c[i],temperature,pbar)
    n_h2o[i] = n_water(n_o,n_c[i],n_n,temperature,pbar)
    n_h2o_old[i] = n_water_old(n_o,n_c[i],temperature,pbar)
    n_ch4[i] = n_methane(n_o,n_c[i],n_n,temperature,pbar)
    n_ch4_old[i] = n_methane_old(n_o,n_c[i],temperature,pbar)
    n_nh3[i] = n_ammonia(n_o,n_c[i],n_n,temperature,pbar)
    n_n2[i] = n_nitrogen(n_o,n_c[i],n_n,temperature,pbar)
    n_hcn[i] = n_hcyan(n_o,n_c[i],n_n,temperature,pbar)
    n_co2[i] = n_cdio(n_o,n_c[i],n_n,temperature,pbar)
    n_c2h2[i] = n_acet(n_o,n_c[i],n_n,temperature,pbar)
    n_c2h2_old[i] = n_acet_old(n_o,n_c[i],temperature,pbar)
    n_c2h4[i] = n_ethy(n_o,n_c[i],n_n,temperature,pbar)

# PLOT
line1, =plt.plot(CtoO, n_ch4, linewidth=5, color='k', linestyle='--',zorder=0)
line2, =plt.plot(CtoO, n_h2o, linewidth=5, color='m', linestyle=':',zorder=0)
line3, =plt.plot(CtoO, n_co, linewidth=5, color='c', linestyle='-.',zorder=0)
line4, =plt.plot(CtoO, n_nh3, linewidth=5, color='b', linestyle='-',zorder=0)
line5, =plt.plot(CtoO, n_n2, linewidth=5, color='b', linestyle=':',zorder=0)
line6, =plt.plot(CtoO, n_hcn, linewidth=5, color='r', linestyle='-',zorder=0)
line7, =plt.plot(CtoO, n_co2, linewidth=5, color='g', linestyle='--',zorder=0)
line8, =plt.plot(CtoO, n_c2h2, linewidth=5, color='y', linestyle='--',zorder=0)
line9, =plt.plot(CtoO, n_c2h4, linewidth=5, color='y', linestyle=':',zorder=0)
plt.plot(CtoO, n_c2h2_old, linewidth=1, color='y', linestyle='-',zorder=0)
plt.plot(CtoO, n_co_old, linewidth=1, color='c', linestyle='-',zorder=0)
plt.plot(CtoO, n_ch4_old, linewidth=1, color='k', linestyle='-',zorder=0)
plt.plot(CtoO, n_h2o_old, linewidth=1, color='m', linestyle='-',zorder=0)
plt.scatter(c_o_tea, tea[:,ch4_col]/tea[:,h2_col], s=40, facecolor='1.', edgecolors='k',zorder=1)
plt.scatter(c_o_tea, tea[:,h2o_col]/tea[:,h2_col], s=40, facecolor='m', edgecolors='k',zorder=1)
plt.scatter(c_o_tea, tea[:,co_col]/tea[:,h2_col], s=40, facecolor='c', edgecolors='k',zorder=1)
plt.scatter(c_o_tea, tea[:,nh3_col]/tea[:,h2_col], s=40, facecolor='b', edgecolors='k',zorder=1)
plt.scatter(c_o_tea, tea[:,n2_col]/tea[:,h2_col], s=40, facecolor='b', edgecolors='k',zorder=1)
plt.scatter(c_o_tea, tea[:,hcn_col]/tea[:,h2_col], s=40, facecolor='r', edgecolors='k',zorder=1)
plt.scatter(c_o_tea, tea[:,co2_col]/tea[:,h2_col], s=40, facecolor='g', edgecolors='k',zorder=1)
plt.scatter(c_o_tea, tea[:,c2h2_col]/tea[:,h2_col], s=40, facecolor='y', edgecolors='k',zorder=1)
plt.scatter(c_o_tea, tea[:,c2h4_col]/tea[:,h2_col], s=40, facecolor='y', edgecolors='k',zorder=1)

plt.yscale('log')
plt.xscale('log')
plt.xlim([0.1,10.])
plt.ylim([1e-20,1e0])
plt.title(title,fontsize=12)
plt.text(1000, 1e-15, 'thin curves: Heng & Lyons (2016)', fontsize=12)
plt.text(1000, 1e-16, 'solutions without nitrogen', fontsize=12)
plt.text(1000, 1e-18, 'circles: TEA Gibbs free energy minimization code', fontsize=12)
plt.xlabel('C/O', fontsize=18)
plt.ylabel(r'$\tilde{n}_{\rm X}$', fontsize=14)
plt.legend([line1,line2,line3,line4,line5,line6,line7,line8,line9],[r'CH$_4$', r'H$_2$O', r'CO', r'NH$_3$', r'N$_2$', r'HCN', r'CO$_2$', r'C$_2$H$_2$', r'C$_2$H$_4$'],frameon=False,prop={'size':12},loc=4)
plt.savefig('heng_tsai_'+plot_fig+'.eps', format='eps') #save in EPS format
