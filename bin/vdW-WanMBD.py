######################################################################################################################
##################################################### vdW-WanMBD #####################################################
###################################### Diem TX Dang, Dai-Nam Le and Lilia Woods ######################################
########################################## University of South Florida, USA ##########################################
######################################################################################################################

import datetime
start = datetime.datetime.now().replace( microsecond = 0)                       # calculate time of process 

# Import numpy, math
import numpy as np
import math, cmath                                                              # cmath: math for complex number
import os
from scipy.linalg import logm                                                   # logm: logarithm of matrix
from itertools import islice

# Constants 
hbar = 1.05457182E-34                                                           # Dirac's constant (J.s)
e = 1.602176634E-19                                                             # electron charge (C)
epsilon_0 = 8.854187817E-22                                                     # convert 8.854187817E-12 F/m to F/Å
r_au = 0.529177210903                                                           # Bohr radius (Å)
Ry = 13.6057039763                                                              # Convert Rydberg Constant to eV
fac1 = ( 2 * math.pi)**( -1)                                                    # factor to calculate vdW energy
fac2 = -e * 4**( -3) * ( math.pi * epsilon_0)**( -2)                            # factor to calculate induction energy
       
### IMPORT DATA

# Import input.in file
y1 = 'seedname ='
y2 = 'N_WF_at ='
y3 = 'gM ='
y4 = 'pol ='
y5 = 'mono ='
with open( 'input.in', 'r') as f:
    for line in f:
        if y1 in line:
            seedname = line[ line.find( y1) + len( y1) : 34].strip()
        if y2 in line:
            N_WF_at = np.array( line[ line.find( y2) + len( y2) : 34].split()) 
        if y3 in line:
            gM = float( line[ line.find( y3) + len( y3) : 34]) 
        if y4 in line:
            pol = float( line[ line.find( y4) + len( y4) : 34]) 
        if y5 in line:
            mono = int( line[ line.find( y5) + len( y5) : 34])             

# Import number of atoms per cell
y6 = 'number of atoms/cell      ='
y7 = 'number of atomic types    ='
y8 = '!    total energy              ='
with open( 'scf.out', 'r') as f:
    for line in f:
        if y6 in line:
            N_at = int( line[ line.find( y6) + len( y6) : 46])
        if y7 in line:
            N_ty = int( line[ line.find( y7) + len( y7) : 46]) 
        if y8 in line:
            E_tot = float( line[ line.find( y8) + len( y8) : 50])           

# Number of atoms in monolayer
v = int( N_at / 2) 
fac3 = 1000 / N_at                                                              # factor to convert to meV/atom unit

# Import lattice vectors (Å)    
y9 = 'Number of Wannier Functions               :'
y10 = 'Lattice Vectors (Ang)'
y11 = 'Cartesian Coordinate (Ang)'
y12 = 'Final State'
y13 = 'Unit Cell Volume:'
with open( seedname + '.wout', 'r') as f:
    for line in f:
        if y9 in line:
            N_WF = int( line[ line.find( y9) + len( y9) : 66])                  # number of Wannier functions
        if y10 in line: 
            latt1 = "".join( islice( f, 3)) 
        if y11 in line: 
            coor1 = "".join( islice( f, N_at + 1))   
        if y12 in line: 
            S1 = "".join( islice( f, N_WF))   
        if y13 in line:
            V_bi = float( line[ line.find( y13) + len( y13) : 51])              # volume of bilayer cell (Å^3)                

list1 = [ "a_1", "a_2", "a_3"]
for i in list1:
    latt1 = latt1.replace( i, '')
    
latt = np.fromstring( latt1, sep = '    ').reshape( 3, 3)                       # lattice vectors

# Import coordinates 
list2 = [ "+", "-", "|"]
for i in list2:
    coor1 = coor1.replace( i, '')

at1 = ''.join([ i for i in coor1 if not i.isdigit()])
list3 = [ ".", " "]
for i in list3:
    at1 = at1.replace( i, '')

at2 = at1.splitlines()
while( '' in at2):
    at2.remove( '')

at = []                                                                         # types of atoms
for i in at2:
    if i not in at:
        at.append( i) 

for i in at:
    coor1 = coor1.replace( i, '')
    
coor2 = np.fromstring( coor1, sep = '    ').reshape( N_at, 7)  
coor3 = np.delete( coor2, slice( 4), 1)                                         # coordinates of atoms

# Import number of Wannier functions of each atom            
N_WF_s = N_WF**2                                                                # square of number of Wannier functions
coor4 = np.delete( coor2, slice( 1, 7), 1)
N_at_ty = np.zeros(( N_ty), np.dtype( np.int32))                                # number of atoms corresponding to each type 
N_at_ty[ N_ty - 1] = coor4[ N_at - 1] 
pos = np.where( coor4 == 1)[ 0]
for i in range ( N_ty - 1):
    N_at_ty[ i] = coor4[ pos[ i + 1] - 1]  
   
N_WF_mt1 = np.repeat( N_WF_at, N_at_ty, axis = 0)                               
N_WF_mt = np.zeros(( N_at, 1), dtype = ( np.int32))
N_WF_mt[ :, 0] = N_WF_mt1[ :]                                                   # number of Wannier functions for each atom

# Import the largest spread of each atom     
list4 = ["WF centre and spread", "(", ",", ")"]
for i in list4:
    S1 = S1.replace( i, '')
    
S2 = np.fromstring( S1, sep = '    ').reshape( N_WF, 5)
S3 = np.delete( S2, slice( 4), 1) 
S4 = np.sqrt( S3)
S_lar = np.zeros(( N_at, 1))
for i in range( N_at):
    S_lar[ i, 0] = S4[ i * N_WF_mt[ i, 0] : ( i + 1) * N_WF_mt[ i, 0], :].max( axis = 0)  

# Import Loewdin charge
char = np.zeros(( N_at, 1))
if pol == 1:
    y11 = 'Loewdin charge'
    with open( 'CHARGE.lobster', 'r') as f:
        for line in f:
            if y11 in line:
                char1 = "".join( islice( f, N_at)) 

    for i in at:
        char1 = char1.replace( i, '')

    char2 = np.fromstring( char1, sep = '    ').reshape( N_at, 3)  
    char = np.delete( char2, slice( 2), 1)

# Combine coordinates, number of Wannier functions of each atom, largest spread of each atom, and Loewdin charge into an array
coor5 = np.concatenate(( coor3, N_WF_mt, S_lar, char), axis = 1)    

# Import momentum vectors k
y14 = 'number of k points='   
y15 = 'cart. coord. in units 2pi/alat'   
y16 = 'lattice parameter (alat)  =       '     
with open( 'nscf.out', 'r') as f:
    for line in f:
        if y14 in line:
            N_k = int( line[ line.find( y14) + len( y14) : 31])                 # number of kpoints
        if y15 in line: 
            kpt1 = "".join( islice( f, int( N_k))) 
        if y16 in line:
            a = float( line[ line.find( y16) + len( y16) : 45])                 # lattice parameter     
            
list5 = [ "k(", ") = (", "), wk ="]
for i in list5:
    kpt1 = kpt1.replace( i, '')
 
kpt2 = np.fromstring( kpt1, sep ='    ').reshape( N_k, 5) 
k_vec = np.delete( kpt2, 0, 1)  
k_vec[ :, :3] *= 2 * math.pi / ( a * r_au)                                      # cartesian coordinates of kpoints in Å^(-1)
 
# Import matrix element of r  
with open( seedname + '_r.dat', 'r') as fin:
    data = fin.read().splitlines( True)
with open( seedname + '_r_2.dat', 'w') as fout:
    fout.writelines( data[ 4:])
r_data = np.loadtxt( seedname + '_r_2.dat') 
os.remove( seedname + '_r_2.dat')               

# Import optical polarizability at imag. freq. (Siemens/Å)
xx = np.loadtxt( seedname + '-polarizability_xx.dat')           
xy = np.loadtxt( seedname + '-polarizability_xy.dat')          
xz = np.loadtxt( seedname + '-polarizability_xz.dat')          
yy = np.loadtxt( seedname + '-polarizability_yy.dat')          
yz = np.loadtxt( seedname + '-polarizability_yz.dat')           
zz = np.loadtxt( seedname + '-polarizability_zz.dat')  

# Import total energy of monolayer
with open( 'mono/scf.out', 'r') as f:
    for line in f:
        if y8 in line:
            E_mono = float( line[ line.find( y8) + len( y8) : 50])    

### CREATE ARRAYS, VARIABLES
N_L = int( len( r_data) / N_WF_s)                                               # number of translational vectors
N_ome = len( xx)                                                                # number of frequencies
E = xx[ :, 0]                                                                   # energy (eV)
ome = E * e / hbar                                                              # omega (s^-1)
A = np.zeros(( N_ome, 3 * N_at, 3 * N_at))                                      # optical polarizability at imag freq.        
L = np.zeros(( N_L, 6))                                                         # translational vectors array
r = np.zeros(( N_at, N_at, 3))                                                  # distances between atoms 
rL = np.zeros(( N_at, N_at, N_L, 3))                                            # R_ij,L
kL = np.zeros(( N_k, N_L))                                                      # vector k cdot vecor L
R = np.zeros(( N_at, N_at, N_L))                                                # R_ij,L magnetude
T = np.zeros(( 3 * N_at, 3 * N_at, N_L))
f = np.zeros(( N_at, N_at, N_L))
T_k = np.zeros(( N_k, 3 * N_at, 3 * N_at), dtype = ( np.complex128))
T1 = np.zeros(( N_at, N_at, 3, 3))
H = H2 = H4 = H6 = sum_S_cube = Q = 0      

### CALCULATION                    
# Translational vectors
for l in range( N_L):                                                           # index of vector L
    for i in range( 3):                                                         # x, y, z
        L[ l, i] = r_data[ l * N_WF_s, i]
        L[ l, i + 3] = r_data[ l * N_WF_s, 0] * latt[ 0, i] + r_data[ l * N_WF_s, 1] * latt[ 1, i] + r_data[ l * N_WF_s, 2] * latt[ 2, i]

# Sort atoms in ascending order by z coordinate        
coor = coor5[ coor5[ :, 2].argsort()]                                           # columns 0, 1, 2: coordinates; column 3: number of WFs of each atom, column 4: largest spread, column 5: Loewdin charge 

# Sum of cube of spreads
for b in range( N_at):
    sum_S_cube += coor[ b, 3] * coor[ b, 4]**3

# Matrix A (Eq (7))
for b in range( N_at):                                                                                                                                 # index of atoms
    [[ A[ :, 3 * b + 0, 3 * b + 0], A[ :, 3 * b + 0, 3 * b + 1], A[ :, 3 * b + 0, 3 * b + 2] ],
    [ A[ :, 3 * b + 1, 3 * b + 0], A[ :, 3 * b + 1, 3 * b + 1], A[ :, 3 * b + 1, 3 * b + 2] ],
    [ A[ :, 3 * b + 2, 3 * b + 0], A[ :, 3 * b + 2, 3 * b + 1], A[ :, 3 * b + 2, 3 * b + 2] ]] = np.asarray([[ xx[ :, 1], xy[ :, 1], xz[ :, 1] ], [ xy[ :, 1], yy[ :, 1], yz[ :, 1] ],
    [ xz[ :, 1], yz[ :, 1], zz[ :, 1] ]]) * V_bi * coor[ b, 3] * coor[ b, 4]**3 / sum_S_cube                                                           

# Distance of 2 atoms in same cell
for b1 in range( v):
    for b2 in range( v, N_at):
        for i in range( 3):
            r[ b1, b2, i] = coor[ b1, i] - coor[ b2, i]   
             
# Distance between 2 atoms shifted from original unit cell by a lattice translation    
for b1 in range( v):                                                                                                                                   # atom in monolayer 1
    for b2 in range( v, N_at):                                                                                                                         # atom in monolayer 2
        for l in range( N_L):                                                                                                                          # index of vector L
            for i in range( 3):                                                                                                                        # x, y, z 
                rL[ b1, b2, l, i] = r[ b1, b2, i] + L[ l, i + 3]                              
            
# Interaction tensor, matrix T (Eq(9))
for k in range( N_k):                                                                                                                                  # index of vector k                            
    for b1 in range( v):                                                                                                                               # atom in monolayer 1
        for b2 in range( v, N_at):                                                                                                                     # atom in monolayer 2
            for l in range( N_L):                                                                                                                      # index of vector L
                kL[ k, l] = k_vec[ k, 0] * L[ l, 3] + k_vec[ k, 1] * L[ l, 4] + k_vec[ k, 2] * L[ l, 5]                                                # dot product of vector k and vector L
                R[ b1, b2, l] = np.sqrt( rL[ b1, b2, l, 0]**2 + rL[ b1, b2, l, 1]**2 + rL[ b1, b2, l, 2]**2)                                           # magnitude of distance of 2 atoms 
                f[ b1, b2, l] = ( 1 + math.exp( -1 * ( R[ b1, b2, l] / ( gM * ( coor[ b1, 4] + coor[ b2, 4])) - 1)))**( -1)                            # Damping function Eq (10)                    
                for i1 in range( 3):                                                                                                                   # x, y, z
                    for i2 in range( 3):                                                                                                               # x, y, z                    
                        if i2 == i1:
                            T[ 3 * b1 + i1, 3* b2 + i2, l] = ( 3 * rL[ b1, b2, l, i1] * rL[ b1, b2, l, i2] - R[ b1, b2, l]**2) / R[ b1, b2, l]**5      # interaction tensor Eq (4)                      
                        else:
                            T[ 3 * b1 + i1, 3* b2 + i2, l] = 3 * rL[ b1, b2, l, i1] * rL[ b1, b2, l, i2] / R[ b1, b2, l]**5                            # interaction tensor Eq (4)  
                        T_k[ k, 3 * b1 + i1, 3 * b2 + i2] += f[ b1, b2, l] * T[ 3 * b1 + i1, 3* b2 + i2, l] * cmath.exp( -1j * kL[ k, l])                                                         
                        T_k[ k, 3 * b2 + i2, 3 * b1 + i1] += f[ b1, b2, l] * T[ 3 * b1 + i1, 3* b2 + i2, l] * cmath.exp( +1j * kL[ k, l])              # T_k is a hermitian matrices                                                        
                                      
# Integral Eqs (5)-(6)
for k in range( N_k):                                                           # index of vector k              
    for m in range( 1, N_ome):                                                  # index of omega
        A_LR = A[ m, :, :]                                                      # A_LR matrix at given omega
        T_LR = T_k[ k, :, :] * ( 4 * math.pi * epsilon_0)**( -1)                # T_LR matrix at given vector k
        AT = np.dot( A_LR, T_LR)                                                # AT muplitication
        B = logm( np.identity( 3 * N_at) - np.dot( A_LR, T_LR))                 # take natural log of matrix
        B2 = -(1/2) * np.linalg.matrix_power( AT, 2)                            # 2nd order expansion of log of matrix
        B4 = -(1/4) * np.linalg.matrix_power( AT, 4)                            # 4th order expansion of log of matrix
        B6 = -(1/6) * np.linalg.matrix_power( AT, 6)                            # 6th order expansion of log of matrix
        H  += k_vec[ k, 3] * ( E[ m] - E[ m - 1]) * B.trace()                   # sum on omega in eV and k with w_k for many-body
        H2 += k_vec[ k, 3] * ( E[ m] - E[ m - 1]) * B2.trace()                  # sum on omega in eV and k with w_k for 2nd order
        H4 += k_vec[ k, 3] * ( E[ m] - E[ m - 1]) * B4.trace()                  # sum on omega in eV and k with w_k for 4th order
        H6 += k_vec[ k, 3] * ( E[ m] - E[ m - 1]) * B6.trace()                  # sum on omega in eV and k with w_k for 6th order
        
Emb = fac1 * H.real                                                             # many-body dipole-dipole interaction (eV) Eq(6)
E2  = fac1 * H2.real                                                            # 2nd order dipole-dipole interaction (eV) n = 2 Eq (5)
E4  = fac1 * H4.real                                                            # 4th order dipole-dipole interaction (eV) n = 4 Eq (5)
E6  = fac1 * H6.real                                                            # 6th order dipole-dipole interaction (eV) n = 6 Eq (5)

# Induction energy (Eq (11))
E_ind = 0
if pol == 1:
    for b1 in range( v):                                                                                                                               # atom in monolayer 1
        for b2 in range( v, N_at):                                                                                                                     # atom in monolayer 2
            for i1 in range( 3):                                                                                                                       # x, y, z
                for i2 in range( 3):                                                                                                                   # x, y, z              
                    for l in range( N_L):                                                                                                              # index of vector L                     
                        T1[ b1, b2, i1, i2] += f[ b1, b2, l]**2 * rL[ b1, b2, l, i1] * rL[ b1, b2, l, i2] / R[ b1, b2, l]**6                                                                
                    Q += T1[ b1, b2, i1, i2] * ( coor[ b1, 5]**2 * A[ 0, 3 * b2 + i1, 3 * b2 + i2] + coor[ b2, 5]**2 * A[ 0, 3 * b1 + i1, 3 * b1 + i2])  
    E_ind = fac2 * Q                     

### EXTRACT OUTPUT

print( '#####################################################')
print( '#################### vdW-WanMBD #####################')
print( '###### Diem TX Dang, Dai-Nam Le and Lilia Woods #####')
print( '######### University of South Florida, USA ##########')
print( '#####################################################')
print( '')
print( 'The investigating system:', seedname)
print( 'Lattice constant along c axis (Angstrom):', latt[ 2, 2])
print( 'Interlayer distance (Angstrom)          :', latt[ 2, 2] / 2)
print( 'Unit cell volume (Angstrom^3)           :', V_bi)
print( '')
print( '#####################################################')
print( '############## Dispersion energy (eV) ###############')
print( '#####################################################')
print( '1. Dipole-dipole interaction')
print( 'Many-body vdW energy:', Emb)                                                 
print( '2nd order vdW energy:', E2)
print( '4th order vdW energy:', E4)
print( '6th order vdW energy:', E6)
print( '2. Charge-dipole interaction')
print( 'Induction energy    :', E_ind)  
print( '')

if mono == 1:
    E_b_DFT = ( E_tot - 2 * E_mono) * Ry                                        # convert from Ry to eV unit
    print( '#####################################################')
    print( '#### Binding energy of layered system (meV/atom) ####')
    print( '#####################################################')
    print( 'At DFT level       :', E_b_DFT * fac3)
    E_b = E_b_DFT + Emb + E_ind                                                 # Eq (1)   
    print( 'At vdW-WanMBD level:', E_b * fac3)

print( '')    
end = datetime.datetime.now().replace( microsecond = 0)
print( '#####################################################')
print( '######################## END ########################')
print( '################# Duration:', end - start,'#################')
print( '#####################################################')

