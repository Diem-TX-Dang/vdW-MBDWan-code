# vdW-MBDWan-code
A Python code for calculating van der Waals and induction energies in materials using Wannier Functions
----------------------------------------------------------------------------------------------------------------------------------

This code is provided for the paper: "Dissecting van der Waals interactions with Density Functional Theory - Wannier-basis approach" 

                                     https://arxiv.org/abs/...

Authors: Diem Thi-Xuan Dang, Dai-Nam Le, and Lilia Woods

Address: Advanced Materials and Devices Theory Group, Department of Physics, University of South Florida, Tampa, FL, USA

         https://www.amd-woods-group.com/
	 
------------------------------------------------------------------------------------------------------------------------------------

DEPENDENCIES:

python3 with Numpy, scipy.

------------------------------------------------------------------------------------------------------------------------------------

Usage Notes:
To compute the van der Waals and/or induction energies in materials, the following steps must be executed:

- Obtain the electronic structure of a given system using the Quantum ESPRESSO code, without a van der Waals correction.

- Obtain the MLWFs for the system by interpolating the DFT application using wannier90.x

- Determine the optical polarizability of the system in imaginary frequency by copying the executive postw90_vdW-WanMBD.x to the current folder and excuting: ./postw90_vdW-WanMBD.x seedname. This will generate the following output files:'seedname'-polarizability_xx.dat, 'seedname'-polarizability_xy.dat, ..., 'seedname'-polarizability_zz.dat.    

- Create a folder to execute vdW-WanMBD code and copy the following files to this folder

  	all 9 optical polarizability files 'seedname-polarizability_xx.dat', ..., 'seedname-polarizability_zz.dat'
  
	'scf.out' and 'nscf.out' files from the DFT self-consistent calculations for the bulk layered system
  
	'seedname.wout', 'seedname_r.dat' files for the bulk layered system calculated by Wannier90
  
	'CHARGE.lobster' file from partial charge population analysis by LOBSTER code for the bulk layered system  (this file is needed if the investigating system is polarized i.e. pol = 1)

  	'scf.out' file from the DFT self-consistent calculation for the monolayer component into subfolder /mono (this file is needed if calculating binding energy i.e. mono = 1)
  
	'vdW-WanMBD.py' file from 'bin' folder
		
- In this folder, provide the following information into the input file 'input.in':

  	name of the investigating system: seedname = ...                      Example: seedname = MoS2
  
	numbers of Wannier functions for each type of atoms: N_WF_at = ...    Example: N_WF_at = 5 4 for the case of MoS2 (5 Wannier functions for Mo and 4 Wannier Functions for S)

	gamma, the parameter in damping function (Eq. 10): gM = ...           Example: gM = 1 or gM = 1.31 in this study

	for polar systems, such as hBN, MoS2: pol = 1. For non-polar systems, such as graphite: pol = 0.

  	if calculating the binding energy: mono = 1. Otherwise: mono = 0
	
- Execute: python vdW-WanMBD.py > output.out

- The generated output file 'output.out' includes the following information:
  
  the name of investigated system

  lattice constant along c axis in Å unit
  
  interlayer distance in Å unit
      
  unit cell volume in Å^3 unit
  
  vdW energy in eV unit

  induction energy in eV unit (if 'pol = 0', induction energy is zero)

  binding energy in meV/atom unit (if 'mono = 0', this information is not publised)

---------------------------------------------------------------------------------------------------------------------------------------

Examples:

Dispersion and binding energies are included for three representative materials.
     
- examples/graphite: calculation for AB-stacking graphite at interlayer distance 3.4162060 Å
     
- examples/hBN     : calculation for AA'-stacking hBN     at interlayer distance 3.3086120 Å
     
- examples/MoS2    : calculation for AA'-stacking 2H-MoS2 at interlayer distance 5.9922775 Å
    
----------------------------------------------------------------------------------------------------------------------------------------

For assistance with running please email Diem Thi-Xuan Dang at dangt1@usf.edu 

----------------------------------------------------------------------------------------------------------------------------------------

vdW-WANMBD code is released under the GNU General Public License ver. 3. Please consult the included LICENSE file for detailed licensing conditions.


