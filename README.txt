-----------------------------------------------------------------------
04-26-2019

ALADYN_PI mini-app is a simple molecular dynamics (MD) code for performing 
constant energy atomistic simulations using either a straight 
Artificial Neural Network (ANN) interatomic potential, or
Physically Informed Neural Network (PINN) interatomic potential. 

References:
 1. Vesselin I. Yamakov, Edward H. Glaessgen, NASA/TM-2019-220431
    (https://ntrs.nasa.gov/citations/20200000069)
 2. G.P.Purja Pun, V.Yamakov, J.Hickman, E.H.Glaessgen, Y.Mishin, 
    Physical Review Materials, 4, 113807 (2020).

The trained ANN was produced and offered by Ganga P. Pun and Yuri Mishin 
from George Mason University.

=======================================================================

 COMPILATION

-----------------------------------------------------------------------
* FORTRAN compiler (Intel 2003 or newer, or PGI FORTRAN with OpenACC)
-----------------------------------------------------------------------
Source directory: ALADYN_PI.source
-----------------------------------------------------------------------
Compilation: use the provided makefiles with the following options:

Intel compiler:
> make -f makefile.intel            ! compiles with -O3 optimization on 
> make -f makefile.intel DEBUG=TRUE ! check and warning flags on
> make -f makefile.intel OMP=TRUE   ! compile with OpenMP direectives

PGI compiler:
> make -f makefile.pgi              ! compiles with -O3 optimization on 
> make -f makefile.pgi DEBUG=TRUE   ! check and warning flags on
> make -f makefile.pgi OMP=TRUE     ! compile with OpenMP direectives
> make -f makefile.pgi ACC=TRUE  ! compile with OpenMP+OpenACC direectives

Edit the provided makefiles for specific compiler option of your choice
=======================================================================

 EXECUTION

-----------------------------------------------------------------------
Run from the example (test) directory: ALADYN_PI.test

Running a tets case:
> aladyn_pi              ! executes 10 MD steps report at each step
> aladyn_pi -n 100 -m 10 ! executes 100 MD steps, report at each 10-th step

Available PBS scripts used for NASA/LaRC K3 cluster:
ALADYN_PI_intel.job      - run the intel version
ALADYN_PI_pgi_V100.job   - run the pgi version with OpneACC for V100 gpu.

Screen output is riderected to aladyn_pi.out

Example outputs are saved in: aladyn_pi_intel.out and aladyn_pi_pgi.out
=======================================================================

INPUT FILES:

-----------------------------------------------------------------------
PINN.dat - Neural network potential file
structure.plt - input atomic structure file

-----------------------------------------------------------------------
--- Available potential files in directory POT:
(to be linked or copied as PINN.dat in the working directory of aladyn_pi)

ANN_Si.dat - Straight Artificial Neural Network potential for Si
PINN_Si.dat - Physically Informed Neural Network potential for Si
PINN_Al.dat - Physically Informed Neural Network potential for Al (ref.2)

-----------------------------------------------------------------------
--- Available Si single crystal test structures in directory STR_Si:

Si_N4000.plt      -   4000 atoms Si crystal
Si_N8000.plt      -   8000 atoms Si crystal
Si_N16000.plt     -  16000 atoms Si crystal
Si_N32000.plt     -  32000 atoms Si crystal
Si_N64000.plt     -  64000 atoms Si crystal
Si_N128000.plt    - 128000 atoms Si crystal
Si_N192000.plt    - 192000 atoms Si crystal
Si_N256000.plt    - 256000 atoms Si crystal
Si_N512000.plt    - 512000 atoms Si crystal

-----------------------------------------------------------------------
--- Available Al single crystal test structures in directory STR_Al:

Al_N4000.plt      -   4000 atoms Al crystal
Al_N8000.plt      -   8000 atoms Al crystal
Al_N16000.plt     -  16000 atoms Al crystal
Al_N32000.plt     -  32000 atoms Al crystal
Al_N64000.plt     -  64000 atoms Al crystal
Al_N128000.plt    - 128000 atoms Al crystal
Al_N192000.plt    - 192000 atoms Al crystal
Al_N256000.plt    - 256000 atoms Al crystal

Use any of the above structures by linking them to structure.plt, e.g.,
> ln -s STR/Al_N4000.plt structure.plt

=======================================================================

SOURCE FILES:

    aladyn_pi.f       - Main program
    aladyn_pi_sys.f   - system modul
    aladyn_pi_sys_OMP.f    - system modul for OpneMP compilation
    aladyn_pi_sys_NO_OMP.f - system modul without OpneMP compilation
    aladyn_pi_sys_ACC.f    - system modul for OpneACC compilation
    aladyn_pi_mods.f  - contains general purpose modules
    aladyn_pi_IO.f    - I/O operations

    aladyn_pi_ANN_OMP.f  - Artificial Neural Network OpneMP code
    contains:
     subroutine Frc_ANN_OMP  ! ANN force & energy (OpenMP version)

    aladyn_pi_ANN_ACC.f  - Artificial Neural Network OpneACC code
     subroutine Frc_ANN_ACC  ! ANN force & energy (OpenACC version)

    aladyn_pi_PINN_OMP.f - Physically Informed NN OpneMP code
     subroutine Frc_PINN_OMP ! PINN force & energy (OpenMP version)

    aladyn_pi_PINN_ACC.f - Physically Informed NN OpneACC code
     subroutine Frc_PINN_ACC ! PINN force & energy (OpenACC version)

    aladyn_pi_MD.f    - molecular dynamics module
     contains:
      subroutine get_T     ! Calculates current system temperature
      subroutine predict_atoms ! Gear predictor call !
      subroutine correct_atoms ! Gear corrector call !

-----------------------------------------------------------------------
Suggested subroutines for optimization: 
Frc_ANN_OMP and Frc_ANN_ACC (Optimized versions from ALADYN miniapp)

Frc_PINN_OMP and Frc_PINN_ACC
(Optimized versions from ALADYN miniapp)

-----------------------------------------------------------------------
 For further information contact:

 Vesselin Yamakov
 National Institute of Aerospace
 100 Exploration Way,
 Hampton, VA 23666
 phone: (757)-864-2850
 fax:   (757)-864-8911
 e-mail: yamakov@nianet.org

=======================================================================
 Notices:
 Copyright 2020 United States Government as represented by the 
 Administrator of the National Aeronautics and Space Administration. 
 All Rights Reserved.
 
 Disclaimers:
 No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY 
 WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, 
 INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE 
 WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF 
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM 
INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR 
 FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM 
 TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, 
 CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT 
 OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY 
 OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  
 FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES 
 REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, 
 AND DISTRIBUTES IT "AS IS." 

 Waiver and Indemnity:  
 RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES 
 GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR 
 RECIPIENT. IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY 
 LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH 
 USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, 
 RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND 
 HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND 
 SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED 
 BY LAW. RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE 
 IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
=======================================================================

