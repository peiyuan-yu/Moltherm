
Running Job 1 of 1 /scratch1/scratchdirs/ewcss/with_water/243_2674637/rct_1_471171.in
qchem /scratch1/scratchdirs/ewcss/with_water/243_2674637/rct_1_471171.in_43027.0 /tmp/qchem43027/ 1
srun -n 1 /global/common/edison/software/qchem/5.1beta/exe/qcprog.exe /scratch1/scratchdirs/ewcss/with_water/243_2674637/rct_1_471171.in_43027.0 /tmp/qchem43027/

Process 0 of 1 is on nid04204 - thread support 0
initial socket setup ...start
initial socket setup ...done 
now start server 0 ... 
Check license passed

                  Welcome to Q-Chem
     A Quantum Leap Into The Future Of Chemistry


 Q-Chem 5.1, Q-Chem, Inc., Pleasanton, CA (2018)

 Yihan Shao,  Zhengting Gan,  E. Epifanovsky,  A. T. B. Gilbert,  M. Wormit,  
 J. Kussmann,  A. W. Lange,  A. Behn,  Jia Deng,  Xintian Feng,  D. Ghosh,  
 M. Goldey,  P. R. Horn,  L. D. Jacobson,  I. Kaliman,  T. Kus,  A. Landau,  
 Jie Liu,  E. I. Proynov,  R. M. Richard,  R. P. Steele,  E. J. Sundstrom,  
 H. L. Woodcock III,  P. M. Zimmerman,  D. Zuev,  B. Albrecht,  E. Alguire,  
 S. A. Baeppler,  D. Barton,  Z. Benda,  Y. A. Bernard,  E. J. Berquist,  
 K. B. Bravaya,  H. Burton,  D. Casanova,  Chun-Min Chang,  Yunqing Chen,  
 A. Chien,  K. D. Closser,  M. P. Coons,  S. Coriani,  S. Dasgupta,  
 A. L. Dempwolff,  M. Diedenhofen,  Hainam Do,  R. G. Edgar,  Po-Tung Fang,  
 S. Faraji,  S. Fatehi,  Qingguo Feng,  K. D. Fenk,  J. Fosso-Tande,  
 Qinghui Ge,  A. Ghysels,  G. Gidofalvi,  J. Gomes,  J. Gonthier,  A. Gunina,  
 D. Hait,  M. W. D. Hanson-Heine,  P. H. P. Harbach,  A. W. Hauser,  
 J. E. Herr,  E. G. Hohenstein,  Z. C. Holden,  Kerwin Hui,  B. C. Huynh,  
 T.-C. Jagau,  Hyunjun Ji,  B. Kaduk,  K. Khistyaev,  Jaehoon Kim,  
 P. Klunzinger,  K. Koh,  D. Kosenkov,  L. Koulias,  T. Kowalczyk,  
 C. M. Krauter,  A. Kunitsa,  Ka Un Lao,  A. Laurent,  K. V. Lawler,  
 Joonho Lee,  D. Lefrancois,  S. Lehtola,  D. S. Levine,  Yi-Pei Li,  
 You-Sheng Lin,  Fenglai Liu,  E. Livshits,  A. Luenser,  P. Manohar,  
 E. Mansoor,  S. F. Manzer,  Shan-Ping Mao,  Yuezhi Mao,  N. Mardirossian,  
 A. V. Marenich,  T. Markovich,  L. A. Martinez-Martinez,  S. A. Maurer,  
 N. J. Mayhall,  S. C. McKenzie,  J.-M. Mewes,  A. F. Morrison,  
 J. W. Mullinax,  K. Nanda,  T. S. Nguyen-Beck,  R. Olivares-Amaya,  
 J. A. Parkhill,  T. M. Perrine,  F. Plasser,  P. Pokhilko,  S. Prager,  
 A. Prociuk,  E. Ramos,  D. R. Rehn,  F. Rob,  M. Schneider,  N. Sergueev,  
 S. M. Sharada,  S. Sharma,  D. W. Small,  T. Stauch,  T. Stein,  
 Yu-Chuan Su,  A. J. W. Thom,  A. Tkatchenko,  T. Tsuchimochi,  N. M. Tubman,  
 L. Vogt,  M. L. Vidal,  O. Vydrov,  M. A. Watson,  J. Wenzel,  
 M. de Wergifosse,  T. A. Wesolowski,  A. White,  J. Witte,  A. Yamada,  
 Jun Yang,  K. Yao,  S. Yeganeh,  S. R. Yost,  Zhi-Qiang You,  A. Zech,  
 Igor Ying Zhang,  Xing Zhang,  Yan Zhao,  Ying Zhu,  B. R. Brooks,  
 G. K. L. Chan,  C. J. Cramer,  M. S. Gordon,  W. J. Hehre,  A. Klamt,  
 M. W. Schmidt,  C. D. Sherrill,  D. G. Truhlar,  A. Aspuru-Guzik,  R. Baer,  
 A. T. Bell,  N. A. Besley,  Jeng-Da Chai,  R. A. DiStasio Jr.,  A. Dreuw,  
 B. D. Dunietz,  T. R. Furlani,  Chao-Ping Hsu,  Yousung Jung,  Jing Kong,  
 D. S. Lambrecht,  WanZhen Liang,  C. Ochsenfeld,  V. A. Rassolov,  
 L. V. Slipchenko,  J. E. Subotnik,  T. Van Voorhis,  J. M. Herbert,  
 A. I. Krylov,  P. M. W. Gill,  M. Head-Gordon

 Contributors to earlier versions of Q-Chem not listed above: 
 R. D. Adamson,  B. Austin,  J. Baker,  G. J. O. Beran,  K. Brandhorst,  
 S. T. Brown,  E. F. C. Byrd,  A. K. Chakraborty,  C.-L. Cheng,  
 Siu Hung Chien,  D. M. Chipman,  D. L. Crittenden,  H. Dachsel,  
 R. J. Doerksen,  A. D. Dutoi,  L. Fusti-Molnar,  W. A. Goddard III,  
 A. Golubeva-Zadorozhnaya,  S. R. Gwaltney,  G. Hawkins,  A. Heyden,  
 S. Hirata,  G. Kedziora,  F. J. Keil,  C. Kelley,  Jihan Kim,  R. A. King,  
 R. Z. Khaliullin,  P. P. Korambath,  W. Kurlancheek,  A. M. Lee,  M. S. Lee,  
 S. V. Levchenko,  Ching Yeh Lin,  D. Liotard,  R. C. Lochan,  I. Lotan,  
 P. E. Maslen,  N. Nair,  D. P. O'Neill,  D. Neuhauser,  E. Neuscamman,  
 C. M. Oana,  R. Olson,  B. Peters,  R. Peverati,  P. A. Pieniazek,  
 Y. M. Rhee,  J. Ritchie,  M. A. Rohrdanz,  E. Rosta,  N. J. Russ,  
 H. F. Schaefer III,  N. E. Schultz,  N. Shenvi,  A. C. Simmonett,  A. Sodt,  
 D. Stuck,  K. S. Thanthiriwatte,  V. Vanovschi,  Tao Wang,  A. Warshel,  
 C. F. Williams,  Q. Wu,  X. Xu,  W. Zhang

 Please cite Q-Chem as follows:
 Y. Shao et al., Mol. Phys. 113, 184-215 (2015)
 DOI: 10.1080/00268976.2014.952696

 Q-Chem 5.1.0 for Intel X86 EM64T Linux

 Parts of Q-Chem use Armadillo 8.300.2 (Tropical Shenanigans).
 http://arma.sourceforge.net/

 Q-Chem begins on Wed Jul 18 14:58:37 2018  

Host: edison07
0

     Scratch files written to /tmp/qchem43027//
 Apr1918  _RNUM -1
 Parallel job on  1  processors
Processing $rem in /global/common/edison/software/qchem/5.1beta/config/preferences.
 MEM_TOTAL  2000 
Symmetry turned off for PCM/SM12/SMD calculation

Checking the input file for inconsistencies... 	...done.

--------------------------------------------------------------
User input:
--------------------------------------------------------------
$molecule
0 1
C      1.2086456594      0.0000023618      0.0000000000
C      0.2747861808     -0.0000035811     -1.1771999381
C      0.2747861808     -0.0000035811      1.1771999381
C     -0.9941447273     -0.0000007204     -0.7338254040
C     -0.9941447273     -0.0000007204      0.7338254040
H      1.8682158342     -0.8790894329      0.0000000000
H      1.8681921457      0.8791127858      0.0000000000
H      0.6023067980     -0.0000100561     -2.2099834487
H      0.6023067980     -0.0000100561      2.2099834487
H     -1.8888166634     -0.0000007743     -1.3458131664
H     -1.8888166634     -0.0000007743      1.3458131664
$end

$rem
job_type = sp
method = wb97x-d
basis = 6-311++g(d,p)
max_scf_cycles = 200
gen_scfman = True
scf_algorithm = diis
solvent_method = pcm
$end

$pcm
theory iefpcm
$end

$solvent
dielectric 80.4
$end
--------------------------------------------------------------
 ----------------------------------------------------------------
             Standard Nuclear Orientation (Angstroms)
    I     Atom           X                Y                Z
 ----------------------------------------------------------------
    1      C       1.2086456594     0.0000023618     0.0000000000
    2      C       0.2747861808    -0.0000035811    -1.1771999381
    3      C       0.2747861808    -0.0000035811     1.1771999381
    4      C      -0.9941447273    -0.0000007204    -0.7338254040
    5      C      -0.9941447273    -0.0000007204     0.7338254040
    6      H       1.8682158342    -0.8790894329     0.0000000000
    7      H       1.8681921457     0.8791127858     0.0000000000
    8      H       0.6023067980    -0.0000100561    -2.2099834487
    9      H       0.6023067980    -0.0000100561     2.2099834487
   10      H      -1.8888166634    -0.0000007743    -1.3458131664
   11      H      -1.8888166634    -0.0000007743     1.3458131664
 ----------------------------------------------------------------
 Nuclear Repulsion Energy =   156.8067309507 hartrees
 There are       18 alpha and       18 beta electrons
 Requested basis set is 6-311++G(d,p)
 There are 60 shells and 152 basis functions

 Total QAlloc Memory Limit   2000 MB
 Mega-Array Size       188 MB
 MEM_STATIC part       192 MB
 Discretize the solute cavity surface with Lebedev spheres
	Using 302 Lebedev grid points for each H atom
	Using 302 Lebedev grid points for other atoms
	Atomic van der Waals radii will be scaled by 1.20
 Remove points where switching function is < 1.0e-08
 Keep 1445 surface tesserae and discard 1877 interior tesserae
 Molecular Surface Area = 113.495 Angst**2

                       Distance Matrix (Angstroms)
             C (  1)   C (  2)   C (  3)   C (  4)   C (  5)   H (  6)
   C (  2)  1.502629
   C (  3)  1.502629  2.354400
   C (  4)  2.321806  1.344160  2.293949
   C (  5)  2.321806  2.293949  1.344160  1.467651
   H (  6)  1.099016  2.167397  2.167397  3.082921  3.082921
   H (  7)  1.099016  2.167392  2.167392  3.082907  3.082907  1.758202
   H (  8)  2.291653  1.083472  3.402981  2.174327  3.348831  2.694315
   H (  9)  2.291653  3.402981  1.083472  3.348831  2.174327  2.694315
   H ( 10)  3.377201  2.170163  3.323669  1.083959  2.263920  4.086478
   H ( 11)  3.377201  3.323669  2.170163  2.263920  1.083959  4.086478
             H (  7)   H (  8)   H (  9)   H ( 10)
   H (  8)  2.694318
   H (  9)  2.694318  4.419967
   H ( 10)  4.086462  2.636757  4.341588
   H ( 11)  4.086462  4.341588  2.636757  2.691626
 
 A cutoff of  1.0D-08 yielded   1747 shell pairs
 There are     11501 function pairs (     12301 Cartesian)
 
 -------------------------------------------------------
 OpenMP Integral Computing Module                       
 Release: version 1.0, May 2013, Q-Chem Inc. Pittsburgh 
 -------------------------------------------------------
 Integral Job Info:
 Integral job number is                      11
 Integral operator is                         1
 short-range coefficients              22203600
 long-range coefficients              100000000
 Omega coefficients                         200
 if combine SR and LR in K                    1
 Integral screening is                        0
 Integral computing path is                   2
 max size of driver memory is            800000
 size of driver memory is                593474
 size of scratch memory is              2428936
 max col of scratch BK array               1296
 max len of scratch array in speh3          155
 max len of scratch index in speh4           18
 max int batch size is                      520
 min int batch size is                       52
 fixed nKL is                                52
 max L of basis functions is                  2
 order of int derivative is                   0
 number of shells is                         60
 number of basis is                         157
 number of cartesian basis is               157
 number of contracted shell pairs          1747
 number of primitive shell pairs           4242
 maxK2 (contraction) of shell pair           36
 max number of K2 of shell pair               1
 max number of CS2 of shell pair            288
 max number of PS2 of shell pair            540
 mem total for path MDJ                  520922
 -------------------------------------------------------
 Smallest overlap matrix eigenvalue = 7.73E-06

 Q-Chem warning in module stvman/mkSTV.C, line 295:

 Overlap eigenvalue is smaller than square root of threshold


 Scale SEOQF with 1.000000e-01/1.000000e-01/1.000000e-01

 Standard Electronic Orientation quadrupole field applied
 Nucleus-field energy     =    -0.0000000025 hartrees
 Guess from superposition of atomic densities
 Warning:  Energy on first SCF cycle will be non-variational
 SAD guess density has 36.000000 electrons

 -----------------------------------------------------------------------
  General SCF calculation program by
  Eric Jon Sundstrom, Paul Horn, Yuezhi Mao, Dmitri Zuev, Alec White,
  David Stuck, Shaama M.S., Shane Yost, Joonho Lee, David Small,
  Daniel Levine, Susi Lehtola, Hugh Burton, Evgeny Epifanovsky
 -----------------------------------------------------------------------
 Exchange:     0.2220 Hartree-Fock + 1.0000 wB97X-D + LR-HF
 Correlation:  1.0000 wB97X-D
 Using SG-2 standard quadrature grid
 Dispersion:   Grimme D
 using 24 threads for integral computing
 -------------------------------------------------------
 OpenMP Integral computing Module                
 Release: version 1.0, May 2013, Q-Chem Inc. Pittsburgh 
 -------------------------------------------------------
 using 24 threads for integral computing
 -------------------------------------------------------
 OpenMP Integral computing Module                
 Release: version 1.0, May 2013, Q-Chem Inc. Pittsburgh 
 -------------------------------------------------------
 -------------------------------------------------------
 OpenMP BLAS3 based DFT computing Module                
 Release: version 1.0, May 2013, Q-Chem Inc. Pittsburgh 
 -------------------------------------------------------
 A restricted SCF calculation will be
 performed using DIIS
 SCF converges when DIIS error is below 1.0e-05
 ---------------------------------------
  Cycle       Energy         DIIS error
 ---------------------------------------
 Using non-symmetric K matrix
 Warning: not using a symmetric Q
    1    -195.5406243012      5.08e-02  
    2    -194.0347250769      3.93e-03  
    3    -194.0615833427      2.88e-03  
    4    -194.0850185100      2.14e-04  
    5    -194.0852417243      5.80e-05  
    6    -194.0852583028      2.40e-05  
    7    -194.0852607200      1.01e-05  
    8    -194.0852613094      2.54e-06  Convergence criterion met
 ---------------------------------------
 SCF time:   CPU 151.81s  wall 9.00s 
 SCF   energy in the final basis set =     -194.0852613094

************** Final PCM Free Energy Summary **************
 G_electrostatic  =      -0.00466663 hartree =      -2.92835510 kcal/mol
 G_cavitation     =       0.00000000 hartree =       0.00000000 kcal/mol
 G_dispersion     =       0.00000000 hartree =       0.00000000 kcal/mol
 G_repulsion      =       0.00000000 hartree =       0.00000000 kcal/mol
 --------------------------------------------------
 Non-electrostatic Free Energy =       0.00000000 hartree =       0.00000000 kcal/mol
 Total                         =      -0.00466663 hartree =      -2.92835510 kcal/mol
 --------------------------------------------------
 SCF Energy (H0 + V/2)                       =    -194.08526131 
 Solute Internal Energy (H0)                 =    -194.08059468 
 Total Free Energy (H0 + V/2 + non-elec)     =    -194.08526131 hartree
                                             = -121790.33965317 kcal/mol
***********************************************************

 
 --------------------------------------------------------------
 
                    Orbital Energies (a.u.)
 --------------------------------------------------------------
 
 Alpha MOs
 -- Occupied --
-10.286 -10.282 -10.282 -10.281 -10.281  -0.962  -0.806  -0.801
 -0.643  -0.621  -0.605  -0.524  -0.483  -0.459  -0.441  -0.440
 -0.386  -0.304
 -- Virtual --
  0.044   0.046   0.057   0.059   0.075   0.078   0.084   0.124
  0.125   0.136   0.148   0.152   0.159   0.178   0.178   0.188
  0.190   0.191   0.205   0.207   0.210   0.212   0.215   0.216
  0.219   0.264   0.271   0.279   0.293   0.302   0.303   0.323
  0.365   0.388   0.390   0.406   0.464   0.485   0.529   0.558
  0.565   0.595   0.617   0.639   0.639   0.640   0.647   0.692
  0.695   0.708   0.735   0.758   0.761   0.823   0.862   0.864
  0.873   0.898   0.910   0.924   0.943   0.975   1.015   1.067
  1.085   1.104   1.109   1.173   1.213   1.303   1.487   1.510
  1.558   1.646   1.648   1.653   1.659   1.660   1.664   1.732
  1.753   1.801   1.824   1.901   1.912   1.942   1.944   1.979
  2.011   2.032   2.034   2.066   2.116   2.199   2.251   2.326
  2.364   2.409   2.416   2.513   2.544   2.545   2.649   2.678
  2.704   2.716   2.763   2.812   2.825   2.847   2.850   2.873
  2.923   2.938   3.071   3.123   3.129   3.158   3.284   3.441
  3.475   3.583   3.614   3.867   4.145   4.152   4.180   4.246
  4.687  23.807  23.952  24.039  24.365  24.396
 --------------------------------------------------------------
 
          Ground-State Mulliken Net Atomic Charges

     Atom                 Charge (a.u.)
  ----------------------------------------
      1 C                    -0.758612
      2 C                    -0.002176
      3 C                    -0.002176
      4 C                    -0.288955
      5 C                    -0.288955
      6 H                     0.191464
      7 H                     0.191466
      8 H                     0.244070
      9 H                     0.244070
     10 H                     0.234901
     11 H                     0.234901
  ----------------------------------------
  Sum of atomic charges =    -0.000000

 -----------------------------------------------------------------
                    Cartesian Multipole Moments
 -----------------------------------------------------------------
    Charge (ESU x 10^10)
                -0.0000
    Dipole Moment (Debye)
         X       0.7705      Y       0.0000      Z      -0.0000
       Tot       0.7705
    Quadrupole Moments (Debye-Ang)
        XX     -26.2834     XY       0.0000     YY     -33.4364
        XZ       0.0000     YZ       0.0000     ZZ     -28.3047
    Octopole Moments (Debye-Ang^2)
       XXX      -0.9269    XXY       0.0000    XYY       3.7355
       YYY       0.0002    XXZ      -0.0000    XYZ       0.0000
       YYZ      -0.0000    XZZ       0.8952    YZZ      -0.0001
       ZZZ      -0.0000
    Hexadecapole Moments (Debye-Ang^3)
      XXXX    -191.0844   XXXY      -0.0002   XXYY     -41.2305
      XYYY       0.0000   YYYY     -47.3797   XXXZ       0.0000
      XXYZ       0.0000   XYYZ       0.0000   YYYZ       0.0000
      XXZZ     -63.0221   XYZZ      -0.0001   YYZZ     -49.9439
      XZZZ      -0.0000   YZZZ       0.0000   ZZZZ    -190.4592
 -----------------------------------------------------------------
Archival summary:
1\1\nid04204\SP\ProcedureUnspecified\BasisUnspecified\56\ewcss\WedJul1814:58:462018WedJul1814:58:462018\0\\#,ProcedureUnspecified,BasisUnspecified,\\0,1\C\H,1,1.08396\C,1,1.34416,2,126.367\H,3,1.08347,1,126.855,2,-0,0\C,1,1.46765,2,124.374,3,-180,0\H,5,1.08396,1,124.374,2,0,0\C,5,1.34416,1,109.26,2,180,0\H,7,1.08347,5,126.855,1,-180,0\C,7,1.50263,5,109.165,1,0,0\H,9,1.09902,7,111.899,5,-120.447,0\H,9,1.09902,7,111.899,5,120.445,0\\\@

 Total job time:  8.91s(wall), 161.07s(cpu) 
 Wed Jul 18 14:58:46 2018

        *************************************************************
        *                                                           *
        *  Thank you very much for using Q-Chem.  Have a nice day.  *
        *                                                           *
        *************************************************************


0 sent ACK to 0 
now end server 0 ... 
cleanup process ... done
