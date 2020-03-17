---
title: Quick install guide
theme: _config.yml
filename: quick-guides
---

# Quick Guides
## Writing Input
ADPRES input is designed to be self-explanatory. It has several input cards, for example: `%mode`, `%geom`, `%xsec`, and so on. Some cards are mandatory for any problems. While some cards are conditional, depending on the problem being solved and some cards are optional. Comments are marked by `!`. For example, the following is the [IAEA3D input](https://github.com/imronuke/ADPRES/tree/master/smpl/static), where you can find its specification [here](https://engineering.purdue.edu/PARCS/Code/TestSuite/CalculationMode/StandAloneMode/Eigenvalue/IAEA3DPWR).

```
! IAEA3D input data
! NODE SIZE = 10 cm
! PARCS K-EFF  : 1.029096
! ADPRES K-EFF : 1.029082 (ERROR = 1.4 PCM)

! Mode card
%MODE
FORWARD

! Case card
%CASE
IAEA3D
10 CM NODE SIZE

! Cross-sections card
%XSEC
2  5    ! Number of groups and number of materials
! sigtr   siga   nu*sigf sigf   chi   sigs_g1  sigs_g2
0.222222  0.010  0.000  0.000    1.0   0.1922   0.020
0.833333  0.080  0.135  0.135    0.0   0.000    0.7533   ! MAT1 : Outer Fuel
0.222222  0.010  0.000  0.000    1.0   0.1922   0.020
0.833333  0.085  0.135  0.135    0.0   0.000    0.7483   ! MAT2 : Inner Fuel
0.222222  0.0100 0.000  0.000    1.0   0.1922   0.020
0.833333  0.1300 0.135  0.135    0.0   0.000    0.7033   ! MAT3 : Inner Fuel + Control Rod
0.166667  0.000  0.000  0.000    0.0   0.1267   0.040
1.111111  0.010  0.000  0.000    0.0   0.000    1.1011   ! MAT4 : Reflector
0.166667  0.000  0.000  0.000    0.0   0.000    0.040
1.111111  0.055  0.000  0.000    0.0   0.000    0.000    ! MAT5 : Reflector + Control Rod
%GEOM
9 9 19         ! number of assembly in x, y, z directions
10.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0    !x-direction assembly size in cm
1      8    8    8    8    8    8    8    8     !x-direction assembly divided into 2 (10 cm each)
20.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0 10.0    !y-direction assembly size in cm
8      8    8    8    8    8    8    8    1     !y-direction assembly divided into 2 (10 cm each)
19*20.0                                         !z-direction assembly  in cm
19*1                                            !z-direction nodal is not divided
4                                               !np number of planar type
1  13*2  4*3  4                                 !planar assignment (from bottom to top)
! Planar_type_1 (Bottom Reflector)
  4  4  4  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  0
  4  4  4  4  4  4  4  4  0
  4  4  4  4  4  4  4  0  0
  4  4  4  4  4  4  0  0  0
  4  4  4  4  0  0  0  0  0
! Planar_type_2 (Fuel)
  3  2  2  2  3  2  2  1  4
  2  2  2  2  2  2  2  1  4
  2  2  2  2  2  2  1  1  4
  2  2  2  2  2  2  1  4  4
  3  2  2  2  3  1  1  4  0
  2  2  2  2  1  1  4  4  0
  2  2  1  1  1  4  4  0  0
  1  1  1  4  4  4  0  0  0
  4  4  4  4  0  0  0  0  0
! Planar_type_3 (Fuel+Partial Control Rods)
  3  2  2  2  3  2  2  1  4
  2  2  2  2  2  2  2  1  4
  2  2  3  2  2  2  1  1  4
  2  2  2  2  2  2  1  4  4
  3  2  2  2  3  1  1  4  0
  2  2  2  2  1  1  4  4  0
  2  2  1  1  1  4  4  0  0
  1  1  1  4  4  4  0  0  0
  4  4  4  4  0  0  0  0  0
! Planar_type_4 (Top reflectors)
  5  4  4  4  5  4  4  4  4
  4  4  4  4  4  4  4  4  4
  4  4  5  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  4
  5  4  4  4  5  4  4  4  0
  4  4  4  4  4  4  4  4  0
  4  4  4  4  4  4  4  0  0
  4  4  4  4  4  4  0  0  0
  4  4  4  4  0  0  0  0  0
! Boundary conditions
! 0 = zero-flux
! 1 = zero-incoming current
! 2 = reflective
(east), (west), (north), (south), (bottom), (top)
   1       2       2        1        1        1

! NOTE: Writing 19*20.0 is equivalent to write 20.0 nineteen times in a row
```

In the above example, there are
1. Two mandatory cards : `%MODE` and `%GEOM`
2. One conditional card: `%XSEC`
3. One optional card   : `%CASE`

### Card `%MODE`
This is the mode of ADPRES calculation. Since here we want to calculate static forward calculation (eigenvalue problem) the calculation mode is `FORWARD`

### Card `%CASE`
This card is optional. This describes the problem at hand.

### Card `%XSEC`
This card is conditional, needed only if `XTAB` card is not present. This card tells ADPRES the cross sections data for the problem. The cross section data must be given for each group and for each material as shown in the example. The description of the cross sections data can be seen in the comments.

### Card `%GEOM`
This card is describes the geometry of the problem. It quite similar to other reactor core simulator which you can easily understand if you have background on nuclear engineering. The description of the inputs given in the comments.


## Reading Output 

```
           ###########################################################
           #                     ADPRES 1.2                          #
           #        ABU DHABI POLYTECHNIC REACTOR SIMULATOR          #
           ###########################################################


  CALCULATION MODE : FORWARD CALCULATION                                         

  CASE ID : IAEA3D                                                                                              
  10 CM NODE SIZE                                                                                     

  NODAL KERNEL  : SEMI-ANALYTIC NODAL METHOD

  reading input ... done


  ==============================================================================
                           CALCULATION RESULTS
  ==============================================================================

  Itr     k-eff     Fis.Src Error   Inner Error
 ----------------------------------------------------
    1     0.981424    5.47871E-01    8.55259E+03
    2     1.001319    2.41976E-01    8.56379E+00
    3     1.009804    1.67297E-01    6.27992E-01
    4     1.013500    1.33833E-01    1.44559E-01
     ...FISSION SOURCE EXTRAPOLATED...
    5     1.026980    2.65156E+00    1.31348E-01
    6     1.024050    2.69425E-01    1.90647E+00
    7     1.024688    9.31087E-02    4.43462E+01
    8     1.024358    3.99972E-02    1.10592E+00
    9     1.024191    3.47063E-02    1.96477E-01
     ...FISSION SOURCE EXTRAPOLATED...
   10     1.023838    1.12907E+00    5.87407E-02
   11     1.024188    1.12316E-01    9.62168E-01
   12     1.025744    1.03113E-01    4.24920E-01
   13     1.026299    8.09335E-02    1.86532E-01
   14     1.026714    6.40064E-02    8.00816E-02
     ...FISSION SOURCE EXTRAPOLATED...
   15     1.028565    5.05638E-01    6.73943E-02
   16     1.028253    6.31685E-02    4.20531E-01
   17     1.028196    3.63469E-02    4.93252E-01
   18     1.028082    2.21127E-02    2.48688E-01
   19     1.028020    1.29756E-02    6.24972E-02
     ...FISSION SOURCE EXTRAPOLATED...
   20     1.027400    3.54750E-01    2.34432E-02
   21     1.027655    6.04390E-02    3.09183E-01
     .....NODAL COUPLING UPDATED..... 
MAX. CHANGE IN NODAL COUPLING COEF.=  3.16843E-01 AT NODE I =  6, J =  4, K = 19
   22     1.028084    5.88420E-02    1.59061E-01
   23     1.032373    1.30872E-01    4.33336E-01
   24     1.030818    5.51189E-02    1.84604E-01
     ...FISSION SOURCE EXTRAPOLATED...
   25     1.029029    1.07105E-01    1.14532E-01
   26     1.029615    1.28000E-02    7.92164E-02
   27     1.029310    9.21108E-03    7.02347E-02
   28     1.029268    7.70113E-03    2.33171E-02
   29     1.029252    6.77810E-03    1.26452E-02
     ...FISSION SOURCE EXTRAPOLATED...
   30     1.028885    8.85937E-02    8.34156E-03
   31     1.029019    1.01734E-02    8.43263E-02
   32     1.028826    5.42157E-03    8.78422E-02
   33     1.028802    4.30053E-03    2.23840E-02
   34     1.028806    3.96582E-03    4.21410E-03
     ...FISSION SOURCE EXTRAPOLATED...
   35     1.028882    4.00875E-02    4.15830E-03
   36     1.028869    3.16753E-03    3.66234E-02
   37     1.028896    2.89270E-03    2.66432E-02
   38     1.028890    2.29630E-03    7.77722E-03
   39     1.028884    1.75597E-03    4.50688E-03
     ...FISSION SOURCE EXTRAPOLATED...
   40     1.028806    2.11309E-02    2.01927E-03
   41     1.028850    5.04331E-03    1.98071E-02
   42     1.028822    2.16511E-03    2.31885E-02
   43     1.028816    1.32030E-03    7.40873E-03
     .....NODAL COUPLING UPDATED..... 
MAX. CHANGE IN NODAL COUPLING COEF.=  7.77068E-02 AT NODE I = 13, J =  6, K = 19
   44     1.028815    9.88078E-04    1.75021E-03
     ...FISSION SOURCE EXTRAPOLATED...
   45     1.028749    2.70297E-03    1.01361E-01
   46     1.029184    1.79029E-02    2.54126E-02
   47     1.029130    7.71244E-03    1.02339E-02
   48     1.029084    5.16741E-03    7.30180E-03
   49     1.029046    3.73690E-03    4.77112E-03
     ...FISSION SOURCE EXTRAPOLATED...
   50     1.028894    2.60208E-02    3.45389E-03
   51     1.028998    2.72863E-03    2.28899E-02
   52     1.028973    1.53306E-03    1.53902E-02
   53     1.028980    1.13497E-03    1.55657E-03
   54     1.028984    1.08660E-03    1.08660E-03
     ...FISSION SOURCE EXTRAPOLATED...
   55     1.029008    2.33384E-02    1.05436E-03
   56     1.029009    2.73691E-03    2.15387E-02
   57     1.029003    1.84511E-03    9.83932E-03
   58     1.028990    1.28141E-03    2.51189E-03
   59     1.028982    9.58612E-04    1.31166E-03
     ...FISSION SOURCE EXTRAPOLATED...
   60     1.028957    3.62385E-03    7.81147E-04
   61     1.028967    3.17276E-04    2.97034E-03
   62     1.028964    3.10507E-04    2.33945E-03
   63     1.028966    2.64908E-04    3.56434E-04
   64     1.028966    2.34934E-04    2.64160E-04
     ...FISSION SOURCE EXTRAPOLATED...
   65     1.028972    7.96762E-03    2.25512E-04
     .....NODAL COUPLING UPDATED..... 
MAX. CHANGE IN NODAL COUPLING COEF.=  9.59384E-03 AT NODE I = 13, J =  7, K = 19
   66     1.028973    9.35948E-04    7.33196E-03
   67     1.029140    3.29672E-03    3.07920E-02
   68     1.029115    1.48019E-03    4.09482E-03
   69     1.029102    1.18427E-03    2.15202E-03
     ...FISSION SOURCE EXTRAPOLATED...
   70     1.029051    5.24286E-03    1.08493E-03
   71     1.029076    4.68687E-04    4.47135E-03
   72     1.029064    2.97197E-04    4.10696E-03
   73     1.029065    2.40066E-04    7.22833E-04
   74     1.029065    2.25154E-04    5.88166E-04
     ...FISSION SOURCE EXTRAPOLATED...
   75     1.029067    4.52712E-03    2.25768E-04
   76     1.029066    5.41916E-04    4.39262E-03
   77     1.029063    4.49023E-04    4.09320E-03
   78     1.029059    2.86348E-04    1.02676E-03
   79     1.029057    2.05863E-04    2.77957E-04
     ...FISSION SOURCE EXTRAPOLATED...
   80     1.029050    8.85328E-04    1.75320E-04
   81     1.029052    1.19252E-04    7.58061E-04
   82     1.029052    9.36462E-05    5.70339E-04
   83     1.029053    8.00869E-05    1.55104E-04
   84     1.029053    7.01886E-05    7.92477E-05
     ...FISSION SOURCE EXTRAPOLATED...
   85     1.029055    1.39428E-03    6.88354E-05
   86     1.029055    1.88081E-04    1.35060E-03
   87     1.029053    1.08592E-04    1.20189E-03
     .....NODAL COUPLING UPDATED..... 
MAX. CHANGE IN NODAL COUPLING COEF.=  7.06979E-03 AT NODE I = 10, J =  6, K = 19
   88     1.029053    8.09990E-05    2.29796E-04
   89     1.029097    8.74304E-04    7.62786E-03
     ...FISSION SOURCE EXTRAPOLATED...
   90     1.029088    1.00450E-03    1.44001E-03
   91     1.029089    4.13243E-04    8.96275E-04
   92     1.029086    1.82495E-04    8.18932E-04
   93     1.029085    1.44250E-04    3.07788E-04
   94     1.029084    1.25895E-04    1.51294E-04
     ...FISSION SOURCE EXTRAPOLATED...
   95     1.029076    1.68868E-03    1.24464E-04
   96     1.029081    1.94833E-04    1.56586E-03
   97     1.029078    8.03913E-05    1.47759E-03
   98     1.029077    7.00997E-05    1.50746E-04
   99     1.029077    5.89739E-05    8.19896E-05
     ...FISSION SOURCE EXTRAPOLATED...
  100     1.029078    6.02650E-04    1.23325E-04
  101     1.029078    5.60653E-05    5.44002E-04
  102     1.029078    4.29217E-05    5.72034E-04
  103     1.029078    2.94079E-05    7.34477E-05
  104     1.029077    2.58479E-05    3.93366E-05
     ...FISSION SOURCE EXTRAPOLATED...
  105     1.029076    2.66503E-04    2.45627E-05
  106     1.029077    2.64558E-05    2.42296E-04
  107     1.029076    2.70479E-05    2.06213E-04
  108     1.029076    1.91791E-05    3.50766E-05
  109     1.029076    1.56900E-05    2.79154E-05
     ...FISSION SOURCE EXTRAPOLATED...
     .....NODAL COUPLING UPDATED..... 
MAX. CHANGE IN NODAL COUPLING COEF.=  2.95708E-03 AT NODE I = 12, J =  8, K = 19
  110     1.029076    1.39353E-04    2.50027E-05
  111     1.029088    2.72689E-04    3.43600E-03
  112     1.029087    1.33611E-04    4.30031E-04
  113     1.029086    7.84636E-05    2.05477E-04
  114     1.029086    6.42833E-05    8.26350E-05
     ...FISSION SOURCE EXTRAPOLATED...
  115     1.029082    5.28830E-04    6.02213E-05
  116     1.029084    6.59343E-05    4.72776E-04
  117     1.029083    5.33407E-05    4.38936E-04
  118     1.029083    3.60340E-05    9.12013E-05
  119     1.029083    2.87736E-05    8.59804E-05
     ...FISSION SOURCE EXTRAPOLATED...
  120     1.029084    3.82801E-04    4.74233E-05
  121     1.029084    6.93807E-05    4.04142E-04
  122     1.029083    4.11835E-05    4.28494E-04
  123     1.029083    2.61380E-05    9.62027E-05
  124     1.029083    2.04651E-05    3.91481E-05
     ...FISSION SOURCE EXTRAPOLATED...
  125     1.029082    1.08870E-04    2.56609E-05
  126     1.029083    1.24363E-05    9.51748E-05
  127     1.029082    1.20168E-05    8.28493E-05
  128     1.029082    9.48481E-06    2.04727E-05
  129     1.029082    7.46358E-06    9.66768E-06

  MULTIPLICATION EFFECTIVE (K-EFF) =  1.029082


  CPU time breakdown in seconds
    Input reading time   :    0.0078  ( 5.0%)
    XSEC processing time :    0.0003  ( 0.2%)
    CMFD time            :    0.0893  (57.1%)
    Nodal update time    :    0.0590  (37.7%)
    T-H time             :    0.0000  ( 0.0%)
    ------------------------------------------
    Total time           :    0.1565
  ```
