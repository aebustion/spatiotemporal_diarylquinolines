;; 1. Based on: run59.mod (2023_04_19)
;; x1. Author: aebustion
$PROBLEM    rabbit lesion spatiotemporal modeling
$INPUT      ID OLDID=DROP WT DAT2=DROP TIME TAD EVID AMT DVN BLQ DIST
            CMPD DOSE_FACTOR iCL iCLm DV NDOSES HW
$DATA      TBAJ876_restructured_ids_distance.csv IGNORE=@
            IGNORE(CMPD.EQ.2) IGNORE(BLQ.EQ.1)
$SUBROUTINE ADVAN13 TOL=6
$MODEL      NCOMP=8 COMP=(ABS,DEFDOSE) COMP=(CENTRAL,DEFOBS)
            COMP=(PERIPH) COMP(CMPD) COMP(TRANSIT1) COMP(TRANSIT2)
            COMP=(PLESION) COMP=(MLESION)
$PK  
IF(CMPD.EQ.1) L2 = 1
IF(CMPD.EQ.2) L2 = 2

TVCL = 37.26
CL=TVCL*EXP(ETA(1))
V  = THETA(3)
KE = CL/V
KA = THETA(4)
S2 = V
V2 = THETA(5)
Q = THETA(6)
K12 = Q/V
K21 = Q/V2
Fm = THETA(7)
KM   = THETA(8)
VMAX = KM*CL

TVCLm=60.29
CLm=TVCLm*EXP(ETA(2))
Vm=THETA(9)

MTT =THETA(10)
N=2
KTR=N/MTT; transit rate constant

IF(DIST.LE.10) THEN 
DISTN = 10
ELSE
DISTN = DIST
ENDIF

DISTCOVGAMMA = THETA(11)
TVKPL6 = THETA(12) 
KPL6 = TVKPL6
TVPC6 = EXP(THETA(13))/DISTN**DISTCOVGAMMA
PC6 = TVPC6

DISTCOVGAMMAM = THETA(14)
TVKPL7 = THETA(15) 
KPL7 = TVKPL7
TVPC7 = EXP(THETA(16))/DISTN**DISTCOVGAMMAM
PC7 = TVPC7

$DES  
DADT(1)= -KTR*A(1)
DADT(2)=A(6)*KA-(VMAX*A(2)/V)/(KM+A(2)/V)-A(2)*K12+A(3)*K21
DADT(3)=A(2)*K12-A(3)*K21
DADT(4) = Fm*(CL/V)*A(2) - (CLm/Vm)*A(4) 
DADT(5)=KTR*A(1)-KTR*A(5)
DADT(6)=KTR*A(5)-KA*A(6) 
DADT(7)=KPL6 * (PC6*A(2)/V - A(7))  ; dCp,site/dt = KPL(PC * Cp,plasma - Csite)
DADT(8)=KPL7 * (PC7*A(4)/Vm - A(8))

$ERROR         
CP=A(2)/V
CM=A(4)/Vm ; metabolite prediction

DEL = 1E-12
IPRED=LOG(A(7)+DEL) ; parent lesion
IF(L2.EQ.2) IPRED=LOG(A(8)+DEL) ; metabolite lesion

A6 = A(7)
A7 = A(8)

W = SQRT((THETA(1)/EXP(IPRED))**2 + THETA(2))
IF(L2.EQ.2) W = SQRT((THETA(17)/EXP(IPRED))**2 + THETA(18))
IF(W.EQ.0) W = 1

IRES=DV-IPRED 
IWRES=IRES/W 

;Y = IPRED + W*EPS(1) + EPS(2)
;IF(L2.EQ.2) Y = IPRED+W*EPS(3) + EPS(4)

;M3 method
;LLOQ=LOG(1/1000)
;IF(CMPD.EQ.2) LLOQ=LOG(10/1000)
;DUM  = (LLOQ-IPRED)/W
;CUMD = PHI(DUM)
;IF (BLQ.EQ.0) THEN
;  F_FLAG = 0
  Y = IPRED + W*ERR(1)
;ENDIF

;IF(BLQ.EQ.1) THEN
;  F_FLAG = 1
;  Y = CUMD
;ENDIF

$THETA  0 FIX ; 1 Prop error --- went to zero
 (0,0.912115) FIX; 2 Add error Error.
 (0,207.2) FIX ; 3 V L
 (0,1.394) FIX ; 4 KA 1/h
 (0,1324) FIX ; 5 V2
 (0,40.43) FIX ; 6 Q
 (0,1) FIX ; 7 Fm
 (0,0.1949) FIX ; 8 KM
 (0,607.9) FIX ; 9 Vm
 (0,1.434) FIX ; 10 MTT
 (0,1.98691) FIX; 11 DISTCOVGAMMA
 (0,0.356753,10) ; 12 KPL6
 (0,12.73) FIX; 13 MAXPC
 (0,2.81318) FIX ; 14 DISTCOVGAMMAM
 (0,0.5,10) FIX ; 15 KPL7
 (0,19.5109) FIX ; 16 MAXPC7
 0 FIX ; 17 prop error metabolite
 (0,3.965) FIX ; 17 additive error metabolite
$OMEGA  
 0.089401  FIX
 0.123201  FIX
$SIGMA  1  FIX
$SIM (12345) (54321) ONLYSIM NSUB=500
;$ESTIMATION METHOD=1 INTERACTION LAPLACIAN PRINT=1 NSIG=2 SIGL=6
;            MAXEVAL=9999 NOABORT
;$COVARIANCE UNCONDITIONAL
$TABLE      ID TIME TAD DIST CMPD DOSE_FACTOR IPRED PRED DV CP CM A6
            A7 iCL WT AMT WRES IWRES CWRES EVID PC6 PC7 NOPRINT
            ONEHEADER FILE=sdtabD86_sim

