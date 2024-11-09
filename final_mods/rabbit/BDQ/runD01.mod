;; 1. Based on: 16
$PROBLEM    rabbit PK gradient in caseum,BDQ
$INPUT      C ID LESION RABBITID OLDID=DROP TIME TAD TIMED DIST NGML
            NDV AMT DOSE FLAT_DOSE nDOSES EVID MDV II ADDL METABFL BLQ
            CMTFL CMTDONUT TISSUE SPECIES SURF_AREA_UM2 AREA SLICE_NO
            DIST_250GRP DIST_SEQGRP INFECTION WT iCL iV iKA iKE iQ iV2
            iF1 iCLM iVM iPC6 iKPL6 DV
$DATA      1_NMdata_BQDGradient_20220324_IDperDV.csv IGNORE=C
            IGNORE(DIST.GT.1500) IGNORE(BLQ.EQ.1)

;IGNORE(METABFL.EQ.1)
$SUBROUTINE ADVAN6 TOL=2
$MODEL      NCOMP=7 
            COMP=(ABS,DEFDOSE) 
            COMP=(CENTRAL,DEFOBS)
            COMP=(PERIPH) 
            COMP=(METAB)
            COMP=(TRANSIT1)
            COMP=(CASEUMPARENT) 
            COMP=(CASEUMMETAB)
$PK   

;IF(METABFL.EQ.0) L2 = 1
;IF(METABFL.EQ.1) L2 = 2


CL = 11.2
V  = 257*EXP(ETA(2))
KE = CL/V
KA = 0.411
S2 = V

F1 = 1*EXP(ETA(1))

Q = 26
V2 = 301

CLM = 27.2
VM = V
KM = CLM/VM

IF(DIST.LE.10) THEN
DISTN = 10
ELSE
DISTN = DIST
ENDIF

DISTCOV = THETA(3)
DISTCOVGAMMA = THETA(7)
TVKPL6 = THETA(1) 
KPL6 = TVKPL6
TVPC6 = EXP(THETA(2))/DISTN**DISTCOVGAMMA
PC6 = TVPC6

DISTCOVMETAB = THETA(6)
DISTCOVMETABGAMMA = THETA(8)
TVKPL7 = THETA(4) 
KPL7 = TVKPL7
TVPC7 = EXP(THETA(5))/DISTN**DISTCOVMETABGAMMA
PC7 = TVPC7

$DES         
DADT(1)=-KA*A(1)
DADT(2)= KA*A(5)-KE*A(2)-Q/V*A(2)+Q/V2*A(3)
DADT(3)= Q/V*A(2)-Q/V2*A(3)
DADT(4)= KE*A(2)-CLM*A(4)/VM
DADT(5)= KA*A(1)-KA*A(5)
DADT(6)= KPL6 * (PC6*A(2)/V - A(6))
DADT(7)= KPL7 * (PC7*A(4)/VM - A(7))


$ERROR            

CP = A(2)/V
CM = A(4)/VM

DEL = 1E-12
IPRED = LOG(A(6)+DEL)
IF(METABFL.EQ.1) IPRED = LOG(A(7)+DEL)

W = SQRT((THETA(9)/EXP(IPRED))**2 + THETA(10)**2)
IF(METABFL.EQ.1) W = SQRT((THETA(11)/EXP(IPRED))**2 + THETA(12)**2)
IF(W.EQ.0) W = 1


IRES = DV-IPRED
IWRES = IRES/W


LLOQ=LOG(1/1000)
IF(METABFL.EQ.1) LLOQ=LOG(10/1000)
DUM  = (LLOQ-IPRED)/W
CUMD = PHI(DUM)
IF (BLQ.EQ.0) THEN
  F_FLAG = 0
  Y = IPRED + W*ERR(1)
ENDIF

IF(BLQ.EQ.1) THEN
  F_FLAG = 1
  Y = CUMD
ENDIF


$THETA  
 (0,0.00320542) ; 1 KPL7
 (0,5.08718) ; 2 PC7
 (0,462.129) FIX ; 3 DISTCOV
 (0,0.00350731) ; 4 KPL8
 (0,7.04629) ; 5 PC8
 (0,308.283) FIX ; 6 DISTCOV
 (0,0.681802) ; 7 DISTCOVGAMMA
 (0,0.753408) ; 8 DISTCOVMETABGAMMA
 (0,0.0252316) ; 9 prop Error
 (0,1.5107) ; 10 add Error
 (0,0.0572238) ; 11 metab prop Error
 (0,1.17652) ; 12 metab add Error
$OMEGA
0.10099684 FIX
0.32695524 FIX
$SIGMA  1  FIX  ;    1 error
 1  FIX  ;    2 error
$ESTIMATION METHOD=1 INTERACTION LAPLACIAN NUMERICAL PRINT=1 SIGL=6
            NSIG=1 MAXEVAL=9999 NOABORT
$COVARIANCE UNCONDITIONAL
$TABLE      ID RABBITID LESION TIME TAD IPRED PRED DV METABFL AMT DIST
            WRES IWRES CWRES EVID F1 CL V KPL6 PC6 KPL7 PC7 WT
            FILE=sdtabD01 NOPRINT ONEHEADER

