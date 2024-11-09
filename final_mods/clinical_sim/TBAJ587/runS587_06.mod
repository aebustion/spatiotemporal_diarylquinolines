;; 1. Based on: run102
;; 2. Description: Ulrika's S587 popPK model transcribed for lesion simulations
;; x1. Author: aebustion
$PROBLEM    TBAJ587 LESION Simulations

$INPUT ID TIME DV EVID AMT ADDL II DOSE DIST CMT MONTHS

$DATA      ./20240502_TBAJ587_DIST_therapeutic_scenarios.csv 
IGNORE=@
;IGNORE(DIST.NE.300) ; only edge
IGNORE(TIME.GT.5376) ; only eight months

$SUBROUTINE ADVAN13 TOL=9
$MODEL      NCOMP=11
            ;TBAJ587
            COMP  = (ABS,DEFDOSE) ;1
            COMP  = (CENTRAL)     ;2
            COMP  = (PERIPHERAL)  ;3
            COMP  = (PERIPHERAL2) ;4
            COMP  = (M31)         ;5
            COMP  = (M32)         ;6
            COMP  = (M33)         ;7
            COMP  = (M21)         ;8
            COMP  = (M22)         ;9
            COMP  = (PLESION)     ;10
            COMP  = (M3LESION)    ;11

$PK
FED=1
; Make sure dose gets switched to umol

;-------- PK (from Ulrika's group, REF: PAGE poster) --------

; Parent TBAJ-587
TVCL = 5.4 * ((DOSE/200)**(0.298))       ; L/h              
TVV  = 88.2               ; L
TVKA = 0.0866 * exp(-0.000462 * (DOSE-200))         ; 1/H 
TVV2 = 1910                ; L
Q  = 31.1           ; L/h 
TVV3 = 22500                ; L
Q2 = 31.5         ; L/h
IF(FED.EQ.0) TVF1 = 1                ; F1 for fasted
IF(FED.EQ.1) TVF1 = 1 * (1 + 0.688)  ; F1 for fed
IF(FED.EQ.0) TVMTT = 0.837           ; MTT for fasted
IF(FED.EQ.1) TVMTT = 0.958           ; MTT for fed
NN = 2.36 


CL = TVCL*EXP(ETA(1))
V  = TVV*EXP(ETA(2))
KA = TVKA*EXP(ETA(3))
V2 = TVV2*EXP(ETA(4))
V3 = TVV3*EXP(ETA(5))

;S2 = V

; AMT in mg in input file, MW TBAJ587 614.5 g/mol, DV as nmol/L
;(AMT/1000mg)/(614.5g)*1000000000 = AMT*1627.3393 nmol
TVF1mol = 1627.3393*TVF1
F1 = TVF1mol*EXP(ETA(6))

MTT = TVMTT*EXP(ETA(7))

; M3
  IF(FED.EQ.0) Fm3 = 1 * ((DOSE/200)**(-0.418)) * (1) ; for fasted
  IF(FED.EQ.1) Fm3 = 1 * ((DOSE/200)**(-0.418)) * (1 - 0.479) ; for fed
  TVCLM3 = 18.5 ; L/h
  VM3 = 816 ; L
  TVV2M3 = 1630 ; L
  QM3 = 9.09 ; L/h
  V3M3 = 3010 ; L
  Q2M3 = 99.6 ; L/h

  CLM3 = TVCLM3*EXP(ETA(8))
  V2M3 = TVV2M3*EXP(ETA(9))


; M2 ignoring IIV findings for now
  IF(FED.EQ.0) Fm2 = 1 * ((DOSE/200)**(-0.373)) * (1) ; for fasted
  IF(FED.EQ.1) Fm2 = 1 * ((DOSE/200)**(-0.373)) * (1 - 0.547) ; for fed
  CLM2 = 34.1 * ((DOSE/200)**(0.146))  ; L/h
  VM2 = 247 ; L
  V2M2 = 12200 ; L
  QM2 = 174 ; L/h

;LESION from runG54.mod
;parent lesion
DISTN = DIST
MAXPC = 11.29
DISTCOVGAMMA = 1.641
KPL10 = 0.006014 
PC10 = exp(MAXPC)/DISTN**DISTCOVGAMMA

;metabolite lesion
MAXPCMET = 12.97
DISTCOVMETABGAMMA = 1.868
KPL11 = 0.01166 
PC11 = exp(MAXPCMET)/DISTN**DISTCOVMETABGAMMA

; Sterling-Silver
IF(AMT.GT.0.AND.CMT.EQ.1) PD   = AMT
IF(AMT.GT.0.AND.CMT.EQ.1) TDOS = TIME
TAD = TIME-TDOS

L = LOG(2.5066)+(NN+.5)*LOG(NN)-NN+LOG(1+1/(12*NN)) ;transit approximation      
KTR=(NN+1)/MTT     

$DES 
X=.00001
IF(T.GE.TDOS)THEN
DADT(1)  = EXP(LOG(PD+X) + LOG(KTR+X) + NN*LOG(KTR*(T-TDOS)+X) - KTR*(T-TDOS)-L) - KA*A(1)
ELSE
DADT(1)  = EXP(LOG(PD+X) + LOG(KTR+X) + NN*LOG(KTR*T+X)        - KTR*T-L)        - KA*A(1)
ENDIF
DADT(2)  =  KA*A(1) - CL/V*A(2) - Q/V*A(2) + Q/V2*A(3) - Q2/V*A(2) + Q2/V3*A(4) - Fm3*CL/V*A(2) - Fm2*CL/V*A(2)
DADT(3)  =                        Q/V*A(2) - Q/V2*A(3)
DADT(4)  =                                               Q2/V*A(2) - Q2/V3*A(4)
DADT(5)  =                                                                        Fm3*CL/V*A(2) - CLM3/VM3*A(5) - QM3/VM3*A(5) + QM3/V2M3*A(6) - Q2M3/VM3*A(5) + Q2M3/V3M3*A(7)
DADT(6)  =                                                                                                        QM3/VM3*A(5) - QM3/V2M3*A(6)
DADT(7)  =                                                                                                                                       Q2M3/VM3*A(5) - Q2M3/V3M3*A(7)
DADT(8)  =                                                                                        Fm2*CL/V*A(2) - CLM2/VM2*A(8) - QM2/VM2*A(8) + QM2/V2M2*A(9)
DADT(9)  =                                                                                                                        QM2/VM2*A(8) - QM2/V2M2*A(9)

; LESION DES -----------------------
CP = (A(2)/V)*(614.5)/1000000 ; need to covert from nmol/L to mg/L before lesion
CM = (A(5)/VM3)*(600.5)/1000000 ; need to convert from nmol/L to mg/L before lesion

DADT(10)= KPL10 * (PC10 * CP - A(10))
DADT(11)= KPL11 * (PC11 * CM - A(11))


$ERROR 
;TBAJ587 output
A2 = A(2)
A5 = A(5)
A10 = A(10)
A11 = A(11)

DEL= 1E-12
IPRED=LOG(A(2)/V+DEL)

CPP = LOG((A(2)/V * 614.5 / 1000000) + DEL)
CPM = LOG((A(5)/VM3 * 600.5 / 1000000) + DEL) 
AA10= LOG(A(10) + DEL)
AA11= LOG(A(11) + DEL)

; additive error for parent plasma
RUV1=0.185 ; might need to be square root of this value?
W = RUV1

IRES=DV-IPRED 
IWRES=IRES/W

Y = IPRED + W*EPS(1) 

$THETA  
0 FIX ;DUMMY

$OMEGA 
0.178 FIX; CL
0.232 FIX; V
0.036 FIX; KA
0.148 FIX; V2
0.118 FIX; V3
0.176 FIX; F1
0.138 FIX; MTT
0.050 FIX; CLM3
0.099 FIX; V2M3

$SIGMA  
1  FIX 

$SIM (12345) (54321) ONLYSIM NSUB=500

$TABLE ID DOSE DIST TIME DV IPRED PRED CPP CPM A2 A5 A10 A11 AA10 AA11 MONTHS EVID FILE=sdtabS587_06 NOPRINT ONEHEADER
