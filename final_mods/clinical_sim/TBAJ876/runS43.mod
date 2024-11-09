;; 1. Based on: run0444 and S32
;; 2. Description:  
;; x1. Author: Connie and Annamarie
$PROBLEM    TBAJ876 PK HUMAN
$INPUT      ID TIME DV EVID DOSE AMT II ADDL DIST MONTHS
$DATA      ./20240502_TBAJ876_DIST_therapeutic_scenarios200.csv 
		IGNORE=@ 
		;IGNORE(DIST.NE.300) ; only edge
		IGNORE(TIME.GT.5376) ; only eight months
$SUBROUTINE ADVAN6 TOL=5
$MODEL      NCOMP=12 
		COMP=(ABS,DEFDOSE) 
		COMP=(TRNS1) 
		COMP=(TRNS2)
            COMP=(CENTRAL,DEFOBS) 
            COMP=(PERIPHERAL1) 
            COMP=(PERIPHERAL2) 
            COMP=(METAB) 
            COMP=(PERIPHERAL1META) 
            COMP=(PERIPHERAL2META) ; major metabolite 
            COMP=(AUC)
            COMP=(PLESION) 
            COMP=(MLESION)
$PK 
FED=1      
CL = THETA(2)*EXP(ETA(4)) ; Clearance, L/h.
V1 = THETA(1)*EXP(ETA(2)) ; Central Volume, L.
KE = CL/V1 ; Elimination rate, 
KA = THETA(3) ; Absorption rate, 1/h
V2 = THETA(5) ; Peripheral Volume, L.
Q = THETA(4) ; Q
V3 = THETA(8); Peripheral 2 Volume, L.
Q2 = THETA(9)*EXP(ETA(5)) ; Q
K12 = Q/V1
K13 = Q2/V1
K21 = Q/V2
K31 = Q2/V3 

TVF1= THETA(10)*(1+FED*THETA(11)) ; Bioavailability
F1 = TVF1*EXP(ETA(1))

MTT=THETA(12)*EXP(ETA(3)) ; Transit model
N=2
KTR= (1+N)/MTT

FM=1 ; Fraction of drug metabolized

S2 = V1 ;scaling

TVCLm=THETA(13) ; CLm, clearance of metabolite
CLm=TVCLm

TVVm=THETA(14)*EXP(ETA(6)) ; Vm, Central volume of metabolite
Vm=TVVm

Qm = THETA(17) ; Q for metabolite

Vpm = THETA(18)*EXP(ETA(7)) ; Peripheral Volume, L, for metabolite

K14 = Qm/vm
K41 = Qm/Vpm

Q2m = THETA(19)*EXP(ETA(9)) ; Q2 for metabolite 

V2pm = THETA(20)*EXP(ETA(8)) ; Second Peripheral Volume, L, for metabolite
K15 = Q2m/vm
K51 = Q2m/V2pm

;Cmax
;IF(NEWIND.LE.1) THEN
;COM(1) = -1 
;ENDIF

IF(DIST.LE.10) THEN 
DISTN = 10
ELSE
DISTN = DIST
ENDIF

; values estimated from runD86.mod
DISTCOVGAMMAP = 1.987
TVKPL11 = 0.3568
KPL11 = TVKPL11
TVPC11 = EXP(12.73)/DISTN**DISTCOVGAMMAP
PC11 = TVPC11

; values estimated from runD87.mod
DISTCOVGAMMAM = 2.175
TVKPL12 = 0.3489
KPL12 = TVKPL12
TVPC12 = EXP(14.05)/DISTN**DISTCOVGAMMAM
PC12 = TVPC12                                                                                                        

$DES      
DADT(1)= -KTR*A(1)
DADT(2)= KTR*A(1) - KTR*A(2) 
DADT(3)= KTR*A(2) - KA*A(3) 
DADT(4)= KA*A(3) - CL/V1*A(4) - K12*A(4) + K21*A(5) - K13*A(4) + K31*A(6)
DADT(5) = K12*A(4) - K21*A(5) 
DADT(6) = K13*A(4) - K31*A(6)
DADT(7) = FM*KE*A(4)- (CLm/Vm)*A(7) - K14*A(7) + K41*A(8) - K15*A(7) + K51*A(9)
DADT(8) = K14*A(7)- K41*A(8)
DADT(9) = K15*A(7) -K51*A(9)
DADT(10) = A(4)/V1
DADT(11)= KPL11 * (PC11*A(4)/V1 - A(11)) ; dCp,site/dt = KPL(PC * Cp,plasma - Csite) ;parent lesion
DADT(12)= KPL12 * (PC12*A(7)/Vm - A(12)) ; dCp,site/dt = KPL(PC * Cp,plasma - Csite) ; metabolite lesion

;Cmax 
;CP = A(4)/V1
;IF(CP.GT.COM(1)) THEN 
;COM(1) = CP
;ENDIF

$ERROR   
; Residual Error model  
;REPID =IREP
IPRED= A(4)/V1 ; predicting concentration for parent 
;IF(CMT.EQ.4) IPRED=A(7)/Vm ; predicting concentration for metabolite    

WA=THETA(6)  ; ADDITIVE for parent    
WP=THETA(7)   ; PROPORTIONAL for parent 
WAm=THETA(15)  ; ADDITIVE for metabolite    
WPm=THETA(16)   ; PROPORTIONAL for metabolite 

W = SQRT(WA**2+(WP*IPRED)**2)  ;the standard deviation of the residual error, W, for parent  
;IF(CMT.EQ.4) W=SQRT(WAm**2+(WPm*IPRED)**2) ;the standard deviation of the residual error, W, for metabolite

IRES=DV-IPRED
IWRES=IRES/W  ; Can be used to convert the residual to the weighted residual (IWRES) by dividing the residual by W 
Y = IPRED + W*EPS(1) ; EPS(1) represents the random error term assumed to have proportional variability.

;AUC = A(10)/S2 ;AUC
;CMAX = COM(1) ;CMAX

;;;;OLDER ERROR FROM ANNAMARIE
DEL = 1E-12

CPP = LOG((A(4)/V1) + DEL)
CPM = LOG((A(7)/Vm) + DEL) 
AA11= LOG(A(11) + DEL)
AA12= LOG(A(12) + DEL)

$THETA  
 (0,163) FIX; 1 V L
 (0,12.3) FIX; 2 CL L/h
 (0,0.283) FIX; 3 KA /h
 (0,20) FIX; 4 Q
 (0,16800) FIX; 5 v2
 (0,0.00711) FIX; 6 add
 (0,0.153) FIX; 7 prop
 (0,408) FIX; 8 V3
 (0,21.6) FIX; 9 Q2
 1 FIX ; 10 F1_FASTED
 (0,0.401) FIX; 11 F1_FED
 (0,0.878) FIX; 12 MTT
 (0,40.4) FIX; 13 CLm L/h
 (0,1300) FIX; 14 Vm L
 0 FIX ; 15 ADDITIVE for metabolite
 (0,0.135) FIX; 16 PROPORTIONAL for metabolite
 (0,69.2) FIX; 17 Q for metabolite
 (0,6120) FIX; 18 Peripheral Volume,L,for metabolite
 (0,378) FIX; 19 Second Q for metabolite
 (0,2400) FIX; 20 Second Peripheral Volume,L,for metabolite
$OMEGA  
 0.0676 FIX
 0.264 FIX
 0.372 FIX
 0.236 FIX
 0.0374 FIX
 0.393 FIX
 0.428 FIX
 0.182 FIX
 0.447 FIX

$SIGMA  
1  FIX  ; eps residual variance

$SIM (12345) (54321) ONLYSIM NSUB=500
;$COVARIANCE PRINT=E
;$ESTIMATION METHOD=1 INTERACTION PRINT=5 MAXEVAL=9999
;$TABLE      ID TIME TAD DV EVID CMT DOSE BQL LLOQ FED COHORT DAY AMT
;            ADDL II MDV IPRED PRED DV CL V1 KE KA V2 Q V3 Q2 F1 K12
;            K13 K21 K31 IWRES IRES FILE=sdtabS01 NOPRINT ONEHEADER
;$TABLE ID REPID TIME EVID DV PRED IPRED CP CM A8 WTTIMEE AA9 AA10 DIST V 
;		VE VM2 VM2E RES IRES WRES IWRES CWRES NPDE NOPRINT 
;		ONEHEADER NOAPPEND FILE=sdtabS01
$TABLE ID DOSE DIST TIME DV IPRED PRED CPP CPM AA11 AA12 MONTHS EVID FILE=sdtabS43 NOPRINT ONEHEADER


