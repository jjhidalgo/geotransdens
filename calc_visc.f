      SUBROUTINE CALC_VISC
     &          (CCALIT   ,DER_VISC ,ITPTVAR  ,KXX    ,LMXNDL
     &          ,LNNDEL   ,NUMEL    ,NUMNP    ,VAR_REF,VISCOSITY)
 
**************************************************************
* PURPOSE:
*
* Manages the computation (elementwise) of the viscosity and 
* its derivative w. r. t. the state variable.
*
*
* DESCRIPTION
*
* Computes the viscosity and its derivative in every element.
* The expresions used to compute the viscosity depend on the meaning
* of the state variable.
*
* Viscosity is calculated according to the values of temperature
* and mass fraction. (Hassanizadeh and Leijnse, 1988):
*
*
*     VISC = VISC_T * (1 + VIC_C)
*
*     VISC_T = (2.1E-12) * EXP(1808.5/(273.15 + T))
*      VISC_C =	1.85c - 4.1 c^2 + 44.5c^3
*
* The viscosity is measured in (MPa.s).
* and the temperature in(C).
*
*
* Since the state variable can be either the temperature or the mass fraction,
* it is necesary to indicate which one is representated by CCALIT. This can be
* done with the ITPTVAR variable.
*
* VAR_REF represents the state variable tha is not changed, i. e., if solute
* transport is being solved, then ITPTVAR = 0, meaning CCALAN is  mass fraction,
* and VAR_REF is the temperature of the system. Similarly, if energy transport
* is being solved, ITPTVAR = 1, CCALAN is temperature and VAR_REF is the mass
* fraction of the fluid.
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CCALIT                Concetration (actually, mass fraction or temperature)
*                        to compute viscosity.
*  DER_VISC              Derivative of the viscosity w. r. t. the state variable.
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*
*  VAR_REF                Value of the reference variable.
*  VISCOSITY.             Viscosity in every element.
*
* EXTERNAL VARIABLES: SCALARS
*
*  ITPTVAR                Indicates the type of state variable of transport.
*                           0. Mass fraction
*                           1. Temperature.
*  LMXNDL                 Maximum number of nodes per element                   
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  I                        Node counter.
*  L                        Element counter.
*  MEANCONC                 mean CCALITration
*  NNUD                     Number of nodes of current eleement.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  VISC                     Computes DENSITY by nodes.
*
* HISTORY
*
*	JHG           10-2003        First coding.
*
*******************************************************************************


      IMPLICIT NONE

      INTEGER*4::I,ITPTVAR,L,LMXNDL,NNUD,NODE,NUMEL,NUMNP
      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)
     &           

      REAL*8::DERVISCDTRA,MEANCONC,VAR_REF,VISC
      REAL*8::CCALIT(NUMNP),DER_VISC(NUMEL),VISCOSITY(NUMEL)

C------------------------- For each element
      DO L=1,NUMEL

          NNUD=LNNDEL(L)
          MEANCONC=0D0
C------------------------- the mean CCALITration is calculated
          DO I=1,NNUD

              NODE=KXX(I,L)
              MEANCONC=MEANCONC+CCALIT(NODE)

          END DO !Next nope.

C------------------------- and the viscosity evaluated with the mean
C------------------------- concetration value.

          VISCOSITY(L)=VISC(ITPTVAR,MEANCONC,VAR_REF)
          DER_VISC(L)=DERVISCDTRA(ITPTVAR,MEANCONC,VAR_REF)

      END DO ! Next element.

      END SUBROUTINE CALC_VISC
*
**************************************************************
**************************************************************
*



      FUNCTION VISC(ITPTVAR,STVAR,VAR_REF) RESULT (VIS)

**************************************************************
* PURPOSE
*
* To calculate the viscosity according to the values of temperature
* and mass fraction.
*  
* DESCRIPTION
*
* Viscosity is calculated according to the values of temperature
* and mass fraction. (Hassanizadeh and Leijnse, 1988):
*
*
*        VISC = VISC_T * (1 + VIC_C)
*
*        VISC_T = (2.1E-6) * EXP(1808.5/(273.15 + T))
*        VISC_C =	1.85c - 4.1 c^2 + 44.5c^3
*
* The viscosity is measured in (MPa*s).
* and the temperature in(C).
*
*
* Since the state variable can be either the temperature or the mass fraction,
* it is necessary to indicate which one is represented by STVAR. This can be
* done with the ITPTVAR variable.
*
*
* ARGUMENTS
*
*  ITPTVAR     Indicates the type of state variable of transport.
*                           0. Mass fraction
*                           1. Temperature.
*  STVAR       Transport state variable.
*  VAR_REF     Reference variable.
*
*
* INTERNAL VARIABLE: SCALAR
*
*  C           Mass fraction
*  T           Temperature
*  VIS         Used to store the result of the function.
*  VIS_C       Contribution to viscosity due to mass fraction.
*  VIS_T       Contribution to viscosity due to temperature.
*
*
* REFERENCES
*
*  Hassanizadeh, S. M. and Leijnse, T.: 1988, On the modeling of brine transport
*  in porous media. Water Resour. Res. 24, 321-330.
*
* HISTORY
*
*     JHG     10-2003     First coding.
*
**************************************************************
      IMPLICIT NONE
      REAL*8::C,T,VIS,VIS_C, VIS_T,STVAR,VAR_REF
      INTEGER*4::ITPTVAR


      SELECT CASE(ITPTVAR)

          CASE(0) !Mass fraction

              C = STVAR
              T = VAR_REF

          CASE(1) !Temperature

              T = STVAR
              C = VAR_REF

      END SELECT

      VIS_T = (2.1D-6) * DEXP(1808.5D0/(273.15D0 + T))
      VIS_C = 1.85D0*C -4.1D0*C**2 + 44.5D0*C**3
      VIS = VIS_T * (1D0 + VIS_C)

      END FUNCTION VISC

**************************************************************
**************************************************************



      FUNCTION DERVISCDTRA(ITPTVAR,STVAR,VAR_REF) RESULT (DERVIS)

**************************************************************
* PURPOSE
*
* To calculate the derivative of the viscosity w. r. t. transport
* state variable.
*
* DESCRIPTION
*
* Viscosity is calculated according to the values of temperature
* and mass fraction. (Hassanizadeh and Leijnse, 1988):
*
*
*     VISC = VISC_T * (1 + VIC_C)
*
*     VISC_T = (2.1E-6) * EXP(1808.5/(273.15 + T))
*     VISC_C =	1.85c - 4.1 c^2 + 44.5c^3
*
* The viscosity is measured in (MPa*s).
* and the temperature in (C).
*
*
* Since the state variable can be either the temperature or the mass fraction,
* it is necessary to indicate which one is representated by STVAR. This can be
* done with the ITPTVAR.
*
*
*
* ARGUMENTS
*
*  ITPTVAR  Indicates the type of state variable of transport.
*                           0. Mass fraction
*                           1. Temperature.
*  STVAR    Transport state variable.
*
*  VAR_REF  Reference variable.
*
*
* INTERNAL VARIABLE: SCALAR
*
*  AUX    Contribution to viscosity due to VAR_REF.
*  C      Mass fraction
*  T      Temperature
*  VIS    Used to store the result of the function.
*
* REFERENCES
*
*  Hassanizadeh, S. M. and Leijnse, T.: 1988, On the modeling of brine transport
*     in porous media. Water Resour. Res. 24, 321-330.
*
* HISTORY
*
*      JHG         10-2003    First coding.
*
**************************************************************
      IMPLICIT NONE
      REAL*8::AUX,C, DERVIS,STVAR,T,VAR_REF
      INTEGER*4::ITPTVAR

      SELECT CASE(ITPTVAR)

          CASE(0) !Mass fraction

              C = STVAR
              T = VAR_REF

              AUX = (2.1D-6) * DEXP(1808.5D0/(273.15D0 + T))

              DERVIS = AUX * (1.85D0 - 8.2D0*C + 133.5D0*C**2 )

          CASE(1) !Temperature

              C = VAR_REF
              T = STVAR

              AUX = 1D0 + 1.85D0*C - 4.1D0*C**2 + 44.5D0*C**3

              DERVIS = AUX * (2.1D-6) * DEXP(1808.5D0/(273.15D0 + T)) *
     &                (-1808.5D0/(273.15D0 + T)**2)

      END SELECT

      END FUNCTION DERVISCDTRA
