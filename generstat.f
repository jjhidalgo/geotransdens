      SUBROUTINE GENSTAT
     ;(CPERO    ,DOF      ,FOVERDOF ,NPAR     ,NUMTOBS  ,OBJF)
   
********************************************************************************
*
* PURPOSE
*
* This subroutine calculates 3 general statistics over the 
* minimization process, such as the degrees of freedom
* and F/Dof
*
* 
* DESCRIPTION
*
* Step 1: Counts how many estimated parameters have prior info.
* Step 2: Calculates degrees of freedom as number of observations + number of
*         param. with prior information - number of parameters
* Step 3: F / (Nr of observations + Nr of parameters with prior information)
* Step 4: F/(degrees of freedom)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  PARM                   Vector containing measured values for all             
*                         parameters                                            
*
* EXTERNAL VARIABLES: SCALARS
*
*  CPERO                  Mean contribution od data to obj. func
*                         = Obj. Func / Ndata
*  DOF                    Degrees of freedom                     
*  FOVERDOF               Obj. function/ Degrees of freedom
*  NPAR                   Total number of parameters to be estimated            
*  NUMTOBS                Total number of observations                          
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  OBJF                   Objective function
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Dummy counter
*  NPRIOR                 Number of estim. param. with prior info.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

C_______________________________________________________Declaration of variables

      IMPLICIT NONE

                                                  !  External variables: scalars
      INTEGER*4 NUMTOBS, NPAR
      REAL*8 OBJF

                                                   ! Internal variables: scalars
      INTEGER*4 NPRIOR
      REAL*8 DOF, CPERO, FOVERDOF

C___________________________________Step 1: Calc. number of estimated parameters
C_______________________________________    and number of param. with prior info

      NPRIOR = NPAR

C____________________________________________Step 2: Calc. degrees of freedom =
C______number of observations + number of parameters with prior info.- number of
C___________________________________________________________________  parameters

      DOF = NUMTOBS + NPRIOR - NPAR 

C__________________Step 3: F / (Nr of observations + Nr of parameters with prior
C___________________________________________________________________information)

      CPERO = OBJF /(NUMTOBS + NPRIOR)

C_________________________________________________Step 4: F/(degrees of freedom)

      FOVERDOF = OBJF / DOF

      RETURN
      END
