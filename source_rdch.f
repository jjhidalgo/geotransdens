       SUBROUTINE SOURCE_RDCH
     ; (LMXNDL   ,NPPEL    ,NUMEL    ,NUMNP    ,AREA     ,CAUX1
     ; ,KXX      ,LNNDEL   ,PAREL    ,SOURCE   ,WATVOL   ,IOVRWC)

********************************************************************************
*
* PURPOSE
*
*  Manages the computation of the contribution of a radioactive element to 
*  its "son"
*
*
* DESCRIPTION
*
*  Manages the computation of the contribution of a radioactive element to
*  its "son" by calling the routines that use WATVOL by nodes or elements 
*  as needed.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  SOURCE_RDCH_I          Computes the contribution with WATVOL by nodes
*  SOURCE_RDCH_L          Computes the contribution with WATVOL by elements
*
* HISTORY
*
*     AMS      7-2002     First coding
*
*******************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       IF (IOVRWC.LE.1) THEN

          CALL SOURCE_RDCH_L
     ; (LMXNDL   ,NPPEL    ,NUMEL    ,NUMNP    ,AREA     ,CAUX1
     ; ,KXX      ,LNNDEL   ,PAREL    ,SOURCE   ,WATVOL)

       ELSE

          CALL SOURCE_RDCH_I
     ; (LMXNDL   ,NPPEL    ,NUMEL    ,NUMNP    ,AREA     ,CAUX1
     ; ,KXX      ,LNNDEL   ,PAREL    ,SOURCE   ,WATVOL)

       ENDIF

       RETURN
       END

*****************************************************************************
*****************************************************************************

       SUBROUTINE SOURCE_RDCH_L
     ; (LMXNDL   ,NPPEL    ,NUMEL    ,NUMNP    ,AREA     ,CAUX1
     ; ,KXX      ,LNNDEL   ,PAREL    ,SOURCE   ,WATVOL)

*******************************************************************************
*
* PURPOSE
*
*     Computes the contribution of a radioactive element to its "son"
*
* DESCRIPTION
*
*     Computes the contribution of a radioactive element to its "son"
*     WATVOL is used by elements
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CAUX1                  Array containing concentrations, weighted by THETAT 
*                         time factor                                           
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*  SOURCE                                                                       
*  WATVOL                 Array containing the water content
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  LMXNDL                 Maximum number of nodes per element                   
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  L                      Current element
*  NNUD                   Number of nodes of the current element                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
*
* HISTORY
*
*     AMS      7-2002     First coding. Previously it was in SIM_JAC
*
*******************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       DIMENSION AREA(NUMEL), CAUX1(NUMNP), KXX(LMXNDL,NUMEL)
     ;      ,LNNDEL(NUMEL) ,PAREL(NUMEL,NPPEL) ,SOURCE(NUMNP)
     ;      ,WATVOL(NUMEL,3)

       CALL ZERO_ARRAY (SOURCE,NUMNP)

       DO L=1,NUMEL
          NNUD=LNNDEL(L)
          VAUX=AREA(L)*PAREL(L,13)*(WATVOL(L,2)+PAREL(L,14))/NNUD
          DO J=1,NNUD
             I=KXX(J,L)
             SOURCE(I)=SOURCE(I)+VAUX*CAUX1(I)
          ENDDO
       ENDDO

       RETURN
       END

*****************************************************************************
*****************************************************************************

       SUBROUTINE SOURCE_RDCH_I
     ; (LMXNDL   ,NPPEL    ,NUMEL    ,NUMNP    ,AREA     ,CAUX1
     ; ,KXX      ,LNNDEL   ,PAREL    ,SOURCE   ,WATVOL)

*******************************************************************************
*
* PURPOSE
*
*     Computes the contribution of a radioactive element to its "son"
*
* DESCRIPTION
*
*     Computes the contribution of a radioactive element to its "son".
*     WATVOL is used by nodes
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CAUX1                  Array containing concentrations, weighted by THETAT 
*                         time factor                                           
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*  SOURCE                                                                       
*  WATVOL                 Array containing the water content
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  LMXNDL                 Maximum number of nodes per element                   
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  L                      Current element
*  NNUD                   Number of nodes of the current element                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
*
* HISTORY
*
*     AMS      7-2002     First coding. Previously it was in SIM_JAC
*
*******************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       DIMENSION AREA(NUMEL), CAUX1(NUMNP), KXX(LMXNDL,NUMEL)
     ;      ,LNNDEL(NUMEL) ,PAREL(NUMEL,NPPEL) ,SOURCE(NUMNP)
     ;      ,WATVOL(LMXNDL,NUMEL,3)

       CALL ZERO_ARRAY (SOURCE,NUMNP)

       DO L=1,NUMEL
          NNUD=LNNDEL(L)
          VAUX=AREA(L)*PAREL(L,13)/NNUD
          DO J=1,NNUD
             I=KXX(J,L)
             SOURCE(I)=SOURCE(I)+VAUX*
     ;                  ( WATVOL(J,L,2)+PAREL(L,14) )*CAUX1(I)
          ENDDO
       ENDDO

       RETURN
       END
