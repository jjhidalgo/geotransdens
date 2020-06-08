      SUBROUTINE WATVOL_INI (INORPAR   ,IOPTS   ,IOVRWC   ,LMXNDL
     &                      ,LXPAREL   ,NPAREL   ,NOPTS
     &                      ,NPBMX ,NPBTP
     &                      ,NTYPAR  ,NUMEL
     &                      ,NZPAR     ,LNNDEL
     &                      ,ACTH      ,CFPAREL
     :                      ,PARC      ,WATVOL)
      
********************************************************************************
*
* PURPOSE
*
*  Computes initial volume of water per unit of aquifer.
*
*
* DESCRIPTION
*
*  Manages de the computation of initial WATVOL vector.
*
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  LNNDEL                 Number of nodes at every element
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter
*  ACTH                   Aquifer thickness of every element. Cross sectional   
*                         area for 1-D elements, thickness for 2-D elements.
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  PARC                   Vector containing calculated values for all           
*                         parameters
*  WATVOL                 Array containing the water content of every element
*                         The array stores the water content in time k+1 and
*                         k+theta.
*                             WATVOL(#,#,1) --> Time k.
*                             WATVOL(#,#,2) --> Time k+theta
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  LMXNDL                 Maximum number of nodes per element                   
*  NUMEL                  Number of elements                                    
*  IOVRWC                 Equal to IOPTS(31).
*                             1.WATVOL calculated elementwise
*                             2.WATVOL calculated nodewise.
*  NPAREL                 Number of element parameters in current problem       
*  NPBMX                  Max(NPBFL,NPBTP). Used to dimension some arrays
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*



* INTERNAL VARIABLES: SCALARS
*
*  H_AVG
*  VAUX
*  I
*  J
*  K
*  L                      Current element
*  NNUD                   Number of nodes of the current element                
*
*

* HISTORY
*
*     JHG      11-2003        First coding.
*                 
*******************************************************************************
      IMPLICIT NONE
      
C------------------------- External       

      INTEGER::LMXNDL, NUMEL, IOVRWC, NOPTS,NPAREL,
     &        NPBMX, NPBTP,NTYPAR, NZPAR


      INTEGER::IOPTS(NOPTS),INORPAR(NTYPAR),LNNDEL(NUMEL)
     &        ,LXPAREL(NUMEL,NPAREL,NPBMX)


      REAL*8::ACTH(NUMEL), CFPAREL(NUMEL,NPAREL),
     &        PARC(NZPAR)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3,NPBTP)

C------------------------- Internal

      INTEGER*4::IPROB,K, L,NNUD,NTP_SIM
      REAL*8:: VAUX

C-------------------------  First executabel estatement.

      NTP_SIM = MAX(1,IOPTS(29)*NPBTP)

C------------------------- Initial water volume

      IF (IOVRWC.LE.1) THEN

          DO IPROB=1,NTP_SIM

              DO L=1,NUMEL

                  WATVOL(1,L,1,IPROB) = CFPAREL(L,7)*
     &            PARC( INORPAR(15)+LXPAREL(L,7,IPROB) )*ACTH(L)

                  WATVOL(1,L,2,IPROB) = WATVOL(1,L,1,IPROB)
                  WATVOL(1,L,3,IPROB) = WATVOL(1,L,1,IPROB)

              END DO !L=1,NUMEL

          END DO !IPROB=1,NTP_SIM

      ELSE

          DO IPROB=1,NTP_SIM

              DO L=1,NUMEL

                  NNUD = LNNDEL(L)
                  VAUX = CFPAREL(L,7)
     &                   *PARC(INORPAR(15)+LXPAREL(L,7,IPROB))*ACTH(L)

                  DO K=1,NNUD
                      WATVOL(K,L,1,IPROB) = VAUX
                      WATVOL(K,L,2,IPROB) = WATVOL(K,L,1,IPROB)
                      WATVOL(K,L,3,IPROB) = WATVOL(K,L,1,IPROB)
                  END DO !K=1,NNUD

              END DO !L=1,NUMEL

          END DO !IPROB=1,NTP_SIM

      END IF !IOVRWC.LE.1

      END SUBROUTINE WATVOL_INI
