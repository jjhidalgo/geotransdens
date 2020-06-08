       SUBROUTINE MODIF_JAC_C_RDCH
     & (IDIMDERC ,INEWT    ,IOLDT    ,LMXNDL   ,NPAR     ,NPPEL
     & ,NUMEL    ,NUMNP    ,THETAT   ,AREA     ,DERC     ,DERCS
     & ,KXX      ,LNNDEL   ,PAREL    ,WATVOL)

*******************************************************************************
*
* PURPOSE
*
*     Computes derivatives of the contribution of parent concentration w.r.t. 
*     estimated parameters
*
* DESCRIPTION
*
*     Computes part of the derivatives of the contribution of parent concentration w.r.t. 
*     estimated parameters to the RHS of the son. The contribution of parent 
*     is FOD*CRD*WATVOL*Cparent. The derivative consists of two parts:
*
*                          d Cparent     d(FOD*CRD*WATVOL)            
*        FOD*CRD*WATVOL*  ----------  +  ----------------- * Cparent
*                             dp                dp                    
*
*     Where Cparent is computed at k+thetat. The second term is computed is added in 
*     subroutines DERFOD, DERCRD, DERPOR when needed. The first term is computed in this
*     subroutine.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters ("father"). 
*  DERCS                  Nodal concentration derivatives with respect to       
*                         estimated parameters ("son").
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*  WATVOL                 Array containing the water content of every element
*
* EXTERNAL VARIABLES: SCALARS
*
*  INEWT                  Third dimension of array DERC used to store           
*                         derivatives of C w.r.t. parameters at next time step  
*  IOLDT                  Third dimension of array DERC used to store           
*                         derivatives of C w.r.t. parameters at previous time   
*                         step
*  LMXNDL                 Maximum number of nodes per element
*  NPAR                   Total number of parameters to be estimated            
*  NPPEL                  Total number of parameters by elements (not confuse
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)   
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* HISTORY
*
*     AMS     11-2003     First coding
*
*******************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       DIMENSION
     ;  DERC(NUMNP,NPAR,IDIMDERC) ,DERCS(NUMNP,NPAR,2)
     & ,PAREL(NUMEL,NPPEL)        ,AREA(NUMEL)
     & ,LNNDEL(NUMEL)             ,KXX(LMXNDL,NUMEL)
     & ,WATVOL(NUMEL)

C------------------------- First term

       DO L=1,NUMEL
          NNUD=LNNDEL(L)   

C------------------------- First order decay coefficient
C------------------------- FOD*CRD*B*AREA

          ALRD=PAREL(L,13)*( PAREL(L,14) + WATVOL(L) )*AREA(L)/NNUD
          DO I=1,NNUD
             I1=KXX(I,L)
             DO IP=1,NPAR
                DERCS(I1,IP,INEWT)=DERCS(I1,IP,INEWT)+ALRD* (
     ;                      THETAT*DERC(I1,IP,INEWT)+(1.D0-THETAT)*
     ;                                             DERC(I1,IP,IOLDT) )

             ENDDO

          ENDDO        ! I=1,NNUD

       ENDDO      ! Elements

       RETURN
       END
