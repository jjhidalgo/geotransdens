       SUBROUTINE DER_ATRA_FP_VD 
     ; (IDIMDQ   ,IODIM    ,LMXNDL   ,NUMEL    ,NUMNP    ,AREA                 `
     ; ,CCAL     ,DAT_VD   ,DVDP     ,GRDFF    ,KXX      ,LDIM                 `
     ; ,LNNDEL   ,RHS_IP)

*******************************************************************************
*
* PURPOSE
*
*   Computes derivatives of ATRA w.r.t flow parameters (only Darcy's velocity
*   dependent part) 
*
* DESCRIPTION
*
*   Computes derivatives of ATRA w.r.t flow parameters (only Darcy's velocity
*   dependent part). This subroutine takes derivatives of ATRA w.r.t. 
*   Darcy's velocity components, DAT_VD, the derivatives of Darcy's velocity
*   w.r.t. flow parameters and join them to form the derivatives of 
*   ATRA w.r.t. flow parameters.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,
*                         volume for 3-D)
*  CCAL                   Computed concentration at every node
*  DAT_VD                 Derivatives of ATRA w.r.t. the components of Darcy's
*                         velocity (qx, qy, qz).
*  DVDP                   Derivatives of Darcy's velocity with respect to
*                         flow parameters. In some cases it is dimensioned
*                         as (IODIM,NUMEL) to reduce storage.
*  GRDFF                  Array containing the product between interpolation
*                         functions integrals and interp. functions gradient
*  KXX                    Node numbers of every element (counterclockwise
*                         order).
*  LDIM                   Vector containing the dimension of each element
*  LNNDEL                 Number of nodes at every element
*  RHS_IP                 Right hand side for the derivatives
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMDQ                 Used to dimension array DAT_VD (second dimension). It
*                         is the number of different terms (ATRA)ij varying
*                         i and  j in the local numeration of one element.
*                         Is equal to LMXNDL*(LMXNDL-1)/2
*  IODIM                  Maximum dimension of the problem
*  LMXNDL                 Maximum number of nodes per element
*  NUMEL                  Number of elements
*  NUMNP                  Number of nodes
*
* INTERNAL VARIABLES: SCALARS
*
*  IC                     Index counter for DAT_VD
*  L                      Element number
*  NNUD                   Number of nodes of the current element                
*  ZERO                   Equal to 0.D0
*
* HISTORY
*
*  First coding AMS 1-02
*
*******************************************************************************


       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION  DAT_VD(IODIM,IDIMDQ,NUMEL),CCAL(NUMNP),
     ;      DVDP(IODIM,NUMEL),KXX(LMXNDL,NUMEL),
     ;      GRDFF(IODIM,LMXNDL,NUMEL),RHS_IP(NUMNP),LDIM(NUMEL),
     ;      AREA(NUMEL),LNNDEL(NUMEL)
       DATA ZERO /0.D0/

       DO L=1,NUMEL
          NNUD=LNNDEL(L)
          IC=1
          DO I=1,NNUD-1
             I1=KXX(I,L)
             DO J=I+1,NNUD
                I2=KXX(J,L)
                T1=ZERO
                T2=ZERO
                T3=ZERO
                DO N=1,LDIM(L)
                   T1=T1+ DAT_VD(N,IC,L) *DVDP(N,L)
                   T2=T2+ AREA(L)*GRDFF(N,J,L)/NNUD *DVDP(N,L)
                   T3=T3+ AREA(L)*GRDFF(N,I,L)/NNUD *DVDP(N,L)
                ENDDO
                S1=(T1+T2)*( CCAL(I2)-CCAL(I1) )
                S2=(T1+T3)*( CCAL(I1)-CCAL(I2) )
                RHS_IP(I1)=RHS_IP(I1)-S1
                RHS_IP(I2)=RHS_IP(I2)-S2
                IC=IC+1
             END DO       ! J
          END DO          ! I
       ENDDO              ! L, ELEMENTS
       RETURN
       END
