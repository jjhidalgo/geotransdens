       SUBROUTINE DER_ATRA_VD 
     ; (DFML     ,DSLL     ,DSTL     ,IDIMBB   ,IDIMDQ   ,IDIMQ                `
     ; ,IODIM    ,L        ,LANI     ,LDIM     ,NNUD     ,NUMEL
     ; ,XN       ,BIBI     ,DAT_VD   ,DISPER   ,QXYZ     ,VD)

*****************************************************************************
*
* PURPOSE
*      Computes derivatives of ATRA matrix w.r.t. Darcy's flow components
*
* DESCRIPTION
*
*      Computes derivatives of ATRA matrix w.r.t. Darcy's flow components
*      and stores the result in DAT_VD (IODIM,IDIMDQ,NUMEL). 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BIBI                   Array containing the product of interpolation
*                         functions gradient, for a given element
*  DAT_VD                 Derivatives of ATRA matrix w.r.t. Darcy's flow 
*                         components
*  DISPER                 Contains the diagonal terms of dispersivity tensor
*  QXYZ                   Products between the different components of
*                         Darcy's velocity divided by its norm
*  VD                     Darcy's velocity
*
* INTERNAL VARIABLES: ARRAYS
*
*  DAUX                   Auxiliar array used to store derivatives of 
*                         dispersion tensor components w.r.t. Darcy's velocity 
*                         components
*
* EXTERNAL VARIABLES: SCALARS
*
*  DFML                   Molecular diffusion of element L
*  DSLL                   Longitudinal dispersivity of element L
*  DSTL                   Transversal dispersivity of element L
*  IDIMBB                 Used to dimension array BIBI. Is equal to IDIMQ times 
*                         the maximum possible anisotropy of the problem
*  IDIMDQ                 Used to dimension array DAT_VD (second dimension). It
*                         is the number of different terms (ATRA)ij varying 
*                         i and  j in the local numeration of one element. 
*                         Is equal to LMXNDL*(LMXNDL-1)/2
*  IDIMQ                  Used to dimension array QXYZ. It is equal to
*                         IODIM*(IODIM+1)/2
*  IODIM                  Maximum dimension of the problem                                        
*  L                      Element number
*  LANI                   Anisotropy of the element
*  LDIM                   Vector containing the dimension of each element       
*  NNUD                   Number of nodes of the current element                
*  NUMEL                  Number of elements                                    
*  XN                     Norm of Darcy's velocity array
*
* INTERNAL VARIABLES: SCALARS
*
*  IC                     Counter to locate DAT_VD array components
*  IND                    Counter to locate BIBI array components
*  XN2                    Square of the norm of Darcy's velocity array
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*       None
*
* HISTORY
*
*       First coding AMS 1-02
*
*****************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       DIMENSION BIBI(IDIMBB,NUMEL),DAT_VD(IODIM,IDIMDQ,NUMEL)
     ;      ,DAUX(6,3),DISPER(3),VD(IODIM,NUMEL),QXYZ(IDIMQ,NUMEL)

C------------------------- Substracts molecular diffusion from the diagonal 
C------------------------- of dispersivity tensor

       DO I=1,LDIM
          DISPER(I)=DISPER(I)-DFML
       ENDDO

C------------------------- Computes derivatives of dispersivity tensor w.r.t 
C------------------------- Darcy's velocity components (qx, qy, qz) and is 
C------------------------- stored in DAUX array

       IF (LDIM.EQ.1) THEN
          DAUX(1,1)=DSLL*VD(1,L)/XN
       ELSE IF (LDIM.EQ.2) THEN
          DAUX(1,1)=VD(1,L)/XN*( 2*DSLL-DISPER(1)/XN )
          DAUX(1,2)=VD(2,L)/XN*( 2*DSTL-DISPER(1)/XN )
          DAUX(2,1)=VD(1,L)/XN*( 2*DSTL-DISPER(2)/XN )
          DAUX(2,2)=VD(2,L)/XN*( 2*DSLL-DISPER(2)/XN )
          DAUX(3,1)=VD(2,L)/XN*(DSLL-DSTL) * ( 1.D0-QXYZ(1,L)/XN )
          DAUX(3,2)=VD(1,L)/XN*(DSLL-DSTL) * ( 1.D0-QXYZ(2,L)/XN )
       ELSE
          XN2=XN*XN
          DAUX(1,1)=VD(1,L)/XN*( 2*DSLL-DISPER(1)/XN )
          DAUX(1,2)=VD(2,L)/XN*( 2*DSTL-DISPER(1)/XN )
          DAUX(1,3)=VD(3,L)/XN*( 2*DSTL-DISPER(1)/XN )

          DAUX(2,1)=VD(1,L)/XN*( 2*DSTL-DISPER(2)/XN )
          DAUX(2,2)=VD(2,L)/XN*( 2*DSLL-DISPER(2)/XN )
          DAUX(2,3)=VD(3,L)/XN*( 2*DSTL-DISPER(2)/XN )

          DAUX(3,1)=VD(1,L)/XN*( 2*DSTL-DISPER(3)/XN )
          DAUX(3,2)=VD(2,L)/XN*( 2*DSTL-DISPER(3)/XN )
          DAUX(3,3)=VD(3,L)/XN*( 2*DSLL-DISPER(3)/XN )

          DAUX(4,1)=VD(2,L)/XN*(DSLL-DSTL) * ( 1.D0-QXYZ(1,L)/XN )
          DAUX(4,2)=VD(1,L)/XN*(DSLL-DSTL) * ( 1.D0-QXYZ(2,L)/XN )
          DAUX(4,3)=-VD(3,L)/XN2*(DSLL-DSTL) * QXYZ(4,L)

          DAUX(5,1)=VD(3,L)/XN*(DSLL-DSTL) * ( 1.D0-QXYZ(1,L)/XN )
          DAUX(5,2)=-VD(2,L)/XN2*(DSLL-DSTL) * QXYZ(5,L)
          DAUX(5,3)=VD(1,L)/XN*(DSLL-DSTL) * ( 1.D0-QXYZ(3,L)/XN )

          DAUX(6,1)=-VD(1,L)/XN2*(DSLL-DSTL) * QXYZ(6,L)
          DAUX(6,2)=VD(3,L)/XN*(DSLL-DSTL) * ( 1.D0-QXYZ(2,L)/XN )
          DAUX(6,3)=VD(2,L)/XN*(DSLL-DSTL) * ( 1.D0-QXYZ(3,L)/XN )
       END IF


C-------------------------           d ATRA   d ATRA   d ATRA
C------------------------- Computes  ------ , ------ , ------  
C-------------------------            d qx     d qy     d qz

       IND=0
       IC=1
       DO I=1,NNUD-1
          DO J=I+1,NNUD
             DO N=1,LDIM
                S=0.D0
                DO K=1,LANI
                   S=S+DAUX(K,N)*BIBI(IND+K,L)
                END DO
                DAT_VD(N,IC,L)=S
             END DO
             IND=IND+LANI
             IC=IC+1
          END DO
       END DO
       RETURN
       END
