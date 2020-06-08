      SUBROUTINE RHS_IN_CHP
     &          (AFLU     ,CFCHP    ,DERH     ,IBCOD    ,IDIMAFLU
     &          ,IDIMDERH ,INEW     ,IVCHP    ,IXCHP    ,KXX
     &          ,LMXNDL   ,LNNDEL   ,NPARF    ,NUMEL    ,NUMNP
     &          ,NZCHP    ,IDIMWGT  ,NZPAR    ,IPNT_PAR)

*****************************************************************************
*
* PURPOSE
*    Adds the correction of the derivatives with respect to presc. head RHS 
*
* DESCRIPTION
*    Does the correction required for the computation of derivatives with
*    respect to prescribed head due to symmetrization of flow system matrix
*    at BC type 1. The needed expression in
*                   d h_l
*    SUM AFLU_il * ------    where the SUM is extended over all nodes l
*                   d Hz
*                            connected to node i and such that IBCOD(l)=1
*                            and z is the zone of H. Of course, if z is not 
*                            the presc. head zone to which node l belongs to, 
*                            the derivative is equal to zero.  
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AFLU                   Matrix of finite elements equations for flow problem  
*                         No boundary conditions are included on it.            
*  CFCHP                  Steady state prescribed head nodal coefficient.       
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  IVCHP                  Presc. head estimation index (0=no estimation)        
*  IXCHP                  Presc. head (steady) zone number at a given node      
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  INEW                   Index that controls which component of DERH are
*                         the RHS (DERH(.,.,1) OR DERH(.,.,2))
*  NPARF                  Number of transient parameters to be estimated        
*  NUMNP                  Number of nodes                                       
*  NZCHP                  Number of prescribed head zones                       
*
* INTERNAL VARIABLES: SCALARS
*
*  I
*  IJPOS
*  INODE
*  J
*  JNODE
*  JP
*  JZON
*  L
*  NNUD
*
* HISTORY
*
*     AMS     11-1998     First coding
*     JHG     09-2005     AFLU stored elemenwise without diagonal
*
********************************************************************************


      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMAFLU ,IDIMDERH ,INEW     ,LMXNDL   ,NPARF
     &          ,NUMEL    ,NUMNP    ,NZCHP,IDIMWGT  ,NZPAR,IPOS

      INTEGER*4::IBCOD(NUMNP)       ,IVCHP(NZCHP)   ,IXCHP(NUMNP)
     &          ,KXX(LMXNDL,NUMEL)  ,LNNDEL(NUMEL)
     ;     ,IPNT_PAR(NZPAR*IDIMWGT)

      REAL*8::AFLU(NUMEL,IDIMAFLU)  ,CFCHP(NUMNP)
     &       ,DERH(NUMNP,NPARF,IDIMDERH)

C------------------------- Internal

      INTEGER*4::I      ,IJPOS  ,INODE  ,IP     ,IZON   ,J      ,JNODE
     &          ,JP     ,JZON   ,L      ,NNUD

C------------------------- First executable statement

      DO L=1,NUMEL

         NNUD = LNNDEL(L)

         DO I=1,NNUD-1

             INODE = KXX(I,L)
             IZON = IXCHP(INODE)


             DO J=I+1,NNUD

                  IJPOS = (I-1)*NNUD + J - I*(I+1)/2

                  JNODE = KXX(J,L)
                  JZON = IXCHP(JNODE)

                 IF (JZON.GT.0) THEN

                     JP = IVCHP(JZON)

                     IF (JP.NE.0.AND.IBCOD(JNODE).EQ.1) THEN

                   IPOS=IPNT_PAR(JP)

                      DERH(INODE,IPOS,INEW) = DERH(INODE,IPOS,INEW)
     &                                      - AFLU(L,IJPOS)*CFCHP(JNODE)
                      
                      DERH(JNODE,IPOS,INEW) = DERH(JNODE,IPOS,INEW)
     &                                      + AFLU(L,IJPOS)*CFCHP(JNODE)

                     END IF !JP.NE.0.AND.IBCOD(JNODE).EQ.1

                  END IF !JZON.GT.0

                  IF (IZON.GT.0) THEN

                     IP = IVCHP(IZON)

                     IF (IP.NE.0.AND.IBCOD(INODE).EQ.1) THEN

                         IPOS=IPNT_PAR(IP)
                         DERH(JNODE,IPOS,INEW) = DERH(JNODE,IPOS,INEW)
     &                                      - AFLU(L,IJPOS)*CFCHP(INODE)

                         DERH(INODE,IPOS,INEW) = DERH(INODE,IPOS,INEW)
     &                                      + AFLU(L,IJPOS)*CFCHP(INODE)

                     END IF !IP.NE.0.AND.IBCOD(INODE).EQ.1

                 END IF !IZON.GT.0

             END DO !J=1,NNUD
         
         END DO !I=1,NNUD

      END DO !L=1,NUMEL


      END SUBROUTINE RHS_IN_CHP
