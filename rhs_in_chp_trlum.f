      SUBROUTINE RHS_IN_CHP_TRLUM
     &          (AFLU     ,CFCHPT   ,DERH     ,DFLU     ,DTIM
     &          ,FNT      ,IBCOD    ,IDIMAFLU ,IDIMDERH ,IDIMDFLU
     &          ,IDIMFNT  ,INEW     ,INTI     ,IVCHP    ,IXCHPT
     &          ,KXX      ,LMXNDL   ,LNNDEL   ,NFTCHP   ,NINT
     &          ,NPARF    ,NUMEL    ,NUMNP    ,NZCHP    ,THETAF
     &          ,TINC     ,IDIMWGT  ,NZPAR    ,IPNT_PAR)

*****************************************************************************
*
* PURPOSE
*    Adds the correction of the derivatives with respect to presc. head RHS 
*    Transient, lumped case.
*
* DESCRIPTION
*    Does the correction required for the computation of derivatives with
*    respect to prescribed head due to symmetrization of flow system matrix
*    at BC type 1. The needed expression in
*                                         d h_l
*    SUM (THETAF*AFLU_il +DFLU_il/TINC)* ------    
*                                         d Hz
*                            where the SUM is extended over all nodes l
*                            connected to node i and such that IBCOD(l)=1
*                            and z is the zone of H. Of course, if z is not 
*                            the presc. head zone to which node l belongs to, 
*                            the derivative is equal to zero. D_il is zero
*                            always but when i=l.
*
* EXTERNAL VARIABLES: ARRAYS
*
*
*  AFLU                   Matrix of finite elements equations for flow problem  
*                         No boundary conditions are included on it.            
*  CFCHPT                 Transient prescribed head nodal coefficient.          
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  DFLU                   Matrix of finite elements equations for flow          
*                         problem related to storage term.                      
*  FNT                    Array containing time functions values                
*  IVCHP                  Presc. head estimation index (0=no estimation)        
*  IXCHPT                 Presc. head (transient) zone number at a given node   
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  DTIM                   Current time for the computation of flow time
*                         functions (counted since the beginning of the
*                         current observation interval, not since the beginning
*                         of the problem) divided by the length of the
*                         observation interval (interval between two
*                         consecutive observation times). Used only to make a
*                         linear interpolation of time function values.
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  INEW                   Index that controls which component of DERH are
*                         the RHS (DERH(.,.,1) OR DERH(.,.,2))
*  NBAND1                 Used to dimension. It is equal to NBAND+1             
*  NINT                   Number of observation times                           
*  NPARF                  Number of transient parameters to be estimated        
*  NUMNP                  Number of nodes                                       
*  NZCHP                  Number of prescribed head zones                       
*  THETAF                 Time weighting parameter for flow problems            
*
* INTERNAL VARIABLES: SCALARS
*
*  I                                                                            
*  IP                                                                           
*  IX                                                                           
*  K                                                                            
*
*
* HISTORY
*
*     AMS     11-1998     First coding
*
*****************************************************************************


      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMAFLU,IDIMDERH,IDIMDFLU,IDIMFNT,INEW,INTI,IPOS
     &          ,LMXNDL,NINT,NPARF,NUMEL,NUMNP,NZCHP,IDIMWGT  ,NZPAR

      REAL*8::DTIM,THETAF,TINC

      INTEGER*4::IBCOD(NUMNP) ,IXCHPT(NUMNP),IVCHP(NZCHP)
     &          ,KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),NFTCHP(NZCHP)
     ;     ,IPNT_PAR(NZPAR*IDIMWGT)

      REAL*8::AFLU(NUMEL,IDIMAFLU)       ,CFCHPT(NUMNP)
     &       ,DERH(NUMNP,NPARF,IDIMDERH) ,DFLU(NUMEL,IDIMDFLU)
     &       ,FNT(IDIMFNT,NINT)
      
C------------------------- Internal

      INTEGER*4::I      ,IJPOS  ,INODE  ,IP     ,IZON   ,J      ,JNODE
     &          ,JP     ,JZON   ,L      ,NFTI   ,NFTJ   ,NFTN   ,NNODE
     &          ,NNUD   ,NP     ,NZON

      REAL*8::DERI   ,DERJ   ,DERN   ,DFACTOR,FACTOR ,FT1I   ,FT1J
     &       ,FT1N   ,TINV


C------------------------- First executable statement

C------------------------- Computes inverse of TINC because multiplications are
C------------------------- performed faster than divisions

      TINV = 1D0/TINC

      DO L=1,NUMEL

         NNUD = LNNDEL(L)

         DO I=1,NNUD-1

             INODE = KXX(I,L)
             IZON = IXCHPT(INODE)

             IF (IZON.GT.0) THEN

                 IP = IVCHP(IZON)

C------------------------- DFLU contribution to diagonal.

                  IF (IP.NE.0.AND.IBCOD(INODE).EQ.1) THEN

                     DERI = CFCHPT(INODE)
                      NFTI = NFTCHP(IZON)

                     IF (NFTI.NE.0) THEN

                          FT1I = FNT(NFTI,INTI)
                          DERI =DERI*((FNT(NFTI,INTI+1)-FT1I)*DTIM+FT1I)

                      END IF !NFTI.NE.0

                     DFACTOR = TINV*DFLU(L,I)*DERI

                     IPOS=IPNT_PAR(IP)
                 DERH(INODE,IPOS,INEW) = DERH(INODE,IPOS,INEW) -DFACTOR

                  END IF !IP.NE.0.AND.IBCOD(INODE).EQ.1

             END IF !IZON.GT.0

C------------------------- AFLU contribution (taking into account AFLU storage).

             DO J=I+1,NNUD

                  IJPOS = (I-1)*NNUD + J - I*(I+1)/2

                  JNODE = KXX(J,L)
                  JZON = IXCHPT(JNODE)

                 IF (JZON.GT.0) THEN

                     JP = IVCHP(JZON)

                     IF (JP.NE.0.AND.IBCOD(JNODE).EQ.1) THEN

                         DERJ = CFCHPT(JNODE)
                          NFTJ = NFTCHP(JZON)

                         IF (NFTJ.NE.0) THEN

                              FT1J = FNT(NFTJ,INTI)
                              DERJ = DERJ*(
     &                                     (FNT(NFTJ,INTI+1)-FT1J)
     &                                     *DTIM+FT1J
     &                                    )

                          END IF !NFT.NE.0

                          FACTOR = THETAF*AFLU(L,IJPOS)*DERJ

                     IPOS=IPNT_PAR(JP)
                         DERH(INODE,IPOS,INEW) = DERH(INODE,IPOS,INEW)
     &                                          - FACTOR
                         DERH(JNODE,IPOS,INEW) = DERH(JNODE,IPOS,INEW)
     &                                          + FACTOR

                     END IF !JP.NE.0.AND.IBCOD(JNODE).EQ.1

                 END IF !JZON.GT.0

                  IF (IZON.GT.0) THEN

                      IF (IP.NE.0.AND.IBCOD(INODE).EQ.1) THEN

C------------------------- DERI has already been computed outside J loop.

                          FACTOR = THETAF*AFLU(L,IJPOS)*DERI

                     IPOS=IPNT_PAR(IP)
                         DERH(JNODE,IPOS,INEW) = DERH(JNODE,IPOS,INEW)
     &                                          - FACTOR
                         DERH(INODE,IPOS,INEW) = DERH(INODE,IPOS,INEW)
     &                                          + FACTOR

                     END IF !IP.NE.0.AND.IBCOD(INODE).EQ.1

                 END IF !IZON.GT.0

             END DO !J=I+1,NNUD
         
         END DO !I=1,NNUD-1

C------------------------- Last diagonal term

              NNODE = KXX(NNUD,L)
             NZON = IXCHPT(NNODE)
             
              IF (NZON.GT.0) THEN

                 NP = IVCHP(NZON)

                  IF (NP.NE.0.AND.IBCOD(NNODE).EQ.1) THEN

                     DERN = CFCHPT(NNODE)
                      NFTN = NFTCHP(NZON)

                     IF (NFTN.NE.0) THEN

                          FT1N = FNT(NFTN,INTI)
                          DERN =DERN*((FNT(NFTN,INTI+1)-FT1N)*DTIM+FT1N)

                      END IF !NFTN.NE.0

                     DFACTOR = TINV*DFLU(L,NNUD)*DERN

                     IPOS=IPNT_PAR(NP)
                     DERH(NNODE,IPOS,INEW) = DERH(NNODE,IPOS,INEW)
     &                                      - DFACTOR

                  END IF !NP.NE.0.AND.IBCOD(NNODE).EQ.1

              END IF !NZON.GT.0

      END DO !L=1,NUMEL


      END SUBROUTINE RHS_IN_CHP_TRLUM
