       SUBROUTINE COMP_IND_CHANGE
     ;      (INALFC      ,IENTRY      ,INDCHANGES  ,INTI     ,NFTPAR
     ;      ,NZALF       ,NZPAR       ,TINC        ,TOLD     ,IOREGIMEN
     &      ,IODENS      ,IOFLLI)

C PURPOSE: Identify whether it is or not necessary to assemb the Right Hand Side 
C of the flow system at this time step. 
C ACRONYM: COMPute the INDicator of CHANGE in the RHS
C DESCRIPTION: There exist two situations which compel to compute the lineal RHS.
C I-the time step increment has changed, respect to the previous one
C II- The leakage coefficient has changed respect to the previous one.
C When verified I or II or both, this routine assign INDCHANGES to 1, which indicates
C that the RHS has to be computed. Otherwise, INDCHANGES=0 and the previous RHS
C can be used in the current time step solution process.
C EXTERNALS
C INTERNALS:
C HYSTORY: Programmed at november 1997 by G.Galarza

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION NFTPAR(NZPAR)     

       INDCHANGES=0                     !no changes, no new RHS computation
       IF(IOREGIMEN.EQ.1 .OR. IODENS.EQ.1 .OR. IOFLLI.NE.0) INDCHANGES=1
       IF ((TINC.NE.TOLD).OR.(IENTRY.EQ.1.AND.INTI.EQ.1)) THEN !the time step increment changed
         INDCHANGES=1
         RETURN
       ELSE
         DO NZ=1,NZALF
           NZP=INALFC+NZ
           NFT=NFTPAR(NZP)
           IF (NFT.NE.0) THEN   !any leakage coefficient changed
             INDCHANGES=1
             RETURN
           ENDIF
         END DO
       ENDIF

       RETURN                    
       END
