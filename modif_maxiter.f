         SUBROUTINE MODIF_MAXITER (MAXITER)
        
         LOGICAL*4 FILEEX
         IF (FILEEX('DETENTE.DAT')) THEN      

           OPEN (UNIT=115,FILE='DETENTE.DAT',STATUS='OLD')
           REWIND(115)
           READ(115,*)IUPDATE
           IF (IUPDATE.EQ.1) READ(115,*)MAXITER
           CLOSE(115)

         END IF

         RETURN
         END
