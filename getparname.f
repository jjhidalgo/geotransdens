       CHARACTER*10 FUNCTION GETPARNAME
     ;(NSTAT,NTYPAR,PARNUMBER,IVPAR,INORPAR2,NZPAR,TYPENAME)
****************************************************************
* PURPOSE
* this function returns the name of the parameter that has the value 
* "parnumber" in ivpar.


      IMPLICIT NONE
      
      INTEGER*4 TYPES, I, NSTAT,NTYPAR
      INTEGER*4 PARNUMBER,NZPAR,IVPAR(NZPAR),INORPAR2(NTYPAR+1)
      CHARACTER*10 TYPENAME(NSTAT+NTYPAR)

      DO I=1,NZPAR
         IF (IVPAR(I) .EQ. PARNUMBER) THEN
            TYPES=0
            DO WHILE(INORPAR2(TYPES+1) .LT. I)
               TYPES=TYPES+1
            ENDDO
            GETPARNAME = TYPENAME(TYPES+NSTAT)
         ENDIF
      ENDDO

      END 
*******************************************************************************
      INTEGER*4 FUNCTION GETPARZNNR
     ;(PARNUMBER,IVPAR,INORPAR2,NZPAR,NTYPAR)
****************************************************************
* PURPOSE
* this function returns the zone of the parameter that has the value 
* "parnumber" in ivpar.

      IMPLICIT NONE
      
      INTEGER*4 PARNUMBER,NZPAR,NTYPAR,IVPAR(NZPAR),INORPAR2(NTYPAR+1) 
      INTEGER*4 TYPES, I
      
      DO I=1,NZPAR
         IF (IVPAR(I) .EQ. PARNUMBER) THEN
            TYPES=0
            DO WHILE(INORPAR2(TYPES+1) .LT. I)
               TYPES=TYPES+1
            ENDDO
            GETPARZNNR = I - INORPAR2(TYPES)
         ENDIF
      ENDDO
      END 
