
DIM CARD A4.2 -> lee número de zonas de formation factor.
DIM A5.1 opción logarímica de FOF
DIM A9.2 pesos FOF

GRI B3 lee zona
NZONE_PAR(19)

PAR_WGT(19)= XLAMFOF

IOLG_PAR(19,1)= IOLGFOF
IOLG_PAR(19,2) > 0 si se estima alguna zona. se asigna en inver.f
INORPAR(22)=INFOF
LXPAREL(1,12,IPROB)
CFPAREL(L,11)
PAREL(L,16)
DPARELDC(16)
DPARELDH(16)

La posición en derc viene  dada por IPNT_PAR(I) 
donde I va de IVPAR(JJ,1),IVPAR(JJ,2) con
      NZFF = NZONE_PAR(19)
   IF (NZFF.GT.0) THEN
       DO IZON=1,NZFF
           JJ = INORPAR(22)+IZON 
           KONT=1 
           DO I=IVPAR(JJ,1),IVPAR(JJ,2) 
               DERIV(KONT)=WGT_PAR(I)*DERIVAUX
               IF (IVPAR(JJ,4).EQ.1) THEN 
                   DERIV(KONT)=DERIV(KONT)*PARAMC*DLOG(10.0D0)
               END IF
                   IPOS(KONT)=IPNT_PAR(I) 
                   KONT=KONT+1 
           ENDDO !I=IVPAR(JJ,1),IVPAR(JJ,2) 
               
           NPTOT=IVPAR(JJ,2)-IVPAR(JJ,1)+1 
       ENDDO !IZON
    ENDIF
        
CCCCCCESOT DEVERIA FUNCIONAR
      NZFOF = NZONE_PAR(19)
      IF (IOLG_PAR(19,2).GT.0 .AND. NZFOF.GT.0) THEN  !se estima alguna zona
          DO IZON=1,NZFOF
              JJ = INORPAR(22)+IZON
              IP = IVPAR(JJ,1)
              IF (IP.GT.0) THEN
                  DO I=IVPAR(JJ,1),IVPAR(JJ,2) 
                      DERIV=DERESDF(IZON)
                      IF (IVPAR(JJ,4).EQ.1) THEN 
                          DERIV=DERIV*PARAMC*DLOG(10.0D0)
                      ENDIF
                      IPOS=IPNT_PAR(I) 
                      VJAC(NO, IPOS) = DERIV 
                  ENDDO !I=IVPAR(JJ,1),IVPAR(JJ,2) 
              ENDIF !IP.GT.0
          ENDDO !IZON
      ENDIF ! IOLG_PAR(19,2).GT.0 .AND. NZFOF.GT.0
