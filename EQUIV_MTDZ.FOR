
      INTEGER*4 IR_VOL,            ! Pointer to VOL_MZ
     .        IR_AN,             ! Pointer to AN_MZ
     .        IR_EXP,            ! Pointer to EXP_MZ
     .        IR_DEXP,           ! Pointer to DEXP_MZ
     .        IR_RHS,            ! Pointer to RHS_MZ
     .        IR_REC0,           ! Pointer to REC0
     .        IR_RECN,           ! Pointer to RECN
     .        IR_DREC0,          ! Pointer to DREC0
     .        IR_DRECN,          ! Pointer to DRECN
     .        II_NUD,            ! Pointer to NUD_MZ
     .        MXRPRV,            ! Max. number of real variables in COMMON/DMTZ/
     .        MXIPRV,            ! Max. number of integer variables in "    "
     .        MXEQIR,            ! Max. number of pointers to RV and IV
     .        IPRV_MZ

      REAL*8 RPRV_MZ

      PARAMETER (MXRPRV=20, ! Max. number of real variables in COMMON/DMTZ/
     .           MXIPRV=20, ! Max. number of integer variables in "    "
     .           MXEQIR=10) ! Max. number of pointers to RV and IV
                            !(IR_.,II_.. variables) saved

      DIMENSION IPRV_MZ(MXIPRV+MXEQIR), RPRV_MZ(MXRPRV)

      EQUIVALENCE (RPRV_MZ(1),ESP_MZ)
      EQUIVALENCE (IPRV_MZ(1),MTDZ_NUM)
      EQUIVALENCE (IPRV_MZ(MXIPRV+1),IR_VOL),
     .            (IPRV_MZ(MXIPRV+2),IR_AN),
     .            (IPRV_MZ(MXIPRV+3),IR_EXP),
     .            (IPRV_MZ(MXIPRV+4),IR_DEXP),
     .            (IPRV_MZ(MXIPRV+5),IR_RHS),
     .            (IPRV_MZ(MXIPRV+6),IR_REC0),
     .            (IPRV_MZ(MXIPRV+7),IR_RECN),
     .            (IPRV_MZ(MXIPRV+8),IR_DREC0),
     .            (IPRV_MZ(MXIPRV+9),IR_DRECN),
     .            (IPRV_MZ(MXIPRV+10),II_NUD)
