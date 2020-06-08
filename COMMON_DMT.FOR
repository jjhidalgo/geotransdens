*      'COMMON_DMT.FOR'

      REAL*8 DELT_LST    ! Last DELTA_T 

      INTEGER  NUM_MTDZ,  ! Number of matrix diffusion zones
     .         IN_IV_DMT, ! Beginning of matrix diffusion data in IV (pointer)
     .         IN_RV_DMT, ! Beginning of matrix diffusion data in IR (pointer)
     .         I_ST_DMT,  ! Index of mat. diff. computation status
     .         IUNIN,  MAINF ! Input and output unit numbers
     ;        ,I_DTRA     ! Beginning of DTRA auxiliar (by nodes)
     ;        ,I_NELBYNUD ! Beginning of auxiliar array to store 
                          ! number of elements by node
     ;        ,I_CMP_NELBNUD  ! 0, number of elements by node not yet computed
                              ! 1, already computed. It is computed ONLY once.
*            I_ST_DMT=
*                        0- Nothing has been done yet
*                        1- Data has already been read (ENTDAT_MTDZ)
*                        2- Preliminary calculations done (COMP_PRE_MTDZ)
*                        3- MTDZ variables have been updated in response to 
*                           change in parameters (UPD_PAR_MTDZ)
*                        4- MTDZ variables have been updated in response to 
*                           change in time increment(UPD_DT_MTDZ)
*                        5- Contribution to transport included (COMP_CNTR_MTDZ)
*                        6- Contribution to derivatives included (COMP_DER_MTDZ)
*                        7- Derivatives of Io and In updated (UPD_DREC_MTDZ)
*                        8- Io and In updated (UPD_REC_MTDZ)

      COMMON /DMT/ DELT_LST, NUM_MTDZ, IN_IV_DMT, IN_RV_DMT, I_ST_DMT,
     .         IUNIN,  MAINF,  I_DTRA, I_NELBYNUD, I_CMP_NELBNUD