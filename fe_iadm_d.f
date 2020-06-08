      SUBROUTINE FE_IADM_D
     ;(KXX,NUMNP,NUMEL,LMXNDL,IADD_D,IAD_D,IADN_D,MAXNB,LNNDEL)
	

!EXTERNALS
	INTEGER*4 NUMNP,NUMEL,LMXNDL
	INTEGER*4 KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),IADD_D(2*NUMNP)
     ;          ,IAD_D(2*MAXNB,2*NUMNP), IADN_D(2*NUMNP)


!INTERNALS
	INTEGER*4 L,I,J,INODE,JNODE,IROW,ICOL,K

C________________________Loops over elements and nodes in element
	DO L=1,NUMEL
		NNUD=LNNDEL(L)
		DO I=1,NNUD
		   INODE = KXX(I,L)
	       DO J=1,I
			  JNODE = KXX(J,L)


C_______________________Insert entry in all 4 blocks.
	          !I,J BLOCK 1
  			  IROW = (INODE-1)*2+1
			  ICOL = (JNODE-1)*2+1
                CALL INSERT (IROW,ICOL,K,IADN_D,IAD_D,2*MAXNB,2*NUMNP)

	          !J,I BLOCK 1
  			  IROW = (JNODE-1)*2+1
			  ICOL = (INODE-1)*2+1
                CALL INSERT (IROW,ICOL,K,IADN_D,IAD_D,2*MAXNB,2*NUMNP)

	          !I,J BLOCK 2
  			  IROW = (INODE-1)*2+1
			  ICOL = JNODE*2
                CALL INSERT (IROW,ICOL,K,IADN_D,IAD_D,2*MAXNB,2*NUMNP)

	          !J,I BLOCK 2
  			  IROW = (JNODE-1)*2+1
			  ICOL = INODE*2
                CALL INSERT (IROW,ICOL,K,IADN_D,IAD_D,2*MAXNB,2*NUMNP)

	          !I,J BLOCK 3
  			  IROW = INODE*2
			  ICOL = (JNODE-1)*2+1
                CALL INSERT (IROW,ICOL,K,IADN_D,IAD_D,2*MAXNB,2*NUMNP)

	          !J,I BLOCK 3
  			  IROW = JNODE*2
			  ICOL = (INODE-1)*2+1
                CALL INSERT (IROW,ICOL,K,IADN_D,IAD_D,2*MAXNB,2*NUMNP)
	
	          !I,J BLOCK 4
  			  IROW = INODE*2
			  ICOL = JNODE*2
                CALL INSERT (IROW,ICOL,K,IADN_D,IAD_D,2*MAXNB,2*NUMNP)

	          !J,I BLOCK 4
  			  IROW = JNODE*2
			  ICOL = INODE*2
                CALL INSERT (IROW,ICOL,K,IADN_D,IAD_D,2*MAXNB,2*NUMNP)
			  	      
		  ENDDO
	   ENDDO
	ENDDO
      

	DO I = 1,2*NUMNP
		DO J= 1, IADN_D(I)
			IF (I.EQ.IAD_D(J,I))THEN
				IADD_D(I)=J
	        ENDIF
	    ENDDO
	ENDDO

	RETURN
	END