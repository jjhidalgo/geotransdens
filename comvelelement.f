      SUBROUTINE COMVEL_ELEMENT
     &(LD        ,IOCALCDEVF,IODIM     ,ISMAX   
     &,LMXNDL    ,NPPEL     ,NNUD      ,NUMEL    ,NUMNP    ,L       
     &,KXX       ,DPARELDH  ,GRDFF    ,HAUX1   ,DVDH_L      
     &,VEL       ,TRACT)


********************************************************************
C     PURPOSE
C     this sub calculates the velocity in an element when
C     there is no density dependent flow.
c     dqdh(j,l,i) = derivative of component j of q in element
c                   l to head in node i
********************************************************************


      IMPLICIT NONE

C------------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::LD,IOCALCDEVF,IODIM,LMXNDL,NNUD,NUMEL
     &          ,NUMNP,L,ISMAX,NPPEL
     &          
C------------------------- EXTERNAL VARIABLES: ARRAYS

      INTEGER*4::KXX(LMXNDL,NUMEL)

      REAL*8::DPARELDH(NPPEL,NUMEL),GRDFF(IODIM,LMXNDL,NUMEL)
     &       ,DVDH_L(LMXNDL,LD),VEL(IODIM),TRACT(9),HAUX1(NUMNP)
                  
C------------------------- INTERNAL VARIABLES: SCALARS
    
	INTEGER*4::I,IDIM,J,JNODE,KDIM

      REAL*8::SUM,DSUM,SUM2,DSUM2

C------------------------- INTERNAL VARIABLES: ARRAYS

      INTEGER*4::IND(3,3,3)



C------------------------- Array IND is used to simplify the coding of 
C------------------------- conductivity matrix times the gradient of FEM 
C------------------------- interpolation functions (shape functions)

      DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,2,6,5,6,3/
      DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/


      DVDH_L = 0D0 

C------------------------- step 1: loop over components of VEL

      DO IDIM=1,LD    

          SUM = 0D0
          DSUM = 0D0

C------------------------- step 2: loop over nodes of element

          DO JNODE=1,NNUD  

              SUM2 = 0D0 
	        DSUM2 = 0D0

C------------------------- step 3:loop over anisotropy directions

              DO KDIM=1,LD

                  IF (IND(IDIM,KDIM,LD).LE.ISMAX) THEN
                      SUM2 = SUM2
     &                     +TRACT(IND(IDIM,KDIM,LD))*GRDFF(KDIM,JNODE,L)
	            END IF

	            IF (IOCALCDEVF.EQ.1)THEN

	                IF (IND(IDIM,KDIM,LD).LE.ISMAX) THEN

                          DSUM2 = DSUM2 + DPARELDH(IND(IDIM,KDIM,LD),L)
     &                                    *GRDFF(KDIM,JNODE,L)
                      END IF

                  END IF !IOCALCDEVF.EQ.1

              END DO !KDIM=1,LD  

              SUM = SUM + SUM2*HAUX1(KXX(JNODE, L))
              DSUM = DSUM + DSUM2*HAUX1(KXX(JNODE, L))
	       
          END DO !JNODE=1,NNUD  

         VEL(IDIM) = -1D0*SUM


      IF (IOCALCDEVF.EQ.1) THEN

          DO JNODE=1,NNUD

                DVDH_L(JNODE,IDIM) = DVDH_L(JNODE,IDIM) - DSUM

              DO KDIM=1,LD

                 IF (IND(IDIM,KDIM,LD).LE.ISMAX) THEN

                    DVDH_L(JNODE,IDIM) = DVDH_L(JNODE,IDIM)
     &                    - TRACT(IND(IDIM,KDIM,LD))*GRDFF(KDIM,JNODE,L)
                   END IF

              END DO !KDIM=1,LD

           END DO !JNODE=1,NNUD

      END IF !IOCALCDEVF.EQ.1

      END DO !DO IDIM=1,LD    
      
      END SUBROUTINE COMVEL_ELEMENT
