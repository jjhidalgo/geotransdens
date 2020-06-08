C   IMSL ROUTINE NAME   - LEQ1PB
C
C---------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - POSITIVE
C                           DEFINITE SYMMETRIC BAND MATRIX - BAND
C                           SYMMETRIC STORAGE MODE - SPACE ECONOMIZER
C                           SOLUTION
C
C   USAGE               - CALL LEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
C
C   ARGUMENTS    A      - THE COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B, WHERE A IS ASSUMED TO BE AN N X N
C                           POSITIVE DEFINITE BAND SYMMETRIC MATRIX. A
C                           IS STORED IN BAND SYMMETRIC STORAGE MODE
C                           AND THEREFORE HAS DIMENSION N BY (NC+1).
C                           (INPUT)
C                         ON OUTPUT, A IS REPLACED BY L WHERE
C                           A = L*L-TRANSPOSE. L IS A LOWER BAND
C                           MATRIX STORED IN BAND FORM AND THEREFORE    
C                           HAS DIMENSION N BY (NC+1). NOTE THAT THE    
C                           DIAGONAL ELEMENTS OF L ARE STORED IN        
C                           RECIPROCAL FORM. (OUTPUT)                   
C                N      - ORDER OF MATRIX A AND NUMBER OF ROWS IN B.    
C                           (INPUT)                                     
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS OF A.   
C                           (INPUT)                                     
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS          
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE 
C                           CALLING PROGRAM. (INPUT)                    
C                B      - INPUT MATRIX OF DIMENSION N X M CONTAINING    
C                           THE M RIGHT-HAND SIDES OF THE EQUATION      
C                           AX = B.                                     
C                         ON OUTPUT, THE N X M SOLUTION MATRIX X        
C                           REPLACES B.                                 
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS          
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE 
C                           CALLING PROGRAM. (INPUT)                    
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).    
C                           (INPUT)                                     
C                IDGT   - THE ELEMENTS OF A ARE ASSUMED TO BE CORRECT   
C                           TO IDGT DECIMAL DIGITS. (INPUT - CURRENTLY  
C                           NOT USED)                                   
C                D1     - COMPONENTS OF THE DETERMINANT OF A.           
C                D2         DETERMINANT(A) = D1*2.**D2. (OUTPUT)        
C                IER    - ERROR PARAMETER. (OUTPUT)                     
C                         TERMINAL ERROR                                
C                           IER = 129 INDICATES THAT THE MATRIX A IS    
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.    
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - LUDAPB,LUELPB,UERTST,UGETIO                   
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE LEQ1PB (A,N,NC,IA,B,IB,M,D1,D2,IER)
C                                                                       
      DOUBLE PRECISION   A(IA,1),B(IB,1),D1,D2                          
C                                  FIRST EXECUTABLE STATEMENT           
      IER = 0                                                           
C                                  DECOMPOSITION OF MATRIX A INTO       
C                                  L*L-TRANSPOSE                        
      CALL LUDAPB(A,N,NC,IA,A,IA,D1,D2,IER)                             
      IF (IER .NE. 0) GO TO 9000                                        
      DO 5 I = 1,M                                                      
C                                  SOLUTION OF AX = B                   
         CALL LUELPB(A,B(1,I),N,NC,IA,B(1,I))                           
    5 CONTINUE                                                          
      GO TO 9005                                                        
 9000 CONTINUE                                                          
      CALL UERTST (IER,'LEQ1PB')                                        
 9005 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - LEQT1B                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - JANUARY 1, 1978                               
C                                                                       
C   PURPOSE             - LINEAR EQUATION SOLUTION - BAND STORAGE       
C                           MODE - SPACE ECONOMIZER SOLUTION            
C                                                                       
C   USAGE               - CALL LEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,   
C                           IER)                                        
C                                                                       
C   ARGUMENTS    A      - INPUT/OUTPUT MATRIX OF DIMENSION N BY         
C                           (NUC+NLC+1). SEE PARAMETER IJOB.            
C                N      - ORDER OF MATRIX A AND THE NUMBER OF ROWS IN   
C                           B. (INPUT)                                  
C                NLC    - NUMBER OF LOWER CODIAGONALS IN MATRIX A.      
C                           (INPUT)                                     
C                NUC    - NUMBER OF UPPER CODIAGONALS IN MATRIX A.      
C                           (INPUT)                                     
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS          
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE 
C                           CALLING PROGRAM. (INPUT)                    
C                B      - INPUT/OUTPUT MATRIX OF DIMENSION N BY M.      
C                           ON INPUT, B CONTAINS THE M RIGHT-HAND SIDES 
C                           OF THE EQUATION AX = B. ON OUTPUT, THE      
C                           SOLUTION MATRIX X REPLACES B. IF IJOB = 1,  
C                           B IS NOT USED.                              
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).    
C                           (INPUT)                                     
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS          
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE 
C                           CALLING PROGRAM. (INPUT)                    
C                IJOB   - INPUT OPTION PARAMETER. IJOB = I IMPLIES WHEN 
C                           I = 0, FACTOR THE MATRIX A AND SOLVE THE    
C                             EQUATION AX = B. ON INPUT, A CONTAINS THE 
C                             COEFFICIENT MATRIX OF THE EQUATION AX = B,
C                             WHERE A IS ASSUMED TO BE AN N BY N BAND   
C                             MATRIX. A IS STORED IN BAND STORAGE MODE  
C                             AND THEREFORE HAS DIMENSION N BY          
C                             (NLC+NUC+1). ON OUTPUT, A IS REPLACED     
C                             BY THE U MATRIX OF THE L-U DECOMPOSITION  
C                             OF A ROWWISE PERMUTATION OF MATRIX A. U   
C                             IS STORED IN BAND STORAGE MODE.           
C                           I = 1, FACTOR THE MATRIX A. A CONTAINS THE  
C                             SAME INPUT/OUTPUT INFORMATION AS IF       
C                             IJOB = 0.                                 
C                           I = 2, SOLVE THE EQUATION AX = B. THIS      
C                             OPTION IMPLIES THAT LEQT1B HAS ALREADY    
C                             BEEN CALLED USING IJOB = 0 OR 1 SO THAT   
C                             THE MATRIX A HAS ALREADY BEEN FACTORED.   
C                             IN THIS CASE, OUTPUT MATRICES A AND XL    
C                             MUST HAVE BEEN SAVED FOR REUSE IN THE     
C                             CALL TO LEQT1B.                           
C                XL     - WORK AREA OF DIMENSION N*(NLC+1). THE FIRST   
C                           NLC*N LOCATIONS OF XL CONTAIN COMPONENTS OF 
C                           THE L MATRIX OF THE L-U DECOMPOSITION OF A  
C                           ROWWISE PERMUTATION OF A. THE LAST N        
C                           LOCATIONS CONTAIN THE PIVOT INDICES.        
C                IER    - ERROR PARAMETER. (OUTPUT)                     
C                         TERMINAL ERROR                                
C                           IER = 129 INDICATES THAT MATRIX A IS        
C                             ALGORITHMICALLY SINGULAR. (SEE THE        
C                             CHAPTER L PRELUDE).                       
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE LEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)             
C                                                                       
      DIMENSION          A(IA,1),XL(N,1),B(IB,1)                        
      DOUBLE PRECISION   A,XL,B,P,Q,ZERO,ONE,RN                         
      DATA               ZERO/0.D0/,ONE/1.0D0/                          
C                                  FIRST EXECUTABLE STATEMENT           
      IER = 0                                                           
      JBEG = NLC+1                                                      
      NLC1 = JBEG                                                       
      IF (IJOB .EQ. 2) GO TO 80                                         
      RN = N                                                            
C                                  RESTRUCTURE THE MATRIX               
C                                  FIND RECIPROCAL OF THE LARGEST       
C                                  ABSOLUTE VALUE IN ROW I              
      I = 1                                                             
      NC = JBEG+NUC                                                     
      NN = NC                                                           
      JEND = NC                                                         
      IF (N .EQ. 1 .OR. NLC .EQ. 0) GO TO 25                            
    5 K = 1                                                             
      P = ZERO                                                          
      DO 10 J = JBEG,JEND                                               
         A(I,K) = A(I,J)                                                
         Q =  DABS(A(I,K))                                              
         IF (Q .GT. P) P = Q                                            
         K = K+1                                                        
   10 CONTINUE                                                          
      IF (P .EQ. ZERO) GO TO 135                                        
      XL(I,NLC1) = ONE/P                                                
      IF (K .GT. NC) GO TO 20                                           
      DO 15 J = K,NC                                                    
         A(I,J) = ZERO                                                  
   15 CONTINUE                                                          
   20 I = I+1                                                           
      JBEG = JBEG-1                                                     
      IF (JEND-JBEG .EQ. N) JEND = JEND-1                               
      IF (I .LE. NLC) GO TO 5                                           
      JBEG = I                                                          
      NN = JEND                                                         
   25 JEND = N-NUC                                                      
      DO 40 I = JBEG,N                                                  
         P = ZERO                                                       
         DO 30 J = 1,NN                                                 
            Q =  DABS(A(I,J))                                           
            IF (Q .GT. P) P = Q                                         
   30    CONTINUE                                                       
         IF (P .EQ. ZERO) GO TO 135                                     
         XL(I,NLC1) = ONE/P                                             
         IF (I .EQ. JEND) GO TO 37                                      
         IF (I .LT. JEND) GO TO 40                                      
         K = NN+1                                                       
         DO 35 J = K,NC                                                 
            A(I,J) = ZERO                                               
   35    CONTINUE                                                       
   37    NN = NN-1                                                      
   40 CONTINUE                                                          
      L = NLC                                                           
C                                  L-U DECOMPOSITION                    
      DO 75 K = 1,N                                                     
         P =  DABS(A(K,1))*XL(K,NLC1)                                   
         I = K                                                          
         IF (L .LT. N) L = L+1                                          
         K1 = K+1                                                       
         IF (K1 .GT. L) GO TO 50                                        
         DO 45 J = K1,L                                                 
            Q = DABS(A(J,1))*XL(J,NLC1)                                 
            IF (Q .LE. P) GO TO 45                                      
            P = Q                                                       
            I = J                                                       
   45    CONTINUE                                                       
   50    XL(I,NLC1) = XL(K,NLC1)                                        
         XL(K,NLC1) = I                                                 
C                                  SINGULARITY FOUND                    
         Q = RN+P                                                       
         IF (Q .EQ. RN) GO TO 135                                       
C                                  INTERCHANGE ROWS I AND K             
         IF (K .EQ. I) GO TO 60                                         
         DO 55 J = 1,NC                                                 
            P = A(K,J)                                                  
            A(K,J) = A(I,J)                                             
            A(I,J) = P                                                  
   55    CONTINUE                                                       
   60    IF (K1 .GT. L) GO TO 75                                        
         DO 70 I = K1,L                                                 
            P = A(I,1)/A(K,1)                                           
            IK = I-K                                                    
            XL(K1,IK) = P                                               
            DO 65 J = 2,NC                                              
               A(I,J-1) = A(I,J)-P*A(K,J)                               
   65    CONTINUE                                                       
         A(I,NC) = ZERO                                                 
   70    CONTINUE                                                       
   75 CONTINUE                                                          
      IF (IJOB .EQ. 1) GO TO 9005                                       
C                                  FORWARD SUBSTITUTION                 
   80 L = NLC                                                           
      DO 105 K = 1,N                                                    
         I = XL(K,NLC1)                                                 
         IF (I .EQ. K) GO TO 90                                         
         DO 85 J = 1,M                                                  
            P = B(K,J)                                                  
            B(K,J) = B(I,J)                                             
            B(I,J) = P                                                  
   85    CONTINUE                                                       
   90    IF (L .LT. N) L = L+1                                          
         K1 = K+1                                                       
         IF (K1 .GT. L) GO TO 105                                       
         DO 100 I = K1,L                                                
            IK = I-K                                                    
            P = XL(K1,IK)                                               
            DO 95 J = 1,M                                               
               B(I,J) = B(I,J)-P*B(K,J)                                 
   95       CONTINUE                                                    
  100    CONTINUE                                                       
  105 CONTINUE                                                          
C                                  BACKWARD SUBSTITUTION                
      JBEG = NUC+NLC                                                    
      DO 125 J = 1,M                                                    
         L = 1                                                          
         K1 = N+1                                                       
         DO 120 I = 1,N                                                 
            K = K1-I                                                    
            P = B(K,J)                                                  
            IF (L .EQ. 1) GO TO 115                                     
            DO 110 KK = 2,L                                             
               IK = KK+K                                                
               P = P-A(K,KK)*B(IK-1,J)                                  
  110       CONTINUE                                                    
  115       B(K,J) = P/A(K,1)                                           
            IF (L .LE. JBEG) L = L+1                                    
  120    CONTINUE                                                       
  125 CONTINUE                                                          
      GO TO 9005                                                        
  135 IER = 129                                                         
      CONTINUE
      CALL UERTST(IER,'LEQT1B')                                         
 9005 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - LEQT1P                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - JANUARY 1, 1978                               
C                                                                       
C   PURPOSE             - LINEAR EQUATION SOLUTION - POSITIVE DEFINITE  
C                           MATRIX - SYMMETRIC STORAGE MODE - SPACE     
C                           ECONOMIZER SOLUTION                         
C                                                                       
C   USAGE               - CALL LEQT1P (A,M,N,B,IB,IDGT,D1,D2,IER)       
C                                                                       
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING THE
C                           N BY N COEFFICIENT MATRIX OF THE EQUATION   
C                           AX = B. A IS A POSITIVE DEFINITE SYMMETRIC  
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE.    
C                         ON OUTPUT, A IS REPLACED BY THE LOWER         
C                           TRIANGULAR MATRIX L WHERE A = L*L-TRANSPOSE.
C                           L IS STORED IN SYMMETRIC STORAGE MODE WITH  
C                           THE DIAGONAL ELEMENTS OF L IN RECIPROCAL    
C                           FORM.                                       
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).    
C                           (INPUT)                                     
C                N      - ORDER OF A AND NUMBER OF ROWS IN B. (INPUT)   
C                B      - INPUT MATRIX OF DIMENSION N BY M CONTAINING   
C                           THE RIGHT-HAND SIDES OF THE EQUATION        
C                           AX = B.                                     
C                         ON OUTPUT, THE N BY M SOLUTION MATRIX X       
C                           REPLACES B.                                 
C                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN    
C                           THE DIMENSION STATEMENT IN THE CALLING      
C                           PROGRAM. (INPUT)                            
C                IDGT   - THE ELEMENTS OF A ARE ASSUMED TO BE CORRECT   
C                           TO IDGT DECIMAL DIGITS. (CURRENTLY NOT USED)
C                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.           
C                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)        
C                IER    - ERROR PARAMETER. (OUTPUT)                     
C                         TERMINAL ERROR                                
C                           IER = 129 INDICATES THAT THE INPUT MATRIX   
C                             A IS ALGORITHMICALLY NOT POSITIVE         
C                             DEFINITE. (SEE THE CHAPTER L PRELUDE).    
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - LUDECP,LUELMP,UERTST,UGETIO                   
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE LEQT1P (A,M,N,B,IB,D1,D2,IER)
C                                                                       
      DIMENSION          A(1),B(IB,1)                                   
      DOUBLE PRECISION   A,B,D1,D2                                      
C                                  FIRST EXECUTABLE STATEMENT           
C                                  INITIALIZE IER                       
      IER = 0                                                           
C                                  DECOMPOSE A                          
      CALL LUDECP (A,A,N,D1,D2,IER)                                     
      IF (IER.NE.0) GO TO 9000                                          
C                                  PERFORM ELIMINATION                  
      DO 5 I = 1,M                                                      
         CALL LUELMP (A,B(1,I),N,B(1,I))                                
    5 CONTINUE                                                          
      GO TO 9005                                                        
 9000 CONTINUE                                                          
      CALL UERTST(IER,'LEQT1P')                                         
 9005 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - LUDAPB                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - JANUARY 1, 1978                               
C                                                                       
C   PURPOSE             - DECOMPOSITION OF A POSITIVE DEFINITE BAND     
C                           SYMMETRIC MATRIX - BAND SYMMETRIC STORAGE   
C                           MODE                                        
C                                                                       
C   USAGE               - CALL LUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)       
C                                                                       
C   ARGUMENTS    A      - N BY N POSITIVE DEFINITE BAND SYMMETRIC       
C                           MATRIX STORED IN BAND SYMMETRIC STORAGE     
C                           MODE. A SHOULD BE DIMENSIONED AT LEAST      
C                           N BY NC+1. (INPUT)                          
C                N      - ORDER OF A. (INPUT)                           
C                NC     - NUMBER OF UPPER OR LOWER CODIAGONALS OF A.    
C                           (INPUT)                                     
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS          
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE 
C                           CALLING PROGRAM. (INPUT)                    
C                UL     - OUTPUT MATRIX L WHERE A = L*L-TRANSPOSE. L IS 
C                           STORED IN BAND STORAGE MODE. UL SHOULD BE   
C                           DIMENSIONED AT LEAST N BY NC+1. NOTE - THE  
C                           DIAGONAL OF UL CONTAINS THE RECIPROCALS     
C                           OF THE ACTUAL DIAGONAL ELEMENTS.            
C                IU     - ROW DIMENSION OF MATRIX UL EXACTLY AS         
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE 
C                           CALLING PROGRAM. (INPUT)                    
C                D1     - COMPONENTS OF THE DETERMINANT OF A.           
C                D2         DETERMINANT(A) = D1*2.**D2. (OUTPUT)        
C                IER    - ERROR PARAMETER. (OUTPUT)                     
C                         TERMINAL ERROR                                
C                           N = 129 INDICATES THAT THE MATRIX A IS      
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.    
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE LUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)                     
C                                                                       
      DOUBLE PRECISION   ZERO,A(IA,1),UL(IU,1),D1,D2,ONE,SUM,
     *                   RN,FOUR,SIXTN,SIXTH                            
      DATA               ZERO/0.0D0/,FOUR/4.0D0/,SIXTN/16.D0/,          
     *                   SIXTH/.0625D0/,ONE/1.0D0/                      
C                                  FIRST EXECUTABLE STATEMENT           
      IER = 0                                                           
      RN = ONE/(N*SIXTN)                                                
      D1 = ONE                                                          
      D2 = ZERO                                                         
      NCP1 = NC+1                                                       
      IF (NC .EQ. 0) GO TO 15                                           
C                                  INITIALIZE ZERO ELEMENTS             
      DO 10 I = 1,NC                                                    
         DO 5 J = I,NC                                                  
            K = NCP1-J                                                  
            UL(I,K) = ZERO                                              
    5    CONTINUE                                                       
   10 CONTINUE                                                          
C                                  I IS ROW INDEX OF ELEMENT BEING      
C                                  COMPUTED                             
   15 DO 60 I = 1,N                                                     
         IMNCP1 = I-NCP1                                                
         I1 = MAX0(1,1-IMNCP1)                                          
C                                  J IS COLUMN INDEX OF ELEMENT BEING   
C                                  COMPUTED                             
         DO 60 J = I1,NCP1                                              
C                                  L IS ROW INDEX OF PREVIOUSLY COMPUTED
C                                  VECTOR BEING USED TO COMPUTE INNER   
C                                  PRODUCT                              
            L = IMNCP1+J                                                
            I2 = NCP1-J                                                 
            SUM = A(I,J)                                                
            JM1 = J-1                                                   
            IF (JM1) 30,30,20                                           
   20       DO 25 K = 1,JM1                                             
C                                  M IS COLUMN INDEX                    
               M = I2+K                                                 
               SUM = SUM-UL(I,K)*UL(L,M)                                
   25       CONTINUE                                                    
   30       IF (J .NE. NCP1) GO TO 55                                   
            IF(A(I,J)+SUM*RN .LE. A(I,J))GO TO 65                       
            UL(I,J) = ONE/DSQRT(SUM)                                    
C                                  UPDATE THE DETERMINANT               
            D1 = D1*SUM                                                 
   35       IF (DABS(D1)-ONE) 45,45,40                                  
   40       D1 = D1*SIXTH                                               
            D2 = D2+FOUR                                                
            GO TO 35                                                    
   45       IF (DABS(D1)-SIXTH) 50,50,60                                
   50       D1 = D1*SIXTN                                               
            D2 = D2-FOUR                                                
            GO TO 45                                                    
   55       UL(I,J) = SUM*UL(L,NCP1)                                    
   60 CONTINUE                                                          
      GO TO 9005                                                        
   65 IER = 129                                                         

C_____________________________________ Writes "dangerous" nodes in case 

      WRITE(6,*) ' ERROR ROUTINE LUDAPB. CHECK RES.OUT I=',I
      WRITE(25,*) ' ERROR ROUTINE LUDAPB. I=',I

      CONTINUE
      CALL UERTST(IER,'LUDAPB')                                         
 9005 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - LUDECP                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - JANUARY 1, 1978                               
C                                                                       
C   PURPOSE             - DECOMPOSITION OF A POSITIVE DEFINITE MATRIX - 
C                           SYMMETRIC STORAGE MODE                      
C                                                                       
C   USAGE               - CALL LUDECP (A,UL,N,D1,D2,IER)                
C                                                                       
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING    
C                           THE N BY N POSITIVE DEFINITE SYMMETRIC      
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE.    
C                UL     - OUTPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING   
C                           THE DECOMPOSED MATRIX L SUCH THAT A = L*    
C                           L-TRANSPOSE. L IS STORED IN SYMMETRIC       
C                           STORAGE MODE. THE DIAGONAL OF L CONTAINS THE
C                           RECIPROCALS OF THE ACTUAL DIAGONAL ELEMENTS.
C                N      - ORDER OF A. (INPUT)                           
C                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.           
C                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)        
C                IER    - ERROR PARAMETER. (OUTPUT)                     
C                         TERMINAL ERROR                                
C                           IER = 129 INDICATES THAT MATRIX A IS        
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.    
C                             (SEE THE CHAPTER L PRELUDE).              
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE LUDECP (A,UL,N,D1,D2,IER)                              
C                                                                       
      DIMENSION          A(1),UL(1)                                     
      DOUBLE PRECISION   A,UL,D1,D2,ZERO,ONE,FOUR,SIXTN,SIXTH,X,RN      
      DATA               ZERO,ONE,FOUR,SIXTN,SIXTH/                     
     *                   0.0D0,1.D0,4.D0,16.D0,.0625D0/                 
C                                  FIRST EXECUTABLE STATEMENT           
      D1=ONE                                                            
      D2=ZERO                                                           
      RN = ONE/(N*SIXTN)                                                
      IP = 1                                                            
      IER=0                                                             
      DO 45 I = 1,N                                                     
         IQ = IP                                                        
         IR = 1                                                         
         DO 40 J = 1,I                                                  
            X = A(IP)                                                   
            IF (J .EQ. 1) GO TO 10                                      
            DO 5  K=IQ,IP1                                              
               X = X - UL(K) * UL(IR)                                   
               IR = IR+1                                                
    5       CONTINUE                                                    
   10       IF (I.NE.J) GO TO 30                                        
            D1 = D1*X                                                   
            IF (A(IP) + X*RN .LE. A(IP)) GO TO 50                       
   15       IF (DABS(D1).LE.ONE) GO TO 20                               
            D1 = D1 * SIXTH                                             
            D2 = D2 + FOUR                                              
            GO TO 15                                                    
   20       IF (DABS(D1) .GE. SIXTH) GO TO 25                           
            D1 = D1 * SIXTN                                             
            D2 = D2 - FOUR                                              
            GO TO 20                                                    
   25       UL(IP) = ONE/DSQRT(X)                                       
            GO TO 35                                                    
   30       UL(IP) = X * UL(IR)                                         
   35       IP1 = IP                                                    
            IP = IP+1                                                   
            IR = IR+1                                                   
   40    CONTINUE                                                       
   45 CONTINUE                                                          
      GO TO 9005                                                        
   50 IER = 129                                                         
      CONTINUE
      CALL UERTST(IER,'LUDECP')                                         
 9005 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - LUELMP                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - JANUARY 1, 1978                               
C                                                                       
C   PURPOSE             - ELIMINATION PART OF THE SOLUTION OF AX=B -    
C                           POSITIVE DEFINITE MATRIX - SYMMETRIC        
C                           STORAGE MODE                                
C                                                                       
C   USAGE               - CALL LUELMP (A,B,N,X)                         
C                                                                       
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING    
C                           THE N BY N MATRIX L WHERE A = L*L-TRANSPOSE.
C                           L IS A LOWER TRIANGULAR MATRIX STORED IN    
C                           SYMMETRIC STORAGE MODE. THE MAIN DIAGONAL   
C                           ELEMENTS OF L ARE STORED IN RECIPROCAL      
C                           FORM. MATRIX L MAY BE OBTAINED FROM IMSL    
C                           ROUTINE LUDECP.                             
C                B      - VECTOR OF LENGTH N CONTAINING THE RIGHT HAND  
C                           SIDE OF THE EQUATION AX = B. (INPUT)        
C                N      - ORDER OF A AND THE LENGTH OF B AND X. (INPUT) 
C                X      - VECTOR OF LENGTH N CONTAINING THE SOLUTION TO 
C                           THE EQUATION AX = B. (OUTPUT)               
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE LUELMP (A,B,N,X)                                       
C                                                                       
      DIMENSION          A(1),B(1),X(1)                                 
      DOUBLE PRECISION   A,B,X,T,ZERO                                   
      DATA               ZERO/0.0D0/                                    
C                                  FIRST EXECUTABLE STATEMENT           
C                                  SOLUTION OF LY = B                   
      IP=1                                                              
      IW = 0                                                            
      DO 15 I=1,N                                                       
         T=B(I)                                                         
         IM1 = I-1                                                      
         IF (IW .EQ. 0) GO TO 9                                         
         IP=IP+IW-1                                                     
         DO 5 K=IW,IM1                                                  
            T = T-A(IP)*X(K)                                            
            IP=IP+1                                                     
    5    CONTINUE                                                       
         GO TO 10                                                       
    9    IF (T .NE. ZERO) IW = I                                        
         IP = IP+IM1                                                    
   10    X(I)=T*A(IP)                                                   
         IP=IP+1                                                        
   15 CONTINUE                                                          
C                                  SOLUTION OF UX = Y                   
      N1 = N+1                                                          
      DO 30 I = 1,N                                                     
         II = N1-I                                                      
         IP=IP-1                                                        
         IS=IP                                                          
         IQ=II+1                                                        
         T=X(II)                                                        
         IF (N.LT.IQ) GO TO 25                                          
         KK = N                                                         
         DO 20 K=IQ,N                                                   
            T = T - A(IS) * X(KK)                                       
            KK = KK-1                                                   
            IS = IS-KK                                                  
   20    CONTINUE                                                       
   25    X(II)=T*A(IS)                                                  
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - LUELPB                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - JANUARY 1, 1978                               
C                                                                       
C   PURPOSE             - ELIMINATION PART OF SOLUTION OF AX=B -        
C                           POSITIVE DEFINITE BAND SYMMETRIC MATRIX -   
C                           BAND SYMMETRIC STORAGE MODE                 
C                                                                       
C   USAGE               - CALL LUELPB (UL,B,N,NC,IA,X)                  
C                                                                       
C   ARGUMENTS    UL     - THE RESULT L COMPUTED IN THE ROUTINE LUDAPB   
C                           WHERE A = L*L-TRANSPOSE. L IS A LOWER BAND  
C                           MATRIX STORED IN BAND STORAGE MODE AND      
C                           THEREFORE HAS DIMENSION N X (NC+1). THE     
C                           MAIN DIAGONAL ELEMENTS OF L ARE STORED IN   
C                           RECIPROCAL FORM. (INPUT)                    
C                B      - VECTOR OF LENGTH N CONTAINING THE RIGHT HAND  
C                           SIDE OF THE EQUATION AX = B. (INPUT)        
C                N      - ORDER OF A AND THE LENGTH OF B AND X. (INPUT) 
C                NC     - NUMBER OF LOWER CODIAGONALS OF A. (INPUT)     
C                IA     - ROW DIMENSION OF MATRIX UL EXACTLY AS         
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE 
C                           CALLING PROGRAM. (INPUT)                    
C                X      - VECTOR OF LENGTH N CONTAINING THE SOLUTION TO 
C                           THE EQUATION AX = B. (OUTPUT)               
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE LUELPB (UL,B,N,NC,IA,X)                                
C                                                                       
      DOUBLE PRECISION   UL(IA,1),B(1),X(1),ZERO,SUM                    
      DATA               ZERO/0.0D0/                                    
C                                  FIRST EXECUTABLE STATEMENT           
C                                  SOLUTION LY = B                      
      NC1 = NC+1                                                        
      IW = 0                                                            
      L = 0                                                             
      DO 15 I = 1,N                                                     
         SUM = B(I)                                                     
         IF (NC .LE. 0) GO TO 10                                        
         IF (IW .EQ. 0) GO TO 9                                         
         L = L+1                                                        
         IF (L .GT. NC) L = NC                                          
         K = NC1-L                                                      
         KL = I-L                                                       
         DO 5 J = K,NC                                                  
            SUM = SUM -X(KL) * UL(I,J)                                  
            KL = KL+1                                                   
    5    CONTINUE                                                       
         GO TO 10                                                       
    9    IF (SUM .NE. ZERO) IW = 1                                      
   10    X(I) = SUM*UL(I,NC1)                                           
   15 CONTINUE                                                          
C                                  SOLUTION UX = Y                      
      X(N) = X(N)*UL(N,NC1)
      IF (N .LE. 1) GO TO 40                                            
      N1 = N+1                                                          
      DO 35 I = 2,N                                                     
         K = N1-I                                                       
         SUM = X(K)                                                     
         IF (NC .LE. 0) GO TO 30                                        
         KL = K+1                                                       
         K1 = MIN0(N,K+NC)                                              
         L = 1                                                          
         DO 25 J = KL,K1                                                
            SUM = SUM -X(J) * UL(J,NC1-L)                               
            L = L+1                                                     
   25    CONTINUE                                                       
   30    X(K) = SUM*UL(K,NC1)                                           
   35 CONTINUE                                                          
   40 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - UERTST                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/SINGLE                                    
C                                                                       
C   LATEST REVISION     - JUNE 1, 1982                                  
C                                                                       
C   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION 
C                                                                       
C   USAGE               - CALL UERTST (IER,NAME)                        
C                                                                       
C   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)                      
C                           IER = I+J WHERE                             
C                             I = 128 IMPLIES TERMINAL ERROR MESSAGE,   
C                             I =  64 IMPLIES WARNING WITH FIX MESSAGE, 
C                             I =  32 IMPLIES WARNING MESSAGE.          
C                             J = ERROR CODE RELEVANT TO CALLING        
C                                 ROUTINE.                              
C                NAME   - A CHARACTER STRING OF LENGTH SIX PROVIDING    
C                           THE NAME OF THE CALLING ROUTINE. (INPUT)    
C                                                                       
C   PRECISION/HARDWARE  - SINGLE/ALL                                    
C                                                                       
C   REQD. IMSL ROUTINES - UGETIO,USPKD                                  
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN        
C                TO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT           
C                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS          
C                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).                   
C                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING       
C                UGETIO AS FOLLOWS..                                    
C                                NIN = 0                                
C                                NOUT = NEW OUTPUT UNIT NUMBER          
C                                CALL UGETIO(3,NIN,NOUT)                
C                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.              
C                                                                       
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE UERTST (IER,NAME)                                      
C                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER            IER                                            
      CHARACTER*6        NAME
*AMS      INTEGER            NAME(1)                                    0530
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER            IEQDF,IOUNIT,LEVEL,LEVOLD,NIN,NMTB             
      CHARACTER*6        NAMEQ,NAMSET,NAMUPK,IEQ                        
      DATA               NAMSET/'UERSET'/                               
      DATA               NAMEQ/'      '/,IEQ/'======'/                  
      DATA               LEVEL/4/,IEQDF/0/                              
C                                  UNPACK NAME INTO NAMUPK              
C                                  FIRST EXECUTABLE STATEMENT           
      CALL USPKD (NAME,6,NAMUPK,NMTB)                                   
C                                  GET OUTPUT UNIT NUMBER               
      CALL UGETIO(1,NIN,IOUNIT)                                         
C                                  CHECK IER                            
      IF (IER.GT.999) GO TO 25                                          
      IF (IER.LT.-32) GO TO 55                                          
      IF (IER.LE.128) GO TO 5                                           
      IF (LEVEL.LT.1) GO TO 30                                          
C                                  PRINT TERMINAL MESSAGE               
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAMUPK             
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAMUPK                       
      GO TO 30                                                          
    5 IF (IER.LE.64) GO TO 10                                           
      IF (LEVEL.LT.2) GO TO 30                                          
C                                  PRINT WARNING WITH FIX MESSAGE       
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAMUPK             
      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAMUPK                       
      GO TO 30                                                          
   10 IF (IER.LE.32) GO TO 15                                           
C                                  PRINT WARNING MESSAGE                
      IF (LEVEL.LT.3) GO TO 30                                          
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAMUPK             
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAMUPK                       
      GO TO 30                                                          
   15 CONTINUE                                                          
C                                  CHECK FOR UERSET CALL                
*AMS      DO 20 I=1,6                                                   0880
         IF (NAMUPK.NE.NAMSET) GO TO 25                                 
*AMS   20 CONTINUE                                                      0900
      LEVOLD = LEVEL                                                    
      LEVEL = IER                                                       
      IER = LEVOLD                                                      
      IF (LEVEL.LT.0) LEVEL = 4                                         
      IF (LEVEL.GT.4) LEVEL = 4                                         
      GO TO 30                                                          
   25 CONTINUE                                                          
      IF (LEVEL.LT.4) GO TO 30                                          
C                                  PRINT NON-DEFINED MESSAGE            
      IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAMUPK             
      IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAMUPK                       
   30 IEQDF = 0                                                         
      RETURN                                                            
   35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,                   
     1       20H) FROM IMSL ROUTINE ,A6,A1,A6)                          
   40 FORMAT(27H *** WARNING WITH FIX ERROR,2X,7H(IER = ,I3,            
     1       20H) FROM IMSL ROUTINE ,A6,A1,A6)                          
   45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,                    
     1       20H) FROM IMSL ROUTINE ,A6,A1,A6)                          
   50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,                   
     1       20H) FROM IMSL ROUTINE ,A6,A1,A6)                          
C                                                                       
C                                  SAVE P FOR P = R CASE                
C                                    P IS THE PAGE NAMUPK               
C                                    R IS THE ROUTINE NAMUPK            
   55 IEQDF = 1                                                         
*AMS      DO 60 I=1,6                                                   1170
*AMS   60 NAMEQ(I) = NAMUPK(I)                                          1180
      NAMEQ = NAMUPK
      RETURN
      END                                                               
C   IMSL ROUTINE NAME   - UGETIO                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/SINGLE                                    
C                                                                       
C   LATEST REVISION     - JUNE 1, 1981                                  
C                                                                       
C   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW     
C                           VALUES FOR INPUT AND OUTPUT UNIT            
C                           IDENTIFIERS.                                
C                                                                       
C   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)                    
C                                                                       
C   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)                     
C                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT     
C                           UNIT IDENTIFIER VALUES ARE RETURNED IN NIN  
C                           AND NOUT, RESPECTIVELY.                     
C                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS     
C                           RESET FOR SUBSEQUENT USE.                   
C                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS    
C                           RESET FOR SUBSEQUENT USE.                   
C                NIN    - INPUT UNIT IDENTIFIER.                        
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.          
C                NOUT   - OUTPUT UNIT IDENTIFIER.                       
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.          
C                                                                       
C   PRECISION/HARDWARE  - SINGLE/ALL                                    
C                                                                       
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT    
C                OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT     
C                IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR  
C                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.    
C                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS. 
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE UGETIO(IOPT,NIN,NOUT)                                  
C                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER            IOPT,NIN,NOUT                                  
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER            NIND,NOUTD                                     
      DATA               NIND/1/,NOUTD/2/                               
C                                  FIRST EXECUTABLE STATEMENT           
      IF (IOPT.EQ.3) GO TO 10                                           
      IF (IOPT.EQ.2) GO TO 5                                            
      IF (IOPT.NE.1) GO TO 9005                                         
      NIN = NIND                                                        
      NOUT = NOUTD                                                      
      GO TO 9005                                                        
    5 NIND = NIN                                                        
      GO TO 9005                                                        
   10 NOUTD = NOUT                                                      
 9005 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - USPKD                                         
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/SINGLE                                    
C                                                                       
C   LATEST REVISION     - NOVEMBER 1, 1984                              
C                                                                       
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES THAT HAVE     
C                           CHARACTER STRING ARGUMENTS                  
C                                                                       
C   USAGE               - CALL USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)     
C                                                                       
C   ARGUMENTS    PACKED - CHARACTER STRING TO BE UNPACKED.(INPUT)       
C                NCHARS - LENGTH OF PACKED. (INPUT)  SEE REMARKS.       
C                UNPAKD - INTEGER ARRAY TO RECEIVE THE UNPACKED         
C                         REPRESENTATION OF THE STRING. (OUTPUT)        
C                NCHMTB - NCHARS MINUS TRAILING BLANKS. (OUTPUT)        
C                                                                       
C   PRECISION/HARDWARE  - SINGLE/ALL                                    
C                                                                       
C   REQD. IMSL ROUTINES - NONE                                          
C                                                                       
C   REMARKS  1.  USPKD UNPACKS A CHARACTER STRING INTO AN INTEGER ARRAY 
C                IN (A1) FORMAT.                                        
C            2.  UP TO 129 CHARACTERS MAY BE USED.  ANY IN EXCESS OF    
C                THAT ARE IGNORED.                                      
C                                                                       
C   COPYRIGHT           - 1984 BY IMSL, INC.  ALL RIGHTS RESERVED.      
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE.  NO OTHER WARRANTY,   
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
      SUBROUTINE USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)                   
C                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER            NC,NCHARS,NCHMTB                               
C                                                                       
      CHARACTER*6 PACKED                       ! LINEA ANYADIDA
      CHARACTER*6         UNPAKD                                        
*AMS      DATA               IBLANK /1H /                               0420
C                                  INITIALIZE NCHMTB                    
      NCHMTB = 0                                                        
C                                  RETURN IF NCHARS IS LE ZERO          
      IF(NCHARS.LE.0) RETURN                                            
*AMS      write(6,*) ' hola, packed es ',packed//'+fin'
C                                  SET NC=NUMBER OF CHARS TO BE DECODED 
      NC = MIN0 (129,NCHARS)                                            
*AMS      READ(PACKED,150) (UNPAKD(I),I=1,NC)     ! CAMBIO EFECTUADO    0490
        UNPAKD=PACKED

**** DECODE (NC,150,PACKED) (UNPAKD(I),I=1,NC)   ! MODIFICACION         

*AMS  150 FORMAT (129A1)
C                                  CHECK UNPAKD ARRAY AND SET NCHMTB    
C                                  BASED ON TRAILING BLANKS FOUND       
      DO 200 N = 1,NC                                                   
         NN = NC - N + 1                                                
         IF (UNPAKD(NN:NN).NE.' ') GO TO 210
  200 CONTINUE                                                          
      NN = 0                                                            
  210 NCHMTB = NN                                                       
      RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - EHOBKS                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - JUNE 1, 1982                                  
C                                                                       
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRS     
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE EHOBKS (A,N,M1,M2,Z,IZ)                                
C                                                                       
      DIMENSION          A(1),Z(IZ,1)                                   
      DOUBLE PRECISION   A,Z,H,S                                        
C                                  FIRST EXECUTABLE STATEMENT           
      IF (N .EQ. 1) GO TO 30                                            
      DO 25 I=2,N                                                       
         L = I-1                                                        
         IA = (I*L)/2                                                   
         H = A(IA+I)                                                    
         IF (H.EQ.0.D0) GO TO 25                                        
C                                  DERIVES EIGENVECTORS M1 TO M2 OF     
C                                  THE ORIGINAL MATRIX FROM EIGENVECTORS
C                                  M1 TO M2 OF THE SYMMETRIC            
C                                  TRIDIAGONAL MATRIX                   
         DO 20 J = M1,M2                                                
            S = 0.0D0                                                   
            DO 10 K = 1,L                                               
               S = S+A(IA+K)*Z(K,J)                                     
   10       CONTINUE                                                    
            S = S/H                                                     
            DO 15 K=1,L                                                 
               Z(K,J) = Z(K,J)-S*A(IA+K)                                
   15       CONTINUE                                                    
   20    CONTINUE                                                       
   25 CONTINUE                                                          
   30 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - EHOUSS                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - NOVEMBER 1, 1984                              
C                                                                       
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRS     
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE EHOUSS (A,N,D,E,E2)                                    
C                                                                       
      DIMENSION          A(1),D(N),E(N),E2(N)                           
      DOUBLE PRECISION   A,D,E,E2,ZERO,H,SCALE,F,G,HH                   
      DATA               ZERO/0.0D0/                                    
C                                  FIRST EXECUTABLE STATEMENT           
      NP1 = N+1                                                         
      NN = (N*NP1)/2-1                                                  
      NBEG = NN+1-N                                                     
      DO 70 II = 1,N                                                    
         I = NP1-II                                                     
         L = I-1                                                        
         H = ZERO                                                       
         SCALE = ZERO                                                   
         IF (L .LT. 1) GO TO 10                                         
C                                  SCALE ROW (ALGOL TOL THEN NOT NEEDED)
         NK = NN                                                        
         DO 5 K = 1,L                                                   
            SCALE = SCALE+DABS(A(NK))                                   
            NK = NK-1                                                   
    5    CONTINUE                                                       
         IF (SCALE .NE. ZERO) GO TO 15                                  
   10    E(I) = ZERO                                                    
         E2(I) = ZERO                                                   
         GO TO 65                                                       
   15    NK = NN                                                        
         DO 20 K = 1,L                                                  
            A(NK) = A(NK)/SCALE                                         
            H = H+A(NK)*A(NK)                                           
            NK = NK-1                                                   
   20    CONTINUE                                                       
         E2(I) = SCALE*SCALE*H                                          
         F = A(NN)                                                      
         G = -DSIGN(DSQRT(H),F)                                         
         E(I) = SCALE*G                                                 
         H = H-F*G                                                      
         A(NN) = F-G                                                    
         IF (L .EQ. 1) GO TO 55                                         
         F = ZERO                                                       
         JK1 = 1                                                        
         DO 40 J = 1,L                                                  
            G = ZERO                                                    
            IK = NBEG+1                                                 
            JK = JK1                                                    
C                                  FORM ELEMENT OF A*U                  
            DO 25 K = 1,J                                               
               G = G+A(JK)*A(IK)                                        
               JK = JK+1                                                
               IK = IK+1                                                
   25       CONTINUE                                                    
            JP1 = J+1                                                   
            IF (L .LT. JP1) GO TO 35                                    
            JK = JK+J-1                                                 
            DO 30 K = JP1,L                                             
               G = G+A(JK)*A(IK)                                        
               JK = JK+K                                                
               IK = IK+1                                                
   30       CONTINUE                                                    
C                                  FORM ELEMENT OF P                    
   35       E(J) = G/H                                                  
            F = F+E(J)*A(NBEG+J)                                        
            JK1 = JK1+J                                                 
   40    CONTINUE                                                       
         HH = F/(H+H)                                                   
C                                  FORM REDUCED A                       
         JK = 1                                                         
         DO 50 J = 1,L                                                  
            F = A(NBEG+J)                                               
            G = E(J)-HH*F                                               
            E(J) = G                                                    
            DO 45 K = 1,J                                               
               A(JK) = A(JK)-F*E(K)-G*A(NBEG+K)                         
               JK = JK+1                                                
   45       CONTINUE                                                    
   50    CONTINUE                                                       
   55    DO 60 K = 1,L                                                  
            A(NBEG+K) = SCALE*A(NBEG+K)                                 
   60    CONTINUE                                                       
   65    D(I) = A(NBEG+I)                                               
         A(NBEG+I) = H*SCALE*SCALE                                      
         NBEG = NBEG-I+1                                                
         NN = NN-I                                                      
   70 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - EIGRS                                         
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - JUNE 1, 1980                                  
C                                                                       
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF  
C                           A REAL SYMMETRIC MATRIX                     
C                                                                       
C   USAGE               - CALL EIGRS (A,N,JOBN,D,Z,IZ,WK,IER)           
C                                                                       
C   ARGUMENTS    A      - INPUT REAL SYMMETRIC MATRIX OF ORDER N,       
C                           WHOSE EIGENVALUES AND EIGENVECTORS          
C                           ARE TO BE COMPUTED. INPUT A IS              
C                           DESTROYED IF IJOB IS EQUAL TO 0 OR 1.       
C                N      - INPUT ORDER OF THE MATRIX A.                  
C                JOBN   - INPUT OPTION PARAMETER.  IF JOBN.GE.10        
C                         A IS ASSUMED TO BE IN FULL STORAGE MODE       
C                         (IN THIS CASE, A MUST BE DIMENSIONED EXACTLY  
C                         N BY N IN THE CALLING PROGRAM).               
C                         IF JOBN.LT.10 THEN A IS ASSUMED TO BE IN      
C                         SYMMETRIC STORAGE MODE.  DEFINE               
C                         IJOB=MOD(JOBN,10).  THEN WHEN                 
C                           IJOB = 0, COMPUTE EIGENVALUES ONLY          
C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-    
C                             VECTORS.                                  
C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS 
C                             AND PERFORMANCE INDEX.                    
C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.   
C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS 
C                           RETURNED IN WK(1). THE ROUTINES HAVE        
C                           PERFORMED (WELL, SATISFACTORILY, POORLY) IF 
C                           WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,   
C                           GREATER THAN 100).                          
C                D      - OUTPUT VECTOR OF LENGTH N,                    
C                           CONTAINING THE EIGENVALUES OF A.            
C                Z      - OUTPUT N BY N MATRIX CONTAINING               
C                           THE EIGENVECTORS OF A.                      
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-    
C                           PONDS TO THE EIGENVALUE D(J).               
C                           IF IJOB = 0, Z IS NOT USED.                 
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS    
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE 
C                           CALLING PROGRAM.                            
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS           
C                           ON THE VALUE OF IJOB, WHEN                  
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST N.   
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST N.   
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST      
C                             N(N+1)/2+N.                               
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.   
C                IER    - ERROR PARAMETER (OUTPUT)                      
C                         TERMINAL ERROR                                
C                           IER = 128+J, INDICATES THAT EQRT2S FAILED   
C                             TO CONVERGE ON EIGENVALUE J. EIGENVALUES  
C                             AND EIGENVECTORS 1,...,J-1 HAVE BEEN      
C                             COMPUTED CORRECTLY, BUT THE EIGENVALUES   
C                             ARE UNORDERED. THE PERFORMANCE INDEX      
C                             IS SET TO 1000.0                          
C                         WARNING ERROR (WITH FIX)                      
C                           IN THE FOLLOWING, IJOB = MOD(JOBN,10).      
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR  
C                             IJOB IS GREATER THAN 3. IJOB SET TO 1.    
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO    
C                             ZERO, AND IZ IS LESS THAN THE ORDER OF    
C                             MATRIX A. IJOB IS SET TO ZERO.            
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - EHOBKS,EHOUSS,EQRT2S,UERTST,UGETIO            
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE EIGRS  (A,N,JOBN,D,Z,IZ,WK,IER)                        
C                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER            N,JOBN,IZ,IER                                  
      DOUBLE PRECISION   A(1),D(1),WK(1),Z(IZ,1)                        
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER            IJOB,IR,JR,IJ,JI,NP1                           
      INTEGER            JER,NA,ND,IIZ,IBEG,IL,KK,LK,I,J,K,L            
      DOUBLE PRECISION   ANORM,ASUM,PI,SUMZ,SUMR,AN,S,TEN,RDELP,ZERO,   
     1                   ONE,THOUS                                      
      DATA               RDELP/2.775557562D-17/                         
      DATA               ZERO,ONE/0.0D0,1.0D0/,TEN/10.0D0/,THOUS/1000.0D
     +0/                                                                
C                                  INITIALIZE ERROR PARAMETERS          
C                                  FIRST EXECUTABLE STATEMENT           
      IER = 0                                                           
      JER = 0                                                           
      IF (JOBN.LT.10) GO TO 15                                          
C                                  CONVERT TO SYMMETRIC STORAGE MODE    
      K = 1                                                             
      JI = N-1                                                          
      IJ = 1                                                            
      DO 10 J=1,N                                                       
         DO 5 I=1,J                                                     
            A(K) = A(IJ)                                                
            IJ = IJ+1                                                   
            K = K+1                                                     
    5    CONTINUE                                                       
         IJ = IJ + JI                                                   
         JI = JI - 1                                                    
   10 CONTINUE                                                          
   15 IJOB = MOD(JOBN,10)                                               
      IF (IJOB.GE.0.AND.IJOB.LE.3) GO TO 20                             
C                                  WARNING ERROR - IJOB IS NOT IN THE   
C                                    RANGE                              
      IER = 66                                                          
      IJOB = 1                                                          
      GO TO 25                                                          
   20 IF (IJOB.EQ.0) GO TO 35                                           
   25 IF (IZ.GE.N) GO TO 30                                             
C                                  WARNING ERROR - IZ IS LESS THAN N    
C                                    EIGENVECTORS CAN NOT BE COMPUTED,  
C                                    IJOB SET TO ZERO                   
      IER = 67                                                          
      IJOB = 0                                                          
   30 IF (IJOB.EQ.3) GO TO 75                                           
   35 NA = (N*(N+1))/2                                                  
      IF (IJOB.NE.2) GO TO 45                                           
      DO 40 I=1,NA                                                      
         WK(I) = A(I)                                                   
   40 CONTINUE                                                          
C                                  SAVE INPUT A IF IJOB = 2             
   45 ND = 1                                                            
      IF (IJOB.EQ.2) ND = NA+1                                          
C                                  REDUCE A TO SYMMETRIC TRIDIAGONAL    
C                                    FORM                               
      CALL EHOUSS (A,N,D,WK(ND),WK(ND))                                 
      IIZ = 1                                                           
      IF (IJOB.EQ.0) GO TO 60                                           
      IIZ = IZ                                                          
C                                  SET Z TO THE IDENTITY MATRIX         
      DO 55 I=1,N                                                       
         DO 50 J=1,N                                                    
            Z(I,J) = ZERO                                               
   50    CONTINUE                                                       
         Z(I,I) = ONE                                                   
   55 CONTINUE                                                          
C                                  COMPUTE EIGENVALUES AND EIGENVECTORS 
   60 CALL EQRT2S (D,WK(ND),N,Z,IIZ,JER)                                
      IF (IJOB.EQ.0) GO TO 9000                                         
      IF (JER.GT.128) GO TO 65                                          
C                                  BACK TRANSFORM EIGENVECTORS          
      CALL EHOBKS (A,N,1,N,Z,IZ)                                        
   65 IF (IJOB.LE.1) GO TO 9000                                         
C                                  MOVE INPUT MATRIX BACK TO A          
      DO 70 I=1,NA                                                      
         A(I) = WK(I)                                                   
   70 CONTINUE                                                          
      WK(1) = THOUS                                                     
      IF (JER.NE.0) GO TO 9000                                          
C                                  COMPUTE 1 - NORM OF A                
   75 ANORM = ZERO                                                      
      IBEG = 1                                                          
      DO 85 I=1,N                                                       
         ASUM = ZERO                                                    
         IL = IBEG                                                      
         KK = 1                                                         
         DO 80 L=1,N                                                    
            ASUM = ASUM+DABS(A(IL))                                     
            IF (L.GE.I) KK = L                                          
            IL = IL+KK                                                  
   80    CONTINUE                                                       
         ANORM = DMAX1(ANORM,ASUM)                                      
         IBEG = IBEG+I                                                  
   85 CONTINUE                                                          
      IF (ANORM.EQ.ZERO) ANORM = ONE                                    
C                                  COMPUTE PERFORMANCE INDEX            
      PI = ZERO                                                         
      DO 100 I=1,N                                                      
         IBEG = 1                                                       
         S = ZERO                                                       
         SUMZ = ZERO                                                    
         DO 95 L=1,N                                                    
            LK = IBEG                                                   
            KK = 1                                                      
            SUMZ = SUMZ+DABS(Z(L,I))                                    
            SUMR = -D(I)*Z(L,I)                                         
            DO 90 K=1,N                                                 
               SUMR = SUMR+A(LK)*Z(K,I)                                 
               IF (K.GE.L) KK = K                                       
               LK = LK+KK                                               
   90       CONTINUE                                                    
            S = S+DABS(SUMR)                                            
            IBEG = IBEG+L                                               
   95    CONTINUE                                                       
         IF (SUMZ.EQ.ZERO) GO TO 100                                    
         PI = DMAX1(PI,S/SUMZ)                                          
  100 CONTINUE                                                          
      AN = N                                                            
      PI = PI/(ANORM*TEN*AN*RDELP)                                      
      WK(1) = PI                                                        
      IF (JOBN.LT.10) GO TO 9000                                        
C                                  CONVERT BACK TO FULL STORAGE MODE    
      NP1 = N+1                                                         
      IJ = (N-1)*NP1 + 2                                                
      K = (N*(NP1))/2                                                   
      DO 110 JR=1,N                                                     
         J = NP1-JR                                                     
         DO 105 IR=1,J                                                  
            IJ = IJ-1                                                   
            A(IJ) = A(K)                                                
            K = K-1                                                     
  105    CONTINUE                                                       
         IJ = IJ-JR                                                     
  110 CONTINUE                                                          
      JI = 0                                                            
      K = N-1                                                           
      DO 120 I=1,N                                                      
         IJ = I-N                                                       
         DO 115 J=1,I                                                   
            IJ = IJ+N                                                   
            JI = JI+1                                                   
            A(IJ) = A(JI)                                               
  115    CONTINUE                                                       
         JI = JI + K                                                    
         K = K-1                                                        
  120 CONTINUE                                                          
 9000 CONTINUE                                                          
      IF (IER.NE.0) CALL UERTST (IER,'EIGRS ')                          
      IF (JER.EQ.0) GO TO 9005                                          
      IER = JER                                                         
      CALL UERTST (IER,'EIGRS ')                                        
 9005 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - EQRT2S                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - NOVEMBER 1, 1984                              
C                                                                       
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF  
C                           A SYMMETRIC TRIDIAGONAL MATRIX USING THE    
C                           QL METHOD.                                  
C                                                                       
C   USAGE               - CALL EQRT2S (D,E,N,Z,IZ,IER)                  
C                                                                       
C   ARGUMENTS    D      - ON INPUT, THE VECTOR D OF LENGTH N CONTAINS   
C                           THE DIAGONAL ELEMENTS OF THE SYMMETRIC      
C                           TRIDIAGONAL MATRIX T.                       
C                           ON OUTPUT, D CONTAINS THE EIGENVALUES OF    
C                           T IN ASCENDING ORDER.                       
C                E      - ON INPUT, THE VECTOR E OF LENGTH N CONTAINS   
C                           THE SUB-DIAGONAL ELEMENTS OF T IN POSITION  
C                           2,...,N. ON OUTPUT, E IS DESTROYED.         
C                N      - ORDER OF TRIDIAGONAL MATRIX T.(INPUT)         
C                Z      - ON INPUT, Z CONTAINS THE IDENTITY MATRIX OF   
C                           ORDER N.                                    
C                           ON OUTPUT, Z CONTAINS THE EIGENVECTORS      
C                           OF T. THE EIGENVECTOR IN COLUMN J OF Z      
C                           CORRESPONDS TO THE EIGENVALUE D(J).         
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS    
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE 
C                           CALLING PROGRAM. IF IZ IS LESS THAN N, THE  
C                           EIGENVECTORS ARE NOT COMPUTED. IN THIS CASE 
C                           Z IS NOT USED.                              
C                IER    - ERROR PARAMETER                               
C                         TERMINAL ERROR                                
C                           IER = 128+J, INDICATES THAT EQRT2S FAILED   
C                             TO CONVERGE ON EIGENVALUE J. EIGENVALUES  
C                             AND EIGENVECTORS 1,...,J-1 HAVE BEEN      
C                             COMPUTED CORRECTLY, BUT THE EIGENVALUES   
C                             ARE UNORDERED.                            
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE EQRT2S (D,E,N,Z,IZ,IER)                                
C                                                                       
      DIMENSION          D(1),E(1),Z(IZ,1)                              
      DOUBLE PRECISION   D,E,Z,B,C,F,G,H,P,R,S,RDELP,ONE,ZERO           
      DATA               RDELP/2.775557562D-17/                         
      DATA               ZERO,ONE/0.0D0,1.0D0/                          
C                                  MOVE THE LAST N-1 ELEMENTS           
C                                  OF E INTO THE FIRST N-1 LOCATIONS    
C                                  FIRST EXECUTABLE STATEMENT           
      IER  = 0                                                          
      IF (N .EQ. 1) GO TO 9005                                          
      DO 5  I=2,N                                                       
         E(I-1) = E(I)                                                  
    5 CONTINUE                                                          
      E(N) = ZERO                                                       
      B = ZERO                                                          
      F = ZERO                                                          
      DO  60  L=1,N                                                     
         J = 0                                                          
         H = RDELP*(DABS(D(L))+DABS(E(L)))                              
         IF (B.LT.H) B = H                                              
C                                  LOOK FOR SMALL SUB-DIAGONAL ELEMENT  
         DO 10  M=L,N                                                   
            K=M                                                         
            IF (DABS(E(K)) .LE. B) GO TO 15                             
   10    CONTINUE                                                       
   15    M = K                                                          
         IF (M.EQ.L) GO TO 55                                           
   20    IF (J .EQ. 30) GO TO 85                                        
         J = J+1                                                        
         L1 = L+1                                                       
         G = D(L)                                                       
         P = (D(L1)-G)/(E(L)+E(L))                                      
         R = DABS(P)                                                    
         IF (RDELP*DABS(P) .LT. 1.0D0) R = DSQRT(P*P+ONE)               
         D(L) = E(L)/(P+DSIGN(R,P))                                     
         H = G-D(L)                                                     
         DO 25 I = L1,N                                                 
            D(I) = D(I)-H                                               
   25    CONTINUE                                                       
         F = F+H                                                        
C                                  QL TRANSFORMATION                    
         P = D(M)                                                       
         C = ONE                                                        
         S = ZERO                                                       
         MM1 = M-1                                                      
         MM1PL = MM1+L                                                  
         IF (L.GT.MM1) GO TO 50                                         
         DO 45 II=L,MM1                                                 
            I = MM1PL-II                                                
            G = C*E(I)                                                  
            H = C*P                                                     
            IF (DABS(P).LT.DABS(E(I))) GO TO 30                         
            C = E(I)/P                                                  
            R = DSQRT(C*C+ONE)                                          
            E(I+1) = S*P*R                                              
            S = C/R                                                     
            C = ONE/R                                                   
            GO TO 35                                                    
   30       C = P/E(I)                                                  
            R = DSQRT(C*C+ONE)                                          
            E(I+1) = S*E(I)*R                                           
            S = ONE/R                                                   
            C = C*S                                                     
   35       P = C*D(I)-S*G                                              
            D(I+1) = H+S*(C*G+S*D(I))                                   
            IF (IZ .LT. N) GO TO 45                                     
C                                  FORM VECTOR                          
            DO 40 K=1,N                                                 
               H = Z(K,I+1)                                             
               Z(K,I+1) = S*Z(K,I)+C*H                                  
               Z(K,I) = C*Z(K,I)-S*H                                    
   40       CONTINUE                                                    
   45    CONTINUE                                                       
   50    E(L) = S*P                                                     
         D(L) = C*P                                                     
         IF (DABS(E(L)) .GT.B) GO TO 20                                 
   55    D(L) = D(L) + F                                                
   60 CONTINUE                                                          
C                                  ORDER EIGENVALUES AND EIGENVECTORS   
      DO  80  I=1,N                                                     
         K = I                                                          
         P = D(I)                                                       
         IP1 = I+1                                                      
         IF (IP1.GT.N) GO TO 70                                         
         DO 65  J=IP1,N                                                 
            IF (D(J) .GE. P) GO TO 65                                   
            K = J                                                       
            P = D(J)                                                    
   65    CONTINUE                                                       
   70    IF (K.EQ.I) GO TO 80                                           
         D(K) = D(I)                                                    
         D(I) = P                                                       
         IF (IZ .LT. N) GO TO 80                                        
         DO 75 J = 1,N                                                  
            P = Z(J,I)                                                  
            Z(J,I) = Z(J,K)                                             
            Z(J,K) = P                                                  
   75    CONTINUE                                                       
   80 CONTINUE                                                          
      GO TO 9005                                                        
   85 IER = 128+L                                                       
      CONTINUE
      CALL UERTST(IER,'EQRT2S')                                         
 9005 RETURN                                                            
      END                                                               
C   IMSL ROUTINE NAME   - LINV1P                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C   COMPUTER            - VAX/DOUBLE                                    
C                                                                       
C   LATEST REVISION     - JANUARY 1, 1978                               
C                                                                       
C   PURPOSE             - INVERSION OF MATRIX - POSITIVE DEFINITE -     
C                           SYMMETRIC STORAGE MODE - SPACE ECONOMIZER   
C                           SOLUTION                                    
C                                                                       
C   USAGE               - CALL LINV1P (A,N,AINV,IDGT,D1,D2,IER)         
C                                                                       
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING    
C                           THE N BY N POSITIVE DEFINITE SYMMETRIC      
C                           MATRIX TO BE INVERTED. A IS STORED IN       
C                           SYMMETRIC STORAGE MODE.                     
C                         ON OUTPUT, A IS REPLACED BY THE LOWER         
C                           TRIANGULAR MATRIX L WHERE A = L*L-TRANSPOSE.
C                           L IS STORED IN SYMMETRIC STORAGE MODE WITH  
C                           THE DIAGONAL ELEMENTS OF L STORED IN        
C                           RECIPROCAL FORM.                            
C                N      - ORDER OF A. (INPUT)                           
C                AINV   - OUTPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING   
C                           THE N BY N INVERSE OF A. AINV IS STORED IN  
C                           SYMMETRIC STORAGE MODE.                     
C                IDGT   - THE ELEMENTS OF A ARE ASSUMED TO BE CORRECT   
C                           TO IDGT DECIMAL DIGITS. (CURRENTLY NOT USED)
C                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.           
C                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)        
C                IER    - ERROR PARAMETER. (OUTPUT)                     
C                         TERMINAL ERROR                                
C                           IER = 129 INDICATES THAT THE ORIGINAL       
C                             MATRIX A IS ALGORITHMICALLY NOT POSITIVE  
C                             DEFINITE.  (SEE THE CHAPTER L PRELUDE).   
C                                                                       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
C                       - SINGLE/H36,H48,H60                            
C                                                                       
C   REQD. IMSL ROUTINES - LUDECP,LUELMP,UERTST,UGETIO                   
C                                                                       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C                                                                       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C                                                                       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE LINV1P (A,N,AINV,D1,D2,IER)
C                                                                       
      DIMENSION          A(1),AINV(1)                                   
      DOUBLE PRECISION   A,AINV,ZERO,ONE,D1,D2                          
      DATA               ZERO,ONE/0.0D0,1.0D0/                          
C                                  FIRST EXECUTABLE STATEMENT           
      IER=0                                                             
      L = 1                                                             
      N1 = N-1                                                          
      LN = N                                                            
C                                  DECOMPOSE A                          
      CALL LUDECP(A,A,N,D1,D2,IER)                                      
      IF (IER .NE. 0) GO TO 9000                                        
      DO 10 I = 1,N                                                     
         DO 5 J = L,LN                                                  
            AINV(J) = ZERO                                              
    5    CONTINUE                                                       
         AINV(L+I-1) = ONE                                              
C                                  OBTAIN THE SOLUTION                  
         CALL LUELMP(A,AINV(L),N,AINV(L))                               
         L = L+I                                                        
         LN = L+N1                                                      
   10 CONTINUE                                                          
      GO TO 9005                                                        
 9000 CONTINUE                                                          
      CALL UERTST(IER,'LINV1P')                                         
 9005 RETURN                                                            
      END                                                               
