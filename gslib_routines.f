********************************************************************************
***                                                                          ***
***                         AUXILIAR GSLIB ROUTINES                          ***
***                                                                          ***
********************************************************************************

      REAL*8 function acorni(ixv)
c-----------------------------------------------------------------------
c
c Fortran implementation of ACORN random number generator of order less
c than or equal to 12 (higher orders can be obtained by increasing the
c parameter value MAXORD).
c
c
c NOTES: 1. The variable idum is a dummy variable. The common block
c           IACO is used to transfer data into the function.
c
c        2. Before the first call to ACORN the common block IACO must
c           be initialised by the user, as follows. The values of
c           variables in the common block must not subsequently be
c           changed by the user.
c
c             KORDEI - order of generator required ( must be =< MAXORD)
c
c             MAXINT - modulus for generator, must be chosen small
c                      enough that 2*MAXINT does not overflow
c
c             ixv(1) - seed for random number generator
c                      require 0 < ixv(1) < MAXINT
c
c             (ixv(I+1),I=1,KORDEI)
c                    - KORDEI initial values for generator
c                      require 0 =< ixv(I+1) < MAXINT
c
c        3. After initialisation, each call to ACORN generates a single
c           random number between 0 and 1.
c
c        4. An example of suitable values for parameters is
c
c             KORDEI   = 10
c             MAXINT   = 2**30
c             ixv(1)   = an odd integer in the (approximate) range 
c                        (0.001 * MAXINT) to (0.999 * MAXINT)
c             ixv(I+1) = 0, I=1,KORDEI
c
c
c
c Author: R.S.Wikramaratna,                           Date: October 1990
c-----------------------------------------------------------------------

C______________________________________________________ Declaration of variables

      IMPLICIT NONE
      INTEGER*4 KORDEI,MAXINT,MAXOP1,IXV(13)
      INTEGER*4 I

      KORDEI=12
      MAXOP1=13
      MAXINT=2**30

      DO I=1,kordei
            IXV(I+1)=(IXV(I+1)+IXV(I))
            IF(IXV(I+1).GE.maxint) IXV(I+1)=IXV(I+1)-maxint
      END DO
      ACORNI=DBLE(IXV(kordei+1))/maxint
      RETURN
      END

********************************************************************************
********************************************************************************
********************************************************************************

      SUBROUTINE CALDRIF
     ;(IDIMCROSS      ,IESTIM         ,IO_CROSS         ,IROWCOVMEASMEAS
     ;,IROWCROSSCOV   ,MAXEQ          ,NACCEPT_SAMPLES  ,NEQUA            
     ;,NESTIM         ,NMAXKRIG       ,SCALE_UNIV       ,COVESTMEAS 
     ;,COVMEASMEAS    ,CROSS_COV      ,IDRIF            ,MEAN_DRIFT 
     ;,XKRIG          ,YKRIG          ,ZKRIG)

********************************************************************************
*
* PURPOSE   Evaluates drift functions and assign calc. values at right
*           positions of kriging matrix and independent term. If so desired,
*           saves cross-covariances (restrictions due to UK in this case)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVESTMEAS             RHS of kriging system
*  COVMEASMEAS            Kriging matrix                              
*  CROSS_COV              Cross covariance matrix Qyz                
*  IDRIF                  Array with options for evaluating drift functions
*  MEAN_DRIFT             Values of drift functions at esimation point
*  XKRIG                  X-coord. of samples used for kriging
*  YKRIG                  Y-coord. of samples used for kriging
*  ZKRIG                  Z-coord. of samples used for kriging
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMCROSS              Used to dimension CROSS_COV
*  IESTIM                 Point/block being estimated
*  IO_CROSS               If <>0, saves restrictions at CROSS_COV
*  IROWCOVMEASMEAS        Pointer to kriging matrix and kriging RHS vector
*  IROWCROSSCOV           Pointer to CROSS_COV
*  MAXEQ                  Maximum number of equations for kriging a block (main
*                         dimension of kriging matrix and RHS)
*  NACCEPT_SAMPLES        Number of samples used for kriging actual point
*  NEQUA                  Number of kriging system equations
*  NESTIM                 Number of estimation points                    
*  NMAXKRIG               Max. number of samples for krig. a estim. point/block
*  SCALE_UNIV             Rescaling factor to make solution more stable
*
* INTERNAL VARIABLES: SCALARS
*
*  ISAMPLE                Dummy counter of samples         
*
* HISTORY    AAR (First coding, Jan 2003)
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE

      INTEGER*4 MAXEQ,NMAXKRIG,IROWCOVMEASMEAS,IROWCROSSCOV,IO_CROSS
     ;         ,NACCEPT_SAMPLES,IESTIM,IDIMCROSS,NESTIM,NEQUA
     ;         ,IDRIF(9)

      REAL*8 SCALE_UNIV
     ;      ,COVMEASMEAS(MAXEQ*MAXEQ),COVESTMEAS(MAXEQ),XKRIG(NMAXKRIG)
     ;      ,YKRIG(NMAXKRIG),ZKRIG(NMAXKRIG),MEAN_DRIFT(9)
     ;      ,CROSS_COV(IDIMCROSS,NESTIM)

      INTEGER*4 ISAMPLE

C_______________________ Step 1: Evaluates drift functions and assign calculated 
C_______________________         values at right positions of kriging matrix and 
C_______________________         independent term

      IF (IDRIF(1).EQ.1) THEN     ! Linear term (X)
         IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
         IROWCROSSCOV = IROWCROSSCOV + 1
         DO ISAMPLE=1,NACCEPT_SAMPLES
            COVMEASMEAS (NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) = 
     ;         XKRIG(ISAMPLE)*SCALE_UNIV
            COVMEASMEAS (NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) = 
     ;         XKRIG(ISAMPLE)*SCALE_UNIV
         END DO
         COVESTMEAS(IROWCOVMEASMEAS) = MEAN_DRIFT(1)
         IF (IO_CROSS.NE.0) CROSS_COV(IROWCROSSCOV,IESTIM)=MEAN_DRIFT(1)
      END IF
  
      IF (IDRIF(2).EQ.1) THEN     ! Linear term (Y)
         IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
         IROWCROSSCOV = IROWCROSSCOV + 1
         DO ISAMPLE=1,NACCEPT_SAMPLES
            COVMEASMEAS(NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) = 
     ;         YKRIG(ISAMPLE)*SCALE_UNIV
            COVMEASMEAS(NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) = 
     ;         YKRIG(ISAMPLE)*SCALE_UNIV
         END DO
         COVESTMEAS(IROWCOVMEASMEAS) = MEAN_DRIFT(2)
         IF (IO_CROSS.NE.0) CROSS_COV(IROWCROSSCOV,IESTIM)=MEAN_DRIFT(2)
      END IF
 
      IF (IDRIF(3).EQ.1) THEN     ! Linear term (Z)
         IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
         IROWCROSSCOV = IROWCROSSCOV + 1
         DO ISAMPLE=1,NACCEPT_SAMPLES
            COVMEASMEAS(NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) = 
     ;         ZKRIG(ISAMPLE)*SCALE_UNIV
            COVMEASMEAS(NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) = 
     ;         ZKRIG(ISAMPLE)*SCALE_UNIV
         END DO
         COVESTMEAS(IROWCOVMEASMEAS) = MEAN_DRIFT(3)
         IF (IO_CROSS.NE.0) CROSS_COV(IROWCROSSCOV,IESTIM)=MEAN_DRIFT(3)
      END IF
 
      IF (IDRIF(4).EQ.1) THEN     ! Cuadratic term (X*X)
         IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
         IROWCROSSCOV = IROWCROSSCOV + 1
         DO ISAMPLE=1,NACCEPT_SAMPLES
            COVMEASMEAS(NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) = 
     ;         XKRIG(ISAMPLE)*XKRIG(ISAMPLE)*SCALE_UNIV
            COVMEASMEAS(NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) = 
     ;         XKRIG(ISAMPLE)*XKRIG(ISAMPLE)*SCALE_UNIV
         END DO
         COVESTMEAS(IROWCOVMEASMEAS) = MEAN_DRIFT(4)
         IF (IO_CROSS.NE.0) CROSS_COV(IROWCROSSCOV,IESTIM)=MEAN_DRIFT(4)
      END IF
  
      IF (IDRIF(5).EQ.1) THEN     ! Cuadratic term (Y*Y)
         IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
         IROWCROSSCOV = IROWCROSSCOV + 1
         DO ISAMPLE=1,NACCEPT_SAMPLES
            COVMEASMEAS(NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) = 
     ;         YKRIG(ISAMPLE)*YKRIG(ISAMPLE)*SCALE_UNIV
            COVMEASMEAS(NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) = 
     ;         YKRIG(ISAMPLE)*YKRIG(ISAMPLE)*SCALE_UNIV
         END DO
         COVESTMEAS(IROWCOVMEASMEAS) = MEAN_DRIFT(5)
         IF (IO_CROSS.NE.0) CROSS_COV(IROWCROSSCOV,IESTIM)=MEAN_DRIFT(5)
      END IF
  
      IF (IDRIF(6).EQ.1) THEN     ! Cuadratic term (Z*Z)
         IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
         IROWCROSSCOV = IROWCROSSCOV + 1
         DO ISAMPLE=1,NACCEPT_SAMPLES
            COVMEASMEAS(NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) = 
     ;         ZKRIG(ISAMPLE)*ZKRIG(ISAMPLE)*SCALE_UNIV
            COVMEASMEAS(NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) = 
     ;         ZKRIG(ISAMPLE)*ZKRIG(ISAMPLE)*SCALE_UNIV
         END DO
         COVESTMEAS(IROWCOVMEASMEAS) = MEAN_DRIFT(6)
         IF(IO_CROSS.NE.0) CROSS_COV(IROWCROSSCOV,IESTIM)=MEAN_DRIFT(6)
      END IF

      IF (IDRIF(7).EQ.1) THEN     ! Cross-Cuadratic term (X*Y)
         IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
         IROWCROSSCOV = IROWCROSSCOV + 1
         DO ISAMPLE=1,NACCEPT_SAMPLES
            COVMEASMEAS(NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) = 
     ;         XKRIG(ISAMPLE)*YKRIG(ISAMPLE)*SCALE_UNIV
            COVMEASMEAS(NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) = 
     ;         XKRIG(ISAMPLE)*YKRIG(ISAMPLE)*SCALE_UNIV
         END DO
         COVESTMEAS(IROWCOVMEASMEAS) = MEAN_DRIFT(7)
         IF (IO_CROSS.NE.0) CROSS_COV(IROWCROSSCOV,IESTIM)=MEAN_DRIFT(7)
      END IF
  
      IF (IDRIF(8).EQ.1) THEN     ! Cross-Cuadratic term (X*Z)
         IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
         IROWCROSSCOV = IROWCROSSCOV + 1
         DO ISAMPLE=1,NACCEPT_SAMPLES
            COVMEASMEAS(NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) = 
     ;         XKRIG(ISAMPLE)*ZKRIG(ISAMPLE)*SCALE_UNIV
            COVMEASMEAS(NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) = 
     ;         XKRIG(ISAMPLE)*ZKRIG(ISAMPLE)*SCALE_UNIV
         END DO
         COVESTMEAS(IROWCOVMEASMEAS) = MEAN_DRIFT(8)
         IF (IO_CROSS.NE.0) CROSS_COV(IROWCROSSCOV,IESTIM)=MEAN_DRIFT(8)
      END IF
  
      IF (IDRIF(9).EQ.1) THEN     ! Cross-Cuadratic term (Y*Z)
         IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
         IROWCROSSCOV = IROWCROSSCOV + 1
         DO ISAMPLE=1,NACCEPT_SAMPLES
            COVMEASMEAS(NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) = 
     ;         YKRIG(ISAMPLE)*ZKRIG(ISAMPLE)*SCALE_UNIV
            COVMEASMEAS(NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) = 
     ;         YKRIG(ISAMPLE)*ZKRIG(ISAMPLE)*SCALE_UNIV
         END DO
         COVESTMEAS(IROWCOVMEASMEAS) = MEAN_DRIFT(9)
         IF (IO_CROSS.NE.0) CROSS_COV(IROWCROSSCOV,IESTIM)=MEAN_DRIFT(9)
      END IF

      RETURN
      END

********************************************************************************
********************************************************************************
********************************************************************************

      subroutine COVARIANCE
     ;(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa,irot,MAXROT
     ;,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c                      5. hole effect model
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested 
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c-----------------------------------------------------------------------

C____________________ Declaration of variables


      IMPLICIT NONE

      INTEGER*4 IVARG,MAXNST,IROT,MAXROT

      INTEGER*4 ISTART,IS,IST,IR

      REAL*8 X1,Y1,Z1,X2,Y2,Z2,CMAX,COVA,H,HR

      REAL*8 PI,EPSLON,PMX

      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-10)
      integer*4   nst(*),it(*)
      real*8      c0(*),cc(*),aa(*)
      real*8    rotmat(MAXROT,3,3),hsqd,sqdist
c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c
      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
      do is=1,nst(ivarg)
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
            h = real(dsqrt(hsqd))
c
c Spherical Variogram Model?
c
            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c
            else if(it(ist).eq.2) then

c           a is range and integral distance. 3a is effective range

                  cova = cova + cc(ist)*exp(-1.0D0*h/aa(ist))
*                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))

c
c Gaussian Variogram Model?
c
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-(3.0*h/aa(ist))
     +                                      *(3.0*h/aa(ist)))
c
c Power Variogram Model?
c
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*(h**aa(ist))
c
c Hole Effect Model?
c
            else if(it(ist).eq.5) then
c                 d = 10.0 * aa(ist)
c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
            endif
      end do
c
c Finished:
c
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************

      subroutine getindx(n,min,siz,loc,index,inflag)
c-----------------------------------------------------------------------
c
c     Gets the coordinate index location of a point within a grid
c     ***********************************************************
c
c
c n       number of "nodes" or "cells" in this coordinate direction
c min     origin at the center of the first cell
c siz     size of the cells
c loc     location of the point being considered
c index   output index within [1,n]
c inflag  true if the location is actually in the grid (false otherwise
c         e.g., if the location is outside then index will be set to
c         nearest boundary
c
c
c
c-----------------------------------------------------------------------

C_____________________________ Declaration of variables


      IMPLICIT NONE

      INTEGER*4 N,INDEX

      REAL*8    MIN,SIZ,LOC
      logical   inflag
c
c Compute the index of "loc":
c
      index = int( (loc-min)/siz + 1.5 )
c
c Check to see if in or out:
c
      if(index.lt.1) then
            index  = 1
            inflag = .false.
      else if(index.gt.n) then
            index  = n
            inflag = .false.
      else
            inflag = .true.
      end if
c
c Return to calling program:
c
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************

      subroutine ktsol(n,ns,nv,a,b,x,ktilt,maxeq)
c-----------------------------------------------------------------------
c
c Solution of a system of linear equations by gaussian elimination with
c partial pivoting.  Several right hand side matrices and several
c variables are allowed.
c
c
c         NOTE: All input matrices must be in double precision
c
c
c INPUT/OUTPUT VARIABLES:
c
c   n                Number of equations
c   ns               Number of right hand side matrices
c   nv               Number of variables.
c   a(n*n*nv)        left hand side matrices versus columnwise.
c   b(n*ns*nv)       input right hand side matrices.
c   x(n*ns*nv)       solution matrices.
c   ktilt            indicator of singularity
c                      =  0  everything is ok.
c                      = -1 n.le.1
c                      =  k  a null pivot appeared at the kth iteration.
c   tol              used in test for null pivot. depends on machine
c                      precision and can also be set for the tolerance
c                      of an ill-defined kriging system.
c
c
c-----------------------------------------------------------------------

C___________________ Declaration of variables


      IMPLICIT NONE

      INTEGER*4 N,NS,NV,KTILT,MAXEQ

      INTEGER*4 NTN,NM1,IV,NVA,NVB,K,KP1,KDIAG,NPIV,IPIV,I1,I,J1,J2,J,I2
     ;         ,NVB1,NVB2,IL,NMK,KB,ITOT


      REAL*8 X(MAXEQ),A(MAXEQ*MAXEQ),B(MAXEQ)
      REAL*8 TOL,T

c
c Make sure there are equations to solve:
c
      if(n.le.1) then
            ktilt = -1
            return
      endif
c
c Initialization:
c
      tol   = 0.1e-10
      ktilt = 0
      ntn   = n*n
      nm1   = n-1
c
c Triangulation is done variable by variable:
c
      do iv=1,nv
c
c Indices of location in vectors a and b:
c
            nva = ntn*(iv-1)
            nvb = n*ns*(iv-1)
c
c Gaussian elimination with partial pivoting:
c
            do k=1,nm1
                  kp1 = k+1
c
c Indice of the diagonal element in the kth row:
c
                  kdiag = nva+(k-1)*n+k
c
c Find the pivot - interchange diagonal element/pivot:
c
                  npiv = kdiag
                  ipiv = k
                  i1   = kdiag
                  do i=kp1,n
                        i1 = i1+1
                        if(abs(a(i1)).gt.abs(a(npiv))) then
                              npiv = i1
                              ipiv = i
                        endif
                  end do
                  t        = a(npiv)
                  a(npiv)  = a(kdiag)
                  a(kdiag) = t
c
c Test for singularity:
c
                  if(abs(a(kdiag)).lt.tol) then
                        ktilt=k
                        return
                  endif
c
c Compute multipliers:
c
                  i1 = kdiag
                  do i=kp1,n
                        i1    = i1+1
                        a(i1) = -a(i1)/a(kdiag)
                  end do
c
c Interchange and eliminate column per column:
c
                  j1 = kdiag
                  j2 = npiv
                  do j=kp1,n
                        j1    = j1+n
                        j2    = j2+n
                        t     = a(j2)
                        a(j2) = a(j1)
                        a(j1) = t
                        i1    = j1
                        i2    = kdiag
                        do i=kp1,n
                              i1    = i1+1
                              i2    = i2+1
                              a(i1) = a(i1)+a(i2)*a(j1)
                        end do
                  end do
c
c Interchange and modify the ns right hand matrices:
c
                  i1 = nvb+ipiv
                  i2 = nvb+k
                  do i=1,ns
                        t     = b(i1)
                        b(i1) = b(i2)
                        b(i2) = t
                        j1    = i2
                        j2    = kdiag
                        do j=kp1,n
                              j1    = j1+1
                              j2    = j2+1
                              b(j1) = b(j1)+b(i2)*a(j2)
                        end do
                        i1 = i1+n
                        i2 = i2+n
                  end do
            end do
c
c Test for singularity for the last pivot:
c
            kdiag = ntn*iv
            if(abs(a(kdiag)).lt.tol) then
                  ktilt = n
                  return
            endif
      end do
c
c End of triangulation. Now, solve back variable per variable:
c
      do iv=1,nv
c
c Indices of location in vectors a and b:
c
            nva  = ntn*iv
            nvb1 = n*ns*(iv-1)+1
            nvb2 = n*ns*iv
c
c Back substitution with the ns right hand matrices:
c
            do il=1,ns
                  do k=1,nm1
                        nmk = n-k
c
c Indice of the diagonal element of the (n-k+1)th row and of
c the (n-k+1)th element of the left hand side.
c
                        kdiag = nva-(n+1)*(k-1)
                        kb    = nvb2-(il-1)*n-k+1
                        b(kb) = b(kb)/a(kdiag)
                        t     = -b(kb)
                        i1    = kb
                        i2    = kdiag
                        do i=1,nmk
                              i1    = i1-1
                              i2    = i2-1
                              b(i1) = b(i1)+a(i2)*t
                        end do
                  end do
                  kdiag = kdiag-n-1
                  kb    = kb-1
                  b(kb) = b(kb)/a(kdiag)
            end do
c
c End of back substitution:
c
      end do
c
c Restitution of the solution:
c
      itot = n*ns*nv
      do i=1,itot
            x(i) = b(i)
      end do
c
c Finished:
c
      return
      end


********************************************************************************
********************************************************************************
********************************************************************************

      subroutine WHICH_SUPER_BLOCK
     ;(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,irot,MAXROT,rotmat
     ;,radsqd,nsbtosr,ixsbtosr,iysbtosr,izsbtosr)

c-----------------------------------------------------------------------
c
c             Establish Which Super Blocks to Search
c             **************************************
c
c This subroutine establishes which super blocks must be searched given
c that a point being estimated/simulated falls within a super block
c centered at 0,0,0.
c
c
c
c INPUT VARIABLES:
c
c   nxsup,xsizsup    Definition of the X super block grid
c   nysup,ysizsup    Definition of the Y super block grid
c   nzsup,zsizsup    Definition of the Z super block grid
c   irot             index of the rotation matrix for searching
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c   radsqd           squared search radius
c
c
c
c OUTPUT VARIABLES:
c
c   nsbtosr          Number of super blocks to search
c   ixsbtosr         X offsets for super blocks to search
c   iysbtosr         Y offsets for super blocks to search
c   izsbtosr         Z offsets for super blocks to search
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c
c
c
c-----------------------------------------------------------------------

C___________________ Declaration of variables

      IMPLICIT NONE

      INTEGER*4 NXSUP,NYSUP,NZSUP,IROT,MAXROT,NSBTOSR

      INTEGER*4 I,J,K,I1,J1,K1,I2,J2,K2

      REAL*8 XSIZSUP,YSIZSUP,ZSIZSUP,RADSQD

      REAL*8 HSQD,SQDIST,SHORTEST,XO,YO,ZO,XDIS,YDIS,ZDIS

C___________________ Declaration of external arrays

      real*8    ROTMAT(MAXROT,3,3)
      integer*4 IXSBTOSR(*),IYSBTOSR(*),IZSBTOSR(*)
c
c MAIN Loop over all possible super blocks:
c
      nsbtosr = 0
      do i=-(nxsup-1),(nxsup-1)
      do j=-(nysup-1),(nysup-1)
      do k=-(nzsup-1),(nzsup-1)
            xo = real(i)*xsizsup
            yo = real(j)*ysizsup
            zo = real(k)*zsizsup
c
c Find the closest distance between the corners of the super blocks:
c
            shortest = 1.0e21
            do i1=-1,1
            do j1=-1,1
            do k1=-1,1
                  do i2=-1,1
                  do j2=-1,1
                  do k2=-1,1
                        if(i1.ne.0.and.j1.ne.0.and.k1.ne.0.and.
     +                     i2.ne.0.and.j2.ne.0.and.k2.ne.0) then
                              xdis = real(i1-i2)*0.5*xsizsup + xo
                              ydis = real(j1-j2)*0.5*ysizsup + yo
                              zdis = real(k1-k2)*0.5*zsizsup + zo
                              hsqd = 
     ;                      sqdist(0.0d0,0.0d0,0.0d0,xdis,ydis,zdis,
     +                                      irot,MAXROT,rotmat)
                              if(hsqd.lt.shortest) shortest = hsqd
                        end if
                  end do
                  end do
                  end do
            end do
            end do
            end do
c
c Keep this super block if it is close enoutgh:
c
            if(real(shortest).le.radsqd) then
                  nsbtosr = nsbtosr + 1
                  ixsbtosr(nsbtosr) = i
                  iysbtosr(nsbtosr) = j
                  izsbtosr(nsbtosr) = k
            end if
      end do
      end do
      end do
c
c Finished:
c
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************

      subroutine setrot(ang1,ang2,ang3,anis1,anis2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c              Sets up an Anisotropic Rotation Matrix
c              **************************************
c
c Sets up the matrix to transform cartesian coordinates to coordinates
c accounting for angles and anisotropy (see manual for a detailed
c definition):
c
c
c INPUT PARAMETERS:
c
c   ang1             Azimuth angle for principal direction
c   ang2             Dip angle for principal direction
c   ang3             Third rotation angle
c   anis1            First anisotropy ratio
c   anis2            Second anisotropy ratio
c   ind              matrix indicator to initialize
c   MAXROT           maximum number of rotation matrices dimensioned
c   rotmat           rotation matrices
c
c
c NO EXTERNAL REFERENCES
c
c
c-----------------------------------------------------------------------

C_________________________ Declaration of variables



      IMPLICIT NONE

      INTEGER*4 IND,MAXROT

      REAL*8 ANG1,ANG2,ANG3,ANIS1,ANIS2

      REAL*8 DEG2RAD,EPSLON,AFAC1,AFAC2,SINA,SINB,SINT,COSA,COSB,COST
     ;      ,ALPHA,BETA,THETA

      parameter(DEG2RAD=3.141592654/180.0,EPSLON=1.e-20)

C_________________________ Array declaration

      REAL*8    ROTMAT(MAXROT,3,3)

c
c Converts the input angles to three angles which make more
c  mathematical sense:
c
c         alpha   angle between the major axis of anisotropy and the
c                 E-W axis. Note: Counter clockwise is positive.
c         beta    angle between major axis and the horizontal plane.
c                 (The dip of the ellipsoid measured positive down)
c         theta   Angle of rotation of minor axis about the major axis
c                 of the ellipsoid.
c
      if(ang1.ge.0.0.and.ang1.lt.270.0) then
            alpha = (90.0   - ang1) * DEG2RAD
      else
            alpha = (450.0  - ang1) * DEG2RAD
      endif
      beta  = -1.0 * ang2 * DEG2RAD
      theta =        ang3 * DEG2RAD
c
c Get the required sines and cosines:
c
      sina  = dble(sin(alpha))
      sinb  = dble(sin(beta))
      sint  = dble(sin(theta))
      cosa  = dble(cos(alpha))
      cosb  = dble(cos(beta))
      cost  = dble(cos(theta))
c
c Construct the rotation matrix in the required memory:
c
      afac1 = 1.0 / dble(max(anis1,EPSLON))
      afac2 = 1.0 / dble(max(anis2,EPSLON))
      rotmat(ind,1,1) =       (cosb * cosa)
      rotmat(ind,1,2) =       (cosb * sina)
      rotmat(ind,1,3) =       (-sinb)
      rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(ind,2,3) = afac1*( sint * cosb)
      rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(ind,3,3) = afac2*(cost * cosb)
c
c Return to calling program:
c
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************

      subroutine SET_SUPERBLOCK_GRID
     ;(nd,x,y,z,vr,tmp,nsec,sec1,sec2,sec3,nisb,nxsup,xmnsup,xsizsup
     ;,nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup)

c-----------------------------------------------------------------------
c
c           Establish Super Block Search Limits and Sort Data
c           *************************************************
c
c This subroutine sets up a 3-D "super block" model and orders the data
c by super block number.  The limits of the super block is set to the
c minimum and maximum limits of the grid; data outside are assigned to
c the nearest edge block.
c
c The idea is to establish a 3-D block network that contains all the
c relevant data.  The data are then sorted by their index location in
c the search network, i.e., the index location is given after knowing
c the block index in each coordinate direction (ix,iy,iz):
c          ii = (iz-1)*nxsup*nysup + (iy-1)*nxsup + ix
c An array, the same size as the number of super blocks, is constructed
c that contains the cumulative number of data in the model.  With this
c array it is easy to quickly check what data are located near any given
c location.
c
c
c
c INPUT VARIABLES:
c
c   nx,xmn,xsiz      Definition of the X grid being considered
c   ny,ymn,ysiz      Definition of the Y grid being considered
c   nz,zmn,zsiz      Definition of the Z grid being considered
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   vr(nd)           Variable at each location.
c   tmp(nd)          Temporary storage to keep track of the super block
c                      index associated to each data (uses the same
c                      storage already allocated for the simulation)
c   nsec             Number of secondary variables to carry with vr
c   sec1(nd)         First secondary variable (if nsec >= 1)
c   sec2(nd)         Second secondary variable (if nsec >= 2)
c   sec3(nd)         Third secondary variable (if nsec = 3)
c   MAXSB[X,Y,Z]     Maximum size of super block network
c
c
c
c OUTPUT VARIABLES:
c
c   nisb()                Array with cumulative number of data in each
c                           super block.
c
c EXTERNAL REFERENCES:
c
c   sortem           Sorting routine to sort the data
c
c
c
c-----------------------------------------------------------------------

C__________________________ Declaration of variables


      IMPLICIT NONE

      INTEGER*4 ND,NSEC,NXSUP,NYSUP,NZSUP

      REAL*8 XMNSUP,XSIZSUP,YMNSUP,YSIZSUP,ZMNSUP,ZSIZSUP

      INTEGER*4 I,IX,IY,IZ,II,NSORT

C__________________________ Array declaration

      REAL*8  X(*),Y(*),Z(*),VR(*),TMP(*),SEC1(*),SEC2(*),SEC3(*)
      INTEGER*4 NISB(*)
      LOGICAL INFLAG

c
c Initialize the extra super block array to zeros:
c
      do i=1,nxsup*nysup*nzsup
            nisb(i) = 0
      end do
c
c Loop over all the data assigning the data to a super block and
c accumulating how many data are in each super block:
c
      do i=1,nd
            call getindx(nxsup,xmnsup,xsizsup,x(i),ix,inflag)
            call getindx(nysup,ymnsup,ysizsup,y(i),iy,inflag)
            call getindx(nzsup,zmnsup,zsizsup,z(i),iz,inflag)
            ii = ix + (iy-1)*nxsup + (iz-1)*nxsup*nysup
            tmp(i)   = ii
            nisb(ii) = nisb(ii) + 1
      end do
c
c Sort the data by ascending super block number:
c
      nsort = 4 + nsec
      call sortem(1,nd,tmp,nsort,x,y,z,vr,sec1,sec2,sec3)
c
c Set up array nisb with the starting address of the block data:
c
      do i=1,(nxsup*nysup*nzsup-1)
            nisb(i+1) = nisb(i) + nisb(i+1)
      end do
c
c Finished:
c
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************

      subroutine sortem(ib,ie,a,iperm,b,c,d,e,f,g,h)
c-----------------------------------------------------------------------
c
c                      Quickersort Subroutine
c                      **********************
c
c This is a subroutine for sorting a real array in ascending order. This
c is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
c in collected algorithms of the ACM.
c
c The method used is that of continually splitting the array into parts
c such that all elements of one part are less than all elements of the
c other, with a third part in the middle consisting of one element.  An
c element with value t is chosen arbitrarily (here we choose the middle
c element). i and j give the lower and upper limits of the segment being
c split.  After the split a value q will have been found such that 
c a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
c performs operations on the two segments (i,q-1) and (q+1,j) as follows
c The smaller segment is split and the position of the larger segment is
c stored in the lt and ut arrays.  If the segment to be split contains
c two or fewer elements, it is sorted and another segment is obtained
c from the lt and ut arrays.  When no more segments remain, the array
c is completely sorted.
c
c
c INPUT PARAMETERS:
c
c   ib,ie        start and end index of the array to be sorteda
c   a            array, a portion of which has to be sorted.
c   iperm        0 no other array is permuted.
c                1 array b is permuted according to array a
c                2 arrays b,c are permuted.
c                3 arrays b,c,d are permuted.
c                4 arrays b,c,d,e are permuted.
c                5 arrays b,c,d,e,f are permuted.
c                6 arrays b,c,d,e,f,g are permuted.
c                7 arrays b,c,d,e,f,g,h are permuted.
c               >7 no other array is permuted.
c
c   b,c,d,e,f,g,h  arrays to be permuted according to array a.
c
c OUTPUT PARAMETERS:
c
c    a      = the array, a portion of which has been sorted.
c
c    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
c
c NO EXTERNAL ROUTINES REQUIRED:
c
c-----------------------------------------------------------------------

C_________________ Declaration of variables

      IMPLICIT NONE

      INTEGER*4 ib,ie,iperm

      INTEGER*4 IRING,i,j,k,m,p,q

      REAL*8 TA,TB,TC,TD,TE,TF,TG,TH,XA,XB,XC,XD,XE,XF,XG,XH

C_________________ Array declaration

      REAL*8 A(*),B(*),C(*),D(*),E(*),F(*),G(*),H(*)
c
c The dimensions for lt and ut have to be at least log (base 2) n
c
      integer*4 lt(64),ut(64)
c
c Initialize:
c
      j     = ie
      m     = 1
      i     = ib
      iring = iperm+1
      if (iperm.gt.7) iring=1
c
c If this segment has more than two elements  we split it
c
 10   if (j-i-1) 100,90,15
c
c p is the position of an arbitrary element in the segment we choose the
c middle element. Under certain circumstances it may be advantageous
c to choose p at random.
c
 15   p    = (j+i)/2
      ta   = a(p)
      a(p) = a(i)
      go to (21,19,18,17,16,161,162,163),iring
 163     th   = h(p)
         h(p) = h(i)
 162     tg   = g(p)
         g(p) = g(i)
 161     tf   = f(p)
         f(p) = f(i)
 16      te   = e(p)
         e(p) = e(i)
 17      td   = d(p)
         d(p) = d(i)
 18      tc   = c(p)
         c(p) = c(i)
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
c
c Start at the beginning of the segment, search for k such that a(k)>t
c
      q = j
      k = i
 20   k = k+1
      if(k.gt.q)     go to 60
      if(a(k).le.ta) go to 20
c
c Such an element has now been found now search for a q such that a(q)<t
c starting at the end of the segment.
c
 30   continue
      if(a(q).lt.ta) go to 40
      q = q-1
      if(q.gt.k)     go to 30
      go to 50
c
c a(q) has now been found. we interchange a(q) and a(k)
c
 40   xa   = a(k)
      a(k) = a(q)
      a(q) = xa
      go to (45,44,43,42,41,411,412,413),iring
 413     xh   = h(k)
         h(k) = h(q)
         h(q) = xh
 412     xg   = g(k)
         g(k) = g(q)
         g(q) = xg
 411     xf   = f(k)
         f(k) = f(q)
         f(q) = xf
 41      xe   = e(k)
         e(k) = e(q)
         e(q) = xe
 42      xd   = d(k)
         d(k) = d(q)
         d(q) = xd
 43      xc   = c(k)
         c(k) = c(q)
         c(q) = xc
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
c
c Update q and search for another pair to interchange:
c
      q = q-1
      go to 20
 50   q = k-1
 60   continue
c
c The upwards search has now met the downwards search:
c
      a(i)=a(q)
      a(q)=ta
      go to (65,64,63,62,61,611,612,613),iring
 613     h(i) = h(q)
         h(q) = th
 612     g(i) = g(q)
         g(q) = tg
 611     f(i) = f(q)
         f(q) = tf
 61      e(i) = e(q)
         e(q) = te
 62      d(i) = d(q)
         d(q) = td
 63      c(i) = c(q)
         c(q) = tc
 64      b(i) = b(q)
         b(q) = tb
 65   continue
c
c The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
c store the position of the largest segment in lt and ut
c
      if (2*q.le.i+j) go to 70
      lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 80
 70   lt(m) = q+1
      ut(m) = j
      j = q-1
c
c Update m and split the new smaller segment
c
 80   m = m+1
      go to 10
c
c We arrive here if the segment has  two elements we test to see if
c the segment is properly ordered if not, we perform an interchange
c
 90   continue
      if (a(i).le.a(j)) go to 100
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      go to (95,94,93,92,91,911,912,913),iring
 913     xh   = h(i)
         h(i) = h(j)
         h(j) = xh
 912     xg   = g(i)
         g(i) = g(j)
         g(j) = xg
 911     xf   = f(i)
         f(i) = f(j)
         f(j) = xf
   91    xe   = e(i)
         e(i) = e(j)
         e(j) = xe
   92    xd   = d(i)
         d(i) = d(j)
         d(j) = xd
   93    xc   = c(i)
         c(i) = c(j)
         c(j) = xc
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
c
c If lt and ut contain more segments to be sorted repeat process:
c
 100  m = m-1
      if (m.le.0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
 110  continue
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************

      real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c    Squared Anisotropic Distance Calculation Given Matrix Indicator
c    ***************************************************************
c
c This routine calculates the anisotropic distance between two points
c  given the coordinates of each point and a definition of the
c  anisotropy.
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         Coordinates of first point
c   x2,y2,z2         Coordinates of second point
c   ind              The rotation matrix to use
c   MAXROT           The maximum number of rotation matrices dimensioned
c   rotmat           The rotation matrices
c
c
c
c OUTPUT VARIABLES:
c
c   sqdist           The squared distance accounting for the anisotropy
c                      and the rotation of coordinates (if any).
c
c
c NO EXTERNAL REFERENCES
c
c
c-----------------------------------------------------------------------

C____________________ Declaration of variables


      IMPLICIT NONE

      INTEGER*4 IND,MAXROT

      INTEGER*4 I

      REAL*8 X1,Y1,Z1,X2,Y2,Z2

      REAL*8 CONT,DX,DY,DZ

C____________________ Declaration of arrays

      REAL*8 ROTMAT (MAXROT,3,3)
c
c Compute component distance vectors and the squared distance:
c
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist = 0.0d0
      do i=1,3
            cont   = rotmat(ind,i,1) * dx
     +             + rotmat(ind,i,2) * dy
     +             + rotmat(ind,i,3) * dz
            sqdist = sqdist + cont * cont
      end do
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************

      subroutine SEARCH_SAMPLES
     ;                   (xloc,yloc,zloc,radsqd,irot,MAXROT,rotmat,
     +                    nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct,nd,
     +                    x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup,
     +                    nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup,
     +                    nclose,close,infoct)
c-----------------------------------------------------------------------
c
c              Search Within Super Block Search Limits
c              ***************************************
c
c
c This subroutine searches through all the data that have been tagged in
c the super block subroutine.  The close data are passed back in the
c index array "close".  An octant search is allowed.
c
c
c
c INPUT VARIABLES:
c
c   xloc,yloc,zloc   location of point being estimated/simulated
c   radsqd           squared search radius
c   irot             index of the rotation matrix for searching
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c   nsbtosr          Number of super blocks to search
c   ixsbtosr         X offsets for super blocks to search
c   iysbtosr         Y offsets for super blocks to search
c   izsbtosr         Z offsets for super blocks to search
c   noct             If >0 then data will be partitioned into octants
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   tmp(nd)          Temporary storage to keep track of the squared
c                      distance associated with each data
c   nisb()                Array with cumulative number of data in each
c                           super block.
c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
c   nysup,ymnsup,ysizsup  Definition of the X super block grid
c   nzsup,zmnsup,zsizsup  Definition of the X super block grid
c
c
c
c OUTPUT VARIABLES:
c
c   nclose           Number of close data
c   close()          Index of close data
c   infoct           Number of informed octants (only computes if
c                      performing an octant search)
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c   sortem           Sorts multiple arrays in ascending order
c
c
c
c-----------------------------------------------------------------------

C_________________ Declaration of variables


      IMPLICIT NONE

      INTEGER*4 IROT,MAXROT,NSBTOSR,NOCT,ND,NXSUP,NYSUP,NZSUP,NCLOSE
     ;         ,INFOCT

      INTEGER*4 IX,IY,IZ,ISUP,IXSUP,IYSUP,IZSUP,II,NUMS,I,NT,NA,J,IQ

      REAL*8 XLOC,YLOC,ZLOC,RADSQD,XMNSUP,XSIZSUP,YMNSUP,YSIZSUP
     ;      ,ZMNSUP,ZSIZSUP


      REAL*8 C,D,E,F,G,H,DX,DY,DZ

C_________________ Declaration of arrays

      REAL*8    X(*),Y(*),Z(*),TMP(*),CLOSE(*)
      REAL*8    ROTMAT(MAXROT,3,3),HSQD,SQDIST
      INTEGER*4 NISB(*),INOCT(8)
      INTEGER*4 IXSBTOSR(*),IYSBTOSR(*),IZSBTOSR(*)
      LOGICAL INFLAG

c
c Determine the super block location of point being estimated:
c
      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
c
c Loop over all the possible Super Blocks:
c
      nclose = 0
      do 1 isup=1,nsbtosr
c
c Is this super block within the grid system:
c
            ixsup = ix + ixsbtosr(isup)
            iysup = iy + iysbtosr(isup)
            izsup = iz + izsbtosr(isup)
            if(ixsup.le.0.or.ixsup.gt.nxsup.or.
     +         iysup.le.0.or.iysup.gt.nysup.or.
     +         izsup.le.0.or.izsup.gt.nzsup) go to 1
c
c Figure out how many samples in this super block:
c
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
            if(ii.eq.1) then
                  nums = nisb(ii)
                  i    = 0
            else
                  nums = nisb(ii) - nisb(ii-1)
                  i    = nisb(ii-1)
            endif
c
c Loop over all the data in this super block:
c
            do 2 ii=1,nums
                  i = i + 1
c
c Check squared distance:
c
                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
     +                          MAXROT,rotmat)
                  if(real(hsqd).gt.radsqd) go to 2
c
c Accept this sample:
c
                  nclose = nclose + 1
                  close(nclose) = real(i)
                  tmp(nclose)  = real(hsqd)
 2          continue
 1    continue
c
c Sort the nearby samples by distance to point being estimated:
c
      call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)
c
c If we aren't doing an octant search then just return:
c
      if(noct.le.0) return
c
c PARTITION THE DATA INTO OCTANTS:
c
      do i=1,8
            inoct(i) = 0
      end do
c
c Now pick up the closest samples in each octant:
c
      nt = 8*noct
      na = 0
      do j=1,nclose
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
 6          continue
            inoct(iq) = inoct(iq) + 1
c
c Keep this sample if the maximum has not been exceeded:
c
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
      end do
c
c End of data selection. Compute number of informed octants and return:
c
 7    nclose = na
      infoct = 0
      do i=1,8
            if(inoct(i).gt.0) infoct = infoct + 1
      end do
c
c Finished:
c
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************
