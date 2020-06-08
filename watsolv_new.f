c********************************************************************
c  modified by Matthias Willmann,19.11.02
c
c  in order to incorporate the subroutine into a modular finite element 
c  flow code.
c  main changes: 1. the two common files ('main_dim.for' and 'watsolv_prm.for'
c  with variables and parameters are defined within the routines passed to 
c  others if necessary.
c  2. the use of maximum array sizes is replaced by paasing the actual 
c  array size (number of unknoens and bandwith) as arguments and its use.
c  (maxnn=nn,maxne=ne)
c======================================================================
c
c                            WATSOLV
c
c              Sparse Matrix Iterative Solver Package
c
c                          Version 1.01
c
c                   GMRES and CGSTAB Acceleration
c               Incomplete lower/upper factorization+
c                       for preconditioning
c                 Compressed banded data structure+
c
c                       Copyright (c) 1995
c
c                        J.E. VanderKwaak
c                          P.A. Forsyth*
c                        K.T.B. MacQuarrie**
c                          E.A. Sudicky
c
c             Waterloo Centre for Groundwater Research
c                      University of Waterloo
c                Waterloo, Ontario, Canada, N2L 3G1
c
c                * Department of Computer Science
c                      University of Waterloo
c
c               ** Department of Civil Engineering
c                    University of New Brunswick
c                         Fredriction, NB
c
c         + other options are available from the author(s)
c            direct inqeries to kwaak@galerkin.uwaterloo.ca
c----------------------------------------------------------------------
c COPYRIGHT NOTICE AND USAGE LIMITATIONS
c
c ALL RIGHTS ARE RESERVED; THE WATSOLV SUBROUTINES AND USER'S GUIDE
c ARE COPYRIGHT. THE DOCUMENTATION AND SOURCE CODE, OR ANY PART
c THEREOF, MAY NOT BE REPRODUCED, DUPLICATED, TRANSLATED, OR
c DISTRIBUTED IN ANY WAY WITHOUT THE EXPRESS WRITTEN PERMISSION OF
c THE COPYRIGHT HOLDER(S). PERMISSION IS GRANTED FOR THE USE OF THIS
c PACKAGE IN NON-PROFIT, SCHOLARLY RESEARCH ONLY: PAPERS OR REPORTS
c PRODUCED USING WATSOLV SHOULD EXPLICITLY ACKNOWLEDGE ITS SOURCE.
c THE WATSOLV PACKAGE, OR ANY OF ITS COMPONENT SUBROUTINES, MUST
c BE SPECIFICALLY LICENSED FOR INCLUSION IN SOFTWARE DISTRIBUTED IN ANY
c MANNER AND/OR SOLD COMMERCIALLY.
c
c----------------------------------------------------------------------
c DISCLAIMER
c
c Although great care has been taken in preparing the WATSOLV
c subroutines and documentation, the author(s) cannot be held
c responsible for any errors or omissions. As such, this code is
c offered `as is'. The author(s) makes no warranty of any kind,
c express or implied. The author(s) shall not be liable for any
c damages arising from a failure of this program to operate in the
c manner desired by the user. The author(s) shall not be liable for
c any damage to data or property which may be caused directly or
c indirectly by use of this program. In no event will the author(s)
c be liable for any damages, including, but not limited to, lost
c profits, lost savings or other incidental or consequential damages
c arising out of the use, or inability to use, this program. Use,
c attempted use, and/or installation of this program shall constitute
c implied acceptance of the above conditions. Authorized users
c encountering problems with the code, or requiring specific
c implementations not supported by this version, are encouraged to
c contact the author(s) for possible assistance.
c
c----------------------------------------------------------------------
c REVISION HISTORY
c
c Version 1.00 (02.05.95)
c   banded data structure, variable level symbolic factorization
c   version sent to C. Voss, USGS for use in SUTRA
c
c Version 1.01 (02.06.95)
c   limited CGSTAB restarts
c   catch level zero factorization
c   free, unlimited distribution
c   cannot be sold or included in commercial software
c
c Version 2.00 (02.17.95)
c   ia,ja data structure, variable level symbolic factorization
c   for sale, limited distribution
c
c Version 3.00 (00.00.00)
c   block ia,ja version
c   for sale, limited distribution
c
c----------------------------------------------------------------------
c Subroutine call summary (consult manual for details):
c
c  Generate node adjacency information:
c
c    Finite elements:
c
c          fe_iadm (in,nn,ne,nln)
c
c    Finite element mimicking finite difference:
c
c          mfe_iadm (in,nn,ne,nln)
c
c    Finite difference or finite volume:
c
c          fd_iadm (nx,ny,nz,nn)
c
c  Locate matrix position for a nodal pair (assembly):
c
c          find (i,j,k)
c
c  Perform symbolic factorization of coefficient matrix up to
c  specified level of fill:
c
c          symfac (nn,level)
c
c  Decompose global matrix for preconditioning:
c
c$          call w_factor (a,nn,maxnb,iad,iadn,iadd,af)
c  Solve the matrix equation using CGSTAB or GMRES:
c
c          watsolv (nn,a,b,x,north,nitmax,idetail,solve_flag,
c                   rtwotol,rmaxtol,smaxtol)
c
c  Include Files (see manual):
c
c  main_dim.for    - defines the array dimensioning parameters
c  main.prm    - arrays to be defined in the source program
c                  (passed to the solver subroutines)
c  watsolv_prm.for - defines solver array space/variables
c
c OUTPUT DIRECTED TO UNIT 6 (defaults to screen if unit not opened)
c----------------------------------------------------------------------
c WATSOLV Termination Possibilities             TEST
c
c            i
c    1(+) ||r ||   <= rmaxtol          absolute residual scale
c               max
c
c            i                   0
c    2(+) ||r ||   <= rtwotol*||r ||   residual scale wrt initial
c               2                   2
c
c    3(-) Happy Breakdown
c
c             i
c    4(+) ||dx ||   <= smaxtol         absolute update scale
c                max
c
c             i
c         ||dx ||
c                max
c    5(+) -------   <= smaxtol         relative update scale
c             i
c          ||x ||
c                max
c
c    6(+) Exceed Maximum Iterations
c
c    7(*) Anticipated divide by zero   trapped in driver
c    |
c    |(+) GMRES and CGSTAB
c    |(-) GMRES only
c    |(*) CSTAB only
c    +--  WATSOLV return code (iterm)
c----------------------------------------------------------------------
c LIBRARY ROUTINES USED
c
c$    w_cputime()                   COMPILER-SPECIFIC CLOCK TIME CALL
c    dcopy(n,dx,dy)    copies a vector, x, to a vector, y
c    ddot(n,dx,dy)     forms the dot product of two vectors
c    daxpy(n,da,dx,dy) constant times a vector plus a vector
c    dnrm2(n,dx,incx)            euclidean norm of dx()
c    dscal(n,da,dx)         scales a vector by a constant
c    idamax(n,dx,incx)           index of max. absolute value entry
c
c======================================================================
c Last Modified: February 1, 1995

      subroutine watsolv(nn,a,b,x,north,nitmax,idetail,solve_flag,
     {                   rtwotol,rmaxtol,smaxtol,iad,iadd,iadn,maxnb,
     &                   iafd,iafdd,iafdn,maxnfb,af,ITERM)


!      include 'main_dim.for'
!      include 'watsolv_prm.for'
c----------------------------------------------------------------------
c Passed Variables
c the variables added by matthias to incorporate the routines into a 
c modular finite element code
c res,af,maxit,mnln,mnorth,numit,iterm,converged,rtnorm,rmnorm,rtcheck,rinit,dx,dxmnorm
c xmnorm,xratio
      real*8
     {                 a(maxnb,nn),        ! stiffness matrix
     {                 x(nn),              ! approx solution         **
     {                 b(nn),              ! forcing vector
     {                 rtwotol,            ! residual two-norm tolerance
     {                 rmaxtol,            ! residual max-norm tolerance
     {                 smaxtol             ! approx sol norm tolerance
      integer*4          maxnb,
     {                 nn,                 ! number of unknowns
     {                 north,              ! max num of orth
     {                 nitmax,             ! max num of iterations
     {                 idetail             ! echo solver info to file
      logical
     {                 solve_flag          ! solution method flag

c                                               modified on return = **
c----------------------------------------------------------------------
c Local variables

      real*8
     {                 time   !$ timing info

      integer*4
     {                 irs                 ! CGSTAB restart counter
c----------------------------------------------------------------------
c matthias: new variables formerly passed throuh common files
      integer*4 iad(maxnb,nn),
     &          iadd(nn),
     &          iadn(nn),
     &          iafd(maxnfb,nn),
     &          iafdd(nn),
     &          iafdn(nn)
      real*8, dimension(:),allocatable:: Res(:)
      real*8 af(maxnfb,nn)
      integer*4 mnorth,numit,iterm,maxnfb
      logical converged
      real*8 rtnorm,rmnorm,rtcheck,rinit,dx,dxmnorm,
     &       xmnorm,xratio
c----------------------------------------------------------------------
      real*8 dnrm2
      integer*4 idamax
      external dcopy,dnrm2,idamax



c----------------------------------------------------------------------
Cmatthias--- set mnorth 10 
      mnorth=25     

c transdens
      allocate(res(nn))
c end transdens


      if (solve_flag.and.(north.gt.mnorth)) then
        write (6,*) 'reset mnorth and recompile (north > mnorth)'
        write (6,*) 'north  = ',north
        write (6,*) 'mnorth = ',mnorth
        stop
      end if

      if (idetail.gt.0) then
c       time = w_cputime() !$
        if (solve_flag) then
          write (6,9) 'GMRES '
          write (6,12) north
        else
          write (6,9) 'CGSTAB'
        end if
      end if

      i = idamax(nn,x(1))
      if (dabs(x(i)).gt.1d-20) then      ! catch initial guess of zero
        call mvmult (res,a,x,nn,iad,iadn,maxnb)
        do i = 1,nn

c matthias 15.11. if statement added to prevent rounding errors
c   with fixed head value (*1.0d20) which are only exact to the 
c   8th value        
          if(b(i).gt.1.0d17) then
            res(i)=0.0
          else  
            res(i) = b(i) - res(i)
          end if
        end do
      else
        call dcopy(nn,b(1),res(1))
      end if

c factor can be called from the main program or from within watsolv.
c if the stiffness matrix doesn't change, decomposition isn't needed.
c
c      call w_factor(a,nn,maxnb,iad,iadn,af,iafd,iafdd,iafdn,maxnfb)                  ! factor stiffness matrix

      numit = 0
      irs   = 0

100   continue                           ! restart location (cgstab)

      rtnorm    = dnrm2(nn,res,1)
      rtcheck   = rtwotol*rtnorm         ! residual conv. criteria
      i         = idamax(nn,res)
      rmnorm    = dabs(res(i))
      rinit     = rtnorm
      dxmnorm   = 0.0
      xratio    = 0.0
      iterm     = 1
      converged = .false.

      if (idetail.gt.0) write (6,7) rtnorm,rmnorm
      if (idetail.gt.1) write (6,10) numit,rtnorm,1.0d0

      if (dabs(rtnorm).gt.0.0) then    ! null initial residual
        if (solve_flag) then
          call gmres(nn,a,x,nitmax,north,rmaxtol,smaxtol,idetail,
     &    maxnb,iad,iadn,res,af,mnorth,numit,iterm,converged,
     &    rtnorm,rmnorm,rtcheck,rinit,dxmnorm,xmnorm,xratio,
     &    iafd,iafdd,iafdn,maxnfb)
        else
          call cgstab(nn,a,x,nitmax,rmaxtol,smaxtol,idetail,maxnb,
     &    res,af,numit,iterm,converged,rmnorm,rtcheck,rinit,dx,
     &    dxmnorm,xmnorm,xratio,iad,iadd,iadn,iafd,iafdd,iafdn,maxnfb)
 
          if (iterm.gt.6) then           ! anticipated divide by zero
            irs = irs + 1
            if (irs.gt.20) then          ! max restarts = 20
              write (16,13)
              stop
            else
              call mvmult (res,a,x,nn,iad,iadn,maxnb)
              do i = 1,nn


c matthias 15.11. if statement added to prevent rounding errors
c   with fixed head value (*1.0d20) which are only exact to the 
c   8th value 
                if(b(i).gt.1.0d17) then
                  res(i)=0.0
                else  
                  res(i) = b(i) - res(i)
                end if
              end do
              goto 100                   ! restart cgstab
            end if
          end if
        end if
      end if

      if ((numit.eq.nitmax).and..not.converged) iterm = 6

      if (idetail.gt.0) then
c       time = w_cputime() - time !$
        write (6,8) time,numit,iterm,rtnorm,rtcheck,rmnorm,rmaxtol,
     {               dxmnorm,smaxtol,xratio,smaxtol
        if (iterm.eq.1) then
          write (6,1)
        else if (iterm.eq.2) then
          write (6,2)
        else if (iterm.eq.3) then
          write (6,3)
        else if (iterm.eq.4) then
          write (6,4)
        else if (iterm.eq.5) then
          write (6,5)
        else
          write (6,6)
        end if
      end if

    1 format (/5x,'convergence due to residual max-norm'/)
    2 format (/5x,'convergence due to residual 2-norm'/)
    3 format (/5x,'convergence due to happy breakdown'/)
    4 format (/5x,'convergence due to absolute solution update scale'/)
    5 format (/5x,'convergence due to relative solution update scale'/)
    6 format (/5x,'maximum iterations exceeded'/)
    7 format (/5x,'initial residual 2-norm---',1pd12.5,
     {        /5x,'initial residual max-norm-',1pd12.5)
    8 format (/5x,'elapsed time-----------------',1pd12.5,
     {        /5x,'number of iterations---------',i5,
     {        /5x,'termination code-------------',i5,
     {        /5x,'final residual 2-norm--------',1pd12.5,
     {        /5x,'residual 2-norm tolerance----',1pd12.5,
     {        /5x,'residual maximum norm--------',1pd12.5,
     {        /5x,'residual max norm tolerance--',1pd12.5,
     {        /5x,'solution update max norm-----',1pd12.5,
     {        /5x,'solution update tolerance----',1pd12.5,
     {        /5x,'update/solution ratio--------',1pd12.5,
     {        /5x,'update/solution tolerance----',1pd12.5)
    9 format (//3x,a6,' Convergence Information')
   10 format (/5x,'iteration   residual      rcurr/rinit',
     {        /5x,i5,2(4x,1pd12.5))
c   11 format (/'Anticipated divide by zero in CGSTAB'
c     {        /'   Restarting Iterations'/)
   12 format (/5x,'Retaining',i3,' basis vectors')
   13 format (//5x,'MAXIMUM CGSTAB RESTART EXCEEDED',
     {         /5x,'      STOPPING WATSOLV',
     {         /5x,' CHECK PROBLEM FORMULATION')
      
      deallocate(res)
       
      return
      end

c======================================================================
c
c NOTE: COMPILER-SPECIFIC CALL IN CPUTIME FUNCTION (real*8)
c
c======================================================================
c Last Modified: February 1, 1995

c      function w_cputime() !$
      !USE DFPORT
c     real cputime
c      real*8 w_cputime !$
c     integer*4 mclock
c
c XLF (IBM RS6000) function call
c
c      cputime = dfloat( mclock() ) / 100.0d0   !seconds
c      cputime = dfloat( mclock() ) / 360000d0  !hours
c
c SALFORD (IBM PC) function call
c
c      call dclock@(cputime)                     !seconds
c     cputime = cputime / 3600d0                !hours
c
c VISUAL FORTRAN (IBM PC)$$

c       REAL*4 TARRAY(2)

c       w_cputime=DTIME(TARRAY)

c      return
c      end

c======================================================================
c
c accelerated iterative solution routines
c
c======================================================================
c Last Modified: February 1, 1995

      subroutine cgstab (nn,a,x,nitmax,rmaxtol,smaxtol,idetail,maxnb,
     &    res,af,numit,iterm,converged,rmnorm,rtcheck,rinit,dx,
     &    dxmnorm,xmnorm,xratio,iad,iadd,iadn,iafd,iafdd,iafdn,maxnfb)

!      include 'main_dim.for'
!      include 'watsolv_prm.for'

      real*8 a(maxnb,nn),x(nn),rmaxtol,smaxtol
      integer*4 nn,idetail,nitmax,maxnfb
c----------------------------------------------------------------------
c Local Variables

      real*8
     {                 alpha,              ! solution change coefficient
     {                 rholst,             ! last rho value
     {                 omega,              ! last omega value
     {                 bbeta,
     {                 rho,                ! current beta value
     {                 step1,              ! solution update
     {                 step2,              !  variables
     {                 rres                ! residual (temp)
      real*8, dimension(:), allocatable::
     {                 res0(:),        ! initial residual vector
     {                 pvec(:),        ! intermediate arrays
     {                 vbar(:),
     {                 avbar(:),
     {                 svec(:),
     {                 zvec(:),
     {                 tvec(:) 

c----------------------------------------------------------------------
c matthias: new variables formerly passed throuh common files
      integer*4 iad(maxnb,nn),
     &          iadd(nn),
     &          iadn(nn),
     &          iafd(maxnfb,nn),
     &          iafdd(nn),
     &          iafdn(nn)
      real*8 res(nn),tiny
      real*8 af(maxnfb,nn)
      integer*4 numit,iterm
      logical converged
      real*8 rtnorm,rmnorm,rtcheck,rinit,dx,dxmnorm
      real*8 xmnorm,xratio
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      real*8 ddot
      external ddot

c----------------------------------------------------------------------
c new transdens : allocate internal arrays
 
      allocate(res0(nn))        ! initial residual vector
      allocate(pvec(nn))        ! intermediate arrays
      allocate(vbar(nn))
      allocate(avbar(nn))
      allocate(svec(nn))
      allocate(zvec(nn))
      allocate(tvec(nn))

c end new transdens

c matthias, set tiny
      tiny=1.0d-300
c initialize variables and vectors

      alpha  = 1.0
      rholst = 1.0
      omega  = 1.0
      do i = 1,nn
        pvec(i)  = 0.0
        avbar(i) = 0.0
        res0(i)  = res(i)   ! save initial residual
      end do

c----
c BEGIN ITERATION LOOP (up to nitmax)

      do while (.not. converged            ! residual tolerance check
     +          .and. numit.lt.nitmax)     ! max iterations
        numit = numit + 1

c------------------------------
c Phases of CG-Stab acceleration
c
c rho    = (res0, res)
c bbeta  = (rho/rholst)*(alpha/omega)
c rholst = rho
c pvec   = res + bbeta*(pvec - omega*Avbar)

        rho = ddot(nn,res0,res)
        if (.not.(dabs(rho).gt.tiny)) then ! anticipated divide by zero
          iterm = 7
          return
        end if
        bbeta = rho/ (rholst + sign(tiny,rholst))
        bbeta = bbeta*(alpha/ (omega + sign(tiny,omega)))
        rholst = rho
        do i = 1,nn
          pvec(i) = res(i) + bbeta* (pvec(i) - omega*avbar(i))
        end do

c--------
c solve (LU) vbar = pvec
c Avbar = A*vbar
c alpha = rholst/(res0,Avbar)

        call lusolv (vbar,pvec,nn,af,
     &       iafd,iafdd,iafdn,maxnfb)
        call mvmult (avbar,a,vbar,nn,iad,iadn,maxnb)

        alpha = ddot(nn,res0,avbar)
        if (.not.(dabs(alpha).gt.tiny)) then  ! divide by zero
          iterm = 7
          return
        end if
        alpha = rho/(alpha + sign(tiny,alpha))

c--------
c svec  = res - alpha*Avbar
c solve  (LU) zvec = svec
c Azvec = A*zvec = tvec
c omega = (A zvec,svec)/(A zvec,A zvec) in CG-STAB, but
c omega = (tvec,svec)/(tvec,tvec) in CG-STAB-P

        do i = 1,nn
          svec(i) = res(i) - alpha*avbar(i)
        end do

        call lusolv (zvec,svec,nn,af,
     &       iafd,iafdd,iafdn,maxnfb)
        call mvmult (tvec,a,zvec,nn,iad,iadn,maxnb)

        omega = ddot(nn,tvec,tvec)
        if (.not.(dabs(omega).gt.tiny)) then  ! divide by zero
          iterm = 7
          return
        end if
        omega = ddot(nn,tvec,svec)/(omega + sign(tiny,omega))

c--------
c x = x + alpha*vbar + omega*zvec
c res = svec - omega*tvec
c
        rtnorm  = 0.0
        rmnorm  = 0.0
        dxmnorm = 0.0
        xmnorm  = 0.0
        do i = 1,nn
          step1   = alpha*vbar(i)
          step2   = omega*zvec(i)
          dx      = step1 + step2
          dxmnorm = dmax1(dabs(dx),dxmnorm)         ! update max norm
          x(i)  = x(i) + dx
          xmnorm  = dmax1(dabs(x(i)),xmnorm)        ! solution max norm
          rres    = svec(i) - omega*tvec(i)
          res(i)  = rres
          rtnorm  = rtnorm + rres*rres              ! residual 2-norm
          rmnorm  = dmax1(dabs(rres),rmnorm)        ! residual max norm
        end do
        

         
        rtnorm = dsqrt(rtnorm)                      ! residual 2-norm
        xratio = dxmnorm/xmnorm

c convergence checks

        converged = .not.(rmnorm.gt.rmaxtol)        ! check (1)
        if (.not.converged) then
          converged = (.not.(rtnorm.gt.rtcheck))    ! check (2)
          if (converged) then
            iterm = 2
          else
            call check_update(smaxtol,converged,dxmnorm,iterm,xratio) ! check (4,5)
          end if
        end if
        if (idetail.gt.1) write (6,2) numit,rtnorm,rtnorm/rinit
      end do                           ! iteration loop

2     format (5x,i5,2(4x,1pd12.5))

c new transdens : deallocate internal arrays
 
      deallocate(res0,pvec,vbar,avbar,svec,zvec,tvec)

c end new transdens




      return
      end

c======================================================================
c Last Modified: February 1, 1995

      subroutine gmres(nn,a,x,nitmax,north,rmaxtol,smaxtol,idetail,
     &    maxnb,iad,iadn,res,af,mnorth,numit,iterm,converged,
     &    rtnorm,rmnorm,rtcheck,rinit,dxmnorm,xmnorm,xratio,
     &    iafd,iafdd,iafdn,maxnfb)

c     include 'main_dim.for'
c     include 'watsolv_prm.for'

c----------------------------------------------------------------------
c change matthias: new variables formerly passed throuh common files
      integer*4 iad(maxnb,nn),iadn(nn),iafd(maxnfb,nn)
     &         ,iafdd(nn), iafdn(nn)
      real*8 res(nn)
      real*8 af(maxnfb,nn)
      integer*4 mnorth,numit,iterm,maxnfb
      logical converged
      real*8 rtnorm,rmnorm,rtcheck,rinit,dxmnorm,
     &       xmnorm,xratio
c----------------------------------------------------------------------

      real*8 a(maxnb,nn),x(nn),rmaxtol,smaxtol
      integer*4 nn,north,idetail,nitmax
c----------------------------------------------------------------------
c Local Variables

      real*8,allocatable, dimension(:)::
     {                 work(:),        ! temp work vector
     {                 hb(:)       ! hess. right hand side

      real*8,allocatable, dimension(:,:)::
     {                 v(:,:),  ! basis vectors
     {                 h(:,:), ! factored hessenburg matrix
     {                 cs(:,:)      ! fact. param. (rotations)

      real*8           c,                  ! cos rotation
     {                 s,                  ! sin rotation
     {                 t1,t2,t3            ! misc. temp. var.
      integer*4
     {                 nk,                 ! factorization counter
     {                 nkp1,               ! fact. counter plus one
     {                 i,j,k               ! misc. index var
c----------------------------------------------------------------------
      real*8 ddot,dnrm2
      integer*4 idamax
      external dcopy,ddot,daxpy,dnrm2,dscal,idamax
c----------------------------------------------------------------------

c new transdens: allocatable local arrays
      allocate(work(nn),hb(mnorth+1) )
      allocate(v(nn,mnorth+1),h(mnorth+1,mnorth), cs(2,mnorth+1))
c  end new transdens

c initialize

      call dcopy(nn,v(1,1),1)

c----------
c BEGIN OUTER ITERATION LOOP (up to nitmax)
c
1     continue                             ! iteration loop (init/restart)
      nk   = 0
      nkp1 = 1
      hb(1) = rtnorm                    ! init right hand side
      t1 = 1.0/rtnorm
      call dscal(nn,t1,v(1,1))
c----
c BEGIN INNER ITERATION LOOP (up to north)

      do while (.not. converged            ! convergence check
     {          .and. numit.lt.nitmax      ! max iterations
     {          .and. nk.lt.north)         ! max factorizations
        nk    = nkp1                       ! fact. counter
        nkp1  = nk + 1
        numit = numit + 1               ! total iterations

c LU solve
        call lusolv (work,v(1,nk),nn,af,
     &       iafd,iafdd,iafdn,maxnfb)
        call mvmult (v(1,nkp1),a,work,nn,iad,iadn,maxnb)

c form basis vector - orthgonalize

        do i = 1,nk
          t1 = ddot(nn,v(1,i),v(1,nkp1))
          h(i,nk) = t1
          call daxpy(nn,-t1,v(1,i),v(1,nkp1))
        end do
        t1 = dnrm2(nn,v(1,nkp1),1)   ! scale the new vector
        h(nkp1,nk) = t1

c test for breakdown (3)

        if (.not.(t1.gt.0.0)) then
          converged = .true.               ! norm of current basis
          iterm = 3                        ! is zero, so we have
          nkp1 = nk                        ! the solution
          nk   = nk - 1
        else
          t1 = 1.0/t1
          call dscal(nn,t1,v(1,nkp1))

c apply previous rotations on new column of h

          if (nk.gt.1) then
            do i = 1,nk - 1
              ip1 = i + 1
              c  = cs(1,i)
              s  = cs(2,i)
              t1 = h(i,nk)
              t2 = h(ip1,nk)
              h(i,nk)   = c*t1 + s*t2
              h(ip1,nk) = c*t2 - s*t1
            end do
          end if

c determine new rotations

          t1 = h(nk,nk)
          t2 = h(nkp1,nk)
          if (.not.(dabs(t2).gt.0.0)) then
            c = 1.0
            s = 0.0
          else if (dabs(t2).gt.dabs(t1)) then
            t3 = t1/t2
            s = 1.0/dsqrt(1.0 + t3*t3)
            c = s*t3
          else
            t3 = t2/t1
            c = 1.0/dsqrt(1.0 + t3*t3)
            s = c*t3
          end if

c   complete the elimination and apply to hb (right hand side)
c   norm of current residual vector is now available as last
c   entry in the right hand side vector

          h(nk,nk) = c*t1 + s*t2
          hb(nkp1) = -s*hb(nk)
          hb(nk)   =  c*hb(nk)
          rtnorm   = dabs(hb(nkp1))
          cs(1,nk) = c                     ! save rotations
          cs(2,nk) = s

c test residual norm relative to initial resitual (2)

          converged = (.not.(rtnorm.gt.rtcheck))
          if (converged) iterm = 2
          if (idetail.gt.1) write (6,2) numit,rtnorm,rtnorm/rinit
        end if
      end do

c END OF INNER ITERATION LOOP
c----
c solve upper triangular system

      hb(nk) = hb(nk)/h(nk,nk)
      if (nk.gt.1) then
        do i = 2,nk
          k   = nkp1 - i
          kp1 = k + 1
          t1 = hb(k)
          do j = kp1,nk
            t1 = t1 - h(k,j)*hb(j)
          end do
          hb(k) = t1/h(k,k)
        end do
      end if

c form solution update right hand side

      do i = 1,nn
        work(i) = v(i,1)*hb(1)
      end do
      do i = 2,nk
        call daxpy(nn,hb(i),v(1,i),work)
      end do

c solve for new solution update (apply precondit1r)

      call lusolv (work,work,nn,af,iafd,
     &             iafdd,iafdn,maxnfb)

      dxmnorm = 0.0
      xmnorm  = 0.0
      do i = 1,nn
        x(i) = x(i) + work(i)                  ! update approx. solution
        dxmnorm = dmax1(dabs(work(i)),dxmnorm) ! max update
        xmnorm  = dmax1(dabs(x(i)),xmnorm)     ! max solution
      end do
      xratio = dxmnorm/xmnorm

c test scale of solution update (4 and 5)

      call check_update(smaxtol,converged,dxmnorm,iterm,xratio)

c build restart residual vector if needed

      if (.not.converged                   ! convergence check
     {   .and. numit.lt.nitmax) then       ! max iteration check
        do i = 1,nk
          j = nkp1 - i
          c = cs(1,j)
          s = cs(2,j)
          jp1 = j + 1
          hb(j)   = -s*hb(jp1)
          hb(jp1) =  c*hb(jp1)
        end do
        t1 = hb(1) - 1.0
        call daxpy(nn,t1,v(1,1),v(1,1))
        do i = 2,nkp1
          call daxpy(nn,hb(i),v(1,i),v(1,1))
        end do

c test residual max norm (1)

        i = idamax(nn,v(1,1))
        rmnorm = dabs(v(i,1))
        converged = (.not.(rmnorm.gt.rmaxtol))
        if (converged) then
          iterm = 1
        else
          goto 1                      ! restart orthogonalizations
        end if
      end if

c END OF OUTER ITERATION LOOP
c----------

c new transdens: allocatable local arrays
      deallocate(work,hb )
      deallocate(v,h, cs)
c  end new transdens


2     format (5x,i5,2(4x,1pd12.5))
      return
      end

c======================================================================
c Last Modified: February 1, 1995

      subroutine check_update(smaxtol,converged,dxmnorm,iterm,xratio)

!      include 'main_dim.for'
!      include 'watsolv_prm.for'
      real*8 smaxtol,dxmnorm,xratio
      logical converged
      integer*4 iterm

c check solution and solution update scale (maximum value)

      if (.not.converged) then
        if (dxmnorm.gt.xratio) then
          converged = (.not.(dxmnorm.gt.smaxtol))
          if (converged) iterm = 4
        else
          converged = (.not.(xratio.gt.smaxtol))
          if (converged) iterm = 5
        end if
      end if
      return
      end

c======================================================================
c
c matrix-vector manipulation routines
c
c======================================================================
c Last Modified: February 1, 1995

      subroutine lusolv (x,b,nn,af,iafd,iafdd,
     &       iafdn,maxnfb)

!      include 'main_dim.for'
!      include 'watsolv_prm.for'

      real*8 x(nn),b(nn),af(maxnfb,nn)
      integer*4 nn,i,j,k,n,maxnfb
      integer*4 iafdn(nn),iafdd(nn),
     &          iafd(maxnfb,nn)
c Lower triangular matrix inversion by forward substitution and upper 
c triangular matrix inversion by backward substitution.
c Lower and upper triangular matrices are in af, and right-hand-side
c vector is in b at start. Solution vector is in x upon exit.

c forward solve:  Lz = b (L has unit diagonal)

      do i = 1,nn
        x(i) = b(i)
        k = iafdd(i)                           ! diagonal entry
        do j = 1,k - 1
          inode = iafd(j,i)                    ! connection
          x(i) = x(i) - af(j,i)*x(inode)
        end do
      end do

c backward solve: Ux = z (U does not have unit diag)

      do i = nn,1,-1
        n = iafdn(i)                            ! num connections
        k = iafdd(i)                            ! diagonal entry
        do j = (k + 1),n
          inode = iafd(j,i)                     ! connection
          x(i) = x(i) - af(j,i)*x(inode)
        end do
        x(i) = x(i)/af(k,i)                     ! diagonal entry
      end do
      return
      end

c======================================================================
c Last Modified: February 11, 1995

      subroutine mvmult (b,a,x,nn,iad,iadn,maxnb)

c Multiply matrix a by vector x to obtain b (Ax = b)

!      include 'main_dim.for'

      real*8 x(nn),b(nn),a(maxnb,nn)
      integer*4 iad(maxnb,nn),iadn(nn)

      do i = 1,nn
        sum = 0.0
        n = iadn(i)                           ! num connections
        do j = 1,n
          inode = iad(j,i)                    ! connection
          sum = sum + a(j,i)*x(inode)
        end do
        b(i) = sum
      end do
      return
      end

c======================================================================
c
c preconditioning routines - level of fill symbolic factorization
c
c======================================================================
c Last Modified: February 1, 1995
            
      subroutine w_factor (a,nn,maxnb,iad,iadn,af,iafd,iafdd,
     &                     iafdn,maxnfb)

c Incomplete lower-upper decomposition of matrix a into af.
c This is general code, organized for arbitrary fill level.
c Diagonal dominance is assumed: no pivoting performed.

!      include 'main_dim.for'
!      include 'watsolv_prm.for'

      real*8 a(maxnb,nn),mult
      real*8,allocatable,dimension(:):: row
      integer*4, allocatable,dimension(:):: list
      integer*4 iad(maxnb,nn),iadn(nn),iafdn(nn),
     &          iafdd(nn),iafd(maxnfb,nn)
      real*8 af(maxnfb,nn)
      integer*4 maxnfb,maxnb,nn

c transdens
      allocate (list(nn))
      allocate(row(nn))
c end transdens

c Initialize

      do i = 1,nn
        row(i)  = 0.0
        list(i) = 0
        do j = 1,maxnfb
          af(j,i) = 0.0
        end do
      end do

c factor

      do i = 1,nn
        do ii = 1,iadn(i)                 ! row entries in a()
          j = iad(ii,i)                      ! location in row
          row(j)  = a(ii,i)                  ! scatter row entries
        end do
        do ii = 1,iafdn(i)                ! row entries in af()
          j = iafd(ii,i)                     ! location in row
          list(j) = i                        ! mark fill locations
        end do
        do ii  = 1,iafdd(i) - 1        ! row entry loop
          id   = iafd(ii,i)                  ! inducing diagonal loc.
          mult = row(id)/af(iafdd(id),id)    ! inducing term
          row(id) = mult                     ! lower fill entry
          do iii = iafdd(id) + 1, iafdn(id)  ! row fill entries
            idd  = iafd(iii,id)              ! fill column location
            if (list(idd).gt.0)          ! include if fill location
     {        row(idd) = row(idd) - mult*af(iii,id)
          end do
        end do
        do ii = 1,iafdn(i)                ! row entries in af()
          j = iafd(ii,i)                     ! location in row
          af(ii,i) = row(j)                  ! gather row entries
          row(j)   = 0.0                   ! zero row entry
          list(j)  = 0                   ! null fill location
        end do
      end do


c transdens
      deallocate(row)
      deallocate(list)
c end transdens

      return
      end

c======================================================================
c Last Modified: February 11, 1995

      subroutine symfac(nn,level,maxnb,iad,iadd,iadn,iafd,iafdd,
     &                  iafdn,maxnfb)

c incomplete LU - brute force factor

!      include 'main_dim.for'
!      include 'watsolv_prm.for'

      integer*4 nn,level,maxnfb,maxnb
      integer*4 iad(maxnb,nn),iadd(nn),iadn(nn),iafd(maxnfb,nn),
     &          iafdd(nn),iafdn(nn)

c local variables
      integer*4, allocatable, dimension(:):: list, lrow
      integer*4, allocatable, dimension(:,:):: levptr
      integer*4 ir,itemp,ii,first,num,next,nnp1,
     {        oldlst,levtemp,nxtlst,irow,jcol

c matthias, set maxint          
           maxint=999999
           ifbmax = 0          ! max bandwidth

c transdens
      allocate(list(nn),lrow(nn),levptr(maxnfb,nn))
c end transdens

      if (level.gt.0) then
        nnp1  = nn + 1       ! end of list marker
        do i = 1,nn
          lrow(i) = maxint
          list(i) = 0
          iafdd(i) = 0      ! factored diagonal location
          iafdn(i) = 0      ! number factored entries in row
          do j = 1,maxnfb
            iafd(j,i) = 0  ! columns
          end do
        end do

        do ir = 1,nn         ! loop through rows
          num = iadn(ir)        ! number of original entries
          do ii = 1,num      ! load row ir of L/U into list
            iafd(ii,ir) = iad(ii,ir)
            j = iad(ii,ir)
            lrow(j) = 0     ! initial level is zero
          end do

          first = iafd(1,ir) ! build linked list
          do ii = 1,num - 1
            list(iafd(ii,ir)) = iafd(ii + 1,ir)
          end do
          list(iafd(num,ir)) = nnp1        ! end of list flag (nn + 1)

          next = first                     ! first entry in linked list
          do while (next.lt.ir)
            oldlst = next
            nxtlst = list(next)            ! next column
            irow   = next                  ! current col/row
            do ii = iafdd(irow),iafdn(irow)! scan row 'irow' of U
              jcol = iafd(ii,irow)
1             continue
              if (jcol.lt.nxtlst) then   ! min fill level of new entry
                levtemp = levptr(ii,irow) + lrow(next) + 1
                levtemp = min0(lrow(jcol),levtemp)
                if (levtemp.le.level) then ! entry <= max level
                  list(oldlst) = jcol    ! add index to list()
                  list(jcol)   = nxtlst  ! row if new level smaller than
                  oldlst       = jcol    ! current level
                  lrow(oldlst) = levtemp
                end if
              elseif (jcol.eq.nxtlst) then! entry in (i,j) and (irow,j)
                oldlst  = nxtlst
                nxtlst  = list(oldlst)
              elseif (jcol.gt.nxtlst) then! entry in (i,j), not (irow,j)
                oldlst = nxtlst
                nxtlst = list(oldlst)
                go to 1                   ! check next entry
              endif
            end do
            next = list(next)             ! next column in linked list
          end do

          itemp = 0
          next = first
          do while (next.lt.nnp1)         ! linked list loop (gather)
            itemp = itemp + 1
            if (itemp.gt.maxnfb) then     ! allocation error
              write (6,2) itemp
              stop
            endif
            iafd(itemp,ir)   = next       ! column
            levptr(itemp,ir) = lrow(next) ! save level
            lrow(next) = maxint           ! reset level
            if (next.eq.ir)
     {        iafdd(ir) = itemp           ! diagonal entry
            next = list(next)             ! next entry in list
          end do
          iafdn(ir) = itemp               ! num entries in row
          ifbmax = max0(ifbmax,itemp)     ! max bandwidth
          if (iafdd(ir).eq.0) then    ! check diagonal
            write(6,*)' no diag in L/U'
            stop
          endif
        end do
      else             ! level zero = only original entries in a
        do i = 1,nn
          iafdd(i) = iadd(i)      ! factored diagonal location
          itemp    = iadn(i)
          iafdn(i) = itemp        ! number factored entries in row
          ifbmax = max0(ifbmax,itemp)     ! max bandwidth
          do j = 1,itemp
            iafd(j,i) = iad(j,i)  ! columns
          end do
        end do
      end if
      write (6,3) level,ifbmax

2     format (/5x,'Factorization bandwidth exceeded'
     {        /5x,' current bandwidth: ',i5)
3     format (/5x,'Symbolic Factorization Level: ',i5,
     {       /5x,'Resulting Maximum Bandwidth:  ',i5)

c transdens
      deallocate(list,lrow,levptr)
c end transdens


      return
      end
c======================================================================
c
c adjacency matrix routines
c
c======================================================================
c Last Modified: February 2, 1995

      subroutine fe_iadm (in,nn,ne,nln,iadd,iad,iadn,maxnb,LNNDEL)

c Generate the adjacency matrix for nodes from a finite element
c incidence matrix. Requires subroutine insert

!      include 'main_dim.for'
!      include 'watsolv_prm.for'

      integer*4 in(nln,ne),iadd(nn),iad(maxnb,nn),iadn(nn),LNNDEL(NE)

c Determine independent adjacency within each element
c - node ie2 is adjacent to node ie1
c - node ie1 is adjacent to node ie2

      do n = 1,ne                    ! loop through ELEMENTS
        NNUD = LNNDEL(N)
        do i = 1,NNUD -1          !(nln - 1)
          do j = (i + 1),NNUD     !NLN
            ie1 = in(i,n)
            ie2 = in(j,n)
            call insert (ie1,ie2,k,iadn,iad,maxnb,nn)
            call insert (ie2,ie1,k,iadn,iad,maxnb,nn)
          end do
        end do
      end do

c Determine self-adjacency terms - store in iadd()

      do i = 1,nn                     ! loop through NODES
        call insert (i,i,k,iadn,iad,maxnb,nn)
        iadd(i) = k
      end do
      return
      end


c======================================================================
c Last Modified: February 1, 1995

      subroutine insert (i,j,k,iadn,iad,maxnb,nn)

c Add j to the adjacency list for i
c Returns the position k where it has been added, or where it was
c already in the list.

!      include 'main_dim.for'
!      include 'watsolv_prm.for'

      integer*4 iad(maxnb,nn),iadn(nn)
c Determine number of nodes already in adjacency list

      n = iadn(i)
      k = n + 1

c Determine whether already in list

      do l = 1,n
        inode = iad(l,i)
        if (inode.ge.j) then
          k = l
          if (inode.eq.j) return
          go to 15
        end if
      end do

   15 continue

c Place in list (numerical order)

      if ((n + 1).gt.maxnb) then
        write (*,1) i,maxnb
        write (6,1) i,maxnb
        stop
      end if

      iadn(i) = n + 1
      do l = (n + 1),(k + 1),(-1)
        iad(l,i) = iad(l - 1,i)
      end do
      iad(k,i) = j

1     format (//5x,'error in iadmake: node ',i5,' has > '
     {         ,i5,' adjacencies')
      return
      end

c======================================================================
c Last Modified: February 1, 1995

      subroutine find (i,j,k,iad,iadn,maxnb,nn)

c For node i, determine the 'band' (k) related to its adjacency to
c node j. If node not adjacent, return 0 as the 'band'

!      include 'main_dim.for'
!      include 'watsolv_prm.for'
      integer*4 iad(maxnb,nn),iadn(nn)

      k = 0
      n = iadn(i)

      do l = 1,n
        inode = iad(l,i)

c Exit the loop if at or past the required position

        if (inode.ge.j) then
          if (inode.eq.j) k = l
          go to 20
        end if
      end do

   20 continue
      return
      end

c======================================================================
c
c Following routines are consistent with BLAS1\ESSL\LAPACK
c  and could be replaced with library calls.
c
c======================================================================
c Last Modified: February 11, 1995

      subroutine dcopy(n,dx,dy)

c copies a vector, x, to a vector, y.

      real*8 DX(N),DY(N)
    !  real*8 dx(*),dy(*)
c     real*8 dx(1),dy(1)
      integer*4 i,n

      do i = 1,n
        dy(i) = dx(i)
      end do
      return
      end

c======================================================================
c Last Modified: February 11, 1995

      real*8 function ddot(n,dx,dy)

c forms the dot product of two vectors.

      real*8 dx(*),dy(*)
c     real*8 dx(1),dy(1),dtemp
      integer*4 i,n

      ddot = 0.0d0
      do i = 1,n
        ddot = ddot + dx(i)*dy(i)
      end do
      return
      end

c======================================================================
c Last Modified: February 11, 1995

      subroutine daxpy(n,da,dx,dy)

c constant times a vector plus a vector.

      real*8 dx(*),dy(*),da
c     real*8 dx(1),dy(1),da
      integer*4 i,n

      if (.not.(dabs(da).gt.0.0d0)) return  ! jvk

      do i = 1,n
        dy(i) = dy(i) + da*dx(i)
      end do
      return
      end

c======================================================================
c Last Modified: February 11, 1995

      subroutine dscal(n,da,dx)

c scales a vector by a constant.

      real*8 da,dx(*)
c     real*8 da,dx(1)
      integer*4 i,n

      do i = 1,n
        dx(i) = da*dx(i)
      end do
      return
      end

c======================================================================
c Last Modified: February 1, 1995

      integer*4 function idamax(n,dx)

c finds the index of element having max. absolute value.

      real*8 dx(*),dmax
c     real*8 dx(1),dmax
      integer*4 i,n

      idamax = 1
      dmax = dabs(dx(1))           ! code for increment equal to 1
      do i = 2,n
        if (dabs(dx(i)).gt.dmax) then
          idamax = i
          dmax = dabs(dx(i))
        end if
      end do
      return
      end

c======================================================================
c Last Modified: February 1, 1995

      real*8 function dnrm2 ( n, dx, incx)
      integer*4          next
      real*8   dx(*), cutlo, cuthi, hitest, sum, xmax,zero,one
c     real*8   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /

      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300

   10 assign 30 to next
      sum = zero
      nn = n * incx

c begin main loop

      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero

c phase 1.  sum is zero

   50 if(.not.(dabs(dx(i)).gt.zero)) go to 200  ! jvk
      if( dabs(dx(i)) .gt. cutlo) go to 85

c prepare for phase 2.

      assign 70 to next
      go to 105

c prepare for phase 4.

  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115

c phase 2.  sum is small.
c scale to avoid destructive underflow.

   70 if( dabs(dx(i)) .gt. cutlo ) go to 75

c common code for phases 2 and 4.
c in phase 4 sum is large.  scale to avoid overflow.

  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200

  115 sum = sum + (dx(i)/xmax)**2
      go to 200

c prepare for phase 3.

   75 sum = (sum * xmax) * xmax

c for real or d.p. set hitest = cuthi/n
c for complex      set hitest = cuthi/(2*n)

   85 hitest = cuthi/dfloat( n )   !jvk

c phase 3.  sum is mid-range.  no scaling.

      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300

  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20

c end of main loop.
c compute square root and adjust for scaling.

      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end

c======================================================================
