        SUBROUTINE simplx(a,m,n,mp,np,m1,m2,m3,icase,izrov,iposv)
        INTEGER icase,m,m1,m2,m3,mp,n,np,iposv(m),izrov(n),MMAX,NMAX
        REAL a(mp,np),EPS
        PARAMETER (MMAX=100,NMAX=100,EPS=1.e-6)

C       USES simp1,simp2,simp3

        INTEGER i,ip,is,k,kh,kp,nl1,l1(NMAX),l3(MMAX)                  ! Simplex method for linear programming. Input parameters a, m, n, mp, np, m1, m2, and m3,
        REAL bmax,q1                                                   ! and output parameters a, icase, izrov, and iposv are described above. Parameters: MMAX is the
        if(m.ne.m1+m2+m3)pause 'bad input constraint counts in simplx' ! maximum number of constraints expected; NMAX is the maximum number of variables expected; 
        nl1=n                                                          ! EPS is the absolute precision, which should be adjusted to the scale of your variables.

        do 11 k=1,n                                                    
            l1(k)=k                                                    ! Initialize index list of columns admissible for exchange.
            izrov(k)=k                                                 ! Initially make all variables right-hand.
        enddo 11

        do 12 i=1,m                                                    ! Constants bi must be non-negative.
            if(a(i+1,1).lt.0.)pause 'bad input tableau in simplx'      ! Initial left-hand variables. m1 type constraints are represented by having their slack variable
                iposv(i)=n+i                                           ! initially left-hand, with no articial variable. m2 type constraints have their slack
        enddo 12                                                       ! variable initially left-hand, with a minus sign, and their articial variable handled implicitly
                                                                       ! during their first exchange. m3 type constraints have their articial variable initially left-hand.
    
        if(m2+m3.eq.0)goto 30                                          ! The origin is a feasible starting solution. Go to phase two.
        do 13 i=1,m2                                                   ! Initialize list of m2 constraints whose slack variables have never
            l3(i)=1                                                    ! been exchanged out of the initial basis.
        enddo 13

        do 15 k=1,n+1                                                  ! Compute the auxiliary objective function.
            q1=0.
            do 14 i=m1+1,m
                q1=q1+a(i+1,k)
            enddo 14
            a(m+2,k)=-q1
        enddo 15

10      call simp1(a,mp,np,m+1,l1,nl1,0,kp,bmax)                       ! Find max. coe. of auxiliary objective fn.

        if(bmax.le.EPS.and.a(m+2,1).lt.-EPS)then
            icase=-1
            return                                                     ! Auxiliary objective function is still negative and can't be improved, hence no feasible solution exists.
        else if(bmax.le.EPS.and.a(m+2,1).le.EPS)then
                                                                       ! Auxiliary objective function is zero and can't be improved; we have a feasible starting vector.
                                                                       ! Clean out the articial variables corresponding to any remaining equality constraints by
                                                                       ! goto 1's and then move on to phase two by goto 30.
            do 16 ip=m1+m2+1,m
                if(iposv(ip).eq.ip+n)then                              ! Found an articial variable for an equality constraint.
                call simp1(a,mp,np,ip,l1,nl1,1,kp,bmax) 
                    if(bmax.gt.EPS)goto 1
                endif                                                  ! Exchange with column corresponding to maximum pivot element in row.
            enddo 16
            do 18 i=m1+1,m1+m2                                         ! Change sign of row for any m2 constraints still present from the initial basis.
                if(l3(i-m1).eq.1)then                                  
                    do 17 k=1,n+1
                        a(i+1,k)=-a(i+1,k)
                    enddo 17
                endif
            enddo 18
        goto 30                                                         ! Go to phase two.
        endif

        call simp2(a,m,n,mp,np,ip,kp)                                   ! Locate a pivot element (phase one).
        if(ip.eq.0)                                                     ! then Maximum of auxiliary objective function is unbounded, so no feasible solution exists.
            icase=-1
            return
        endif
1       call simp3(a,mp,np,m+1,n,ip,kp)
        if(iposv(ip).ge.n+m1+m2+1)                                      ! Exchange a left- and a right-hand variable (phase one), then update lists.
        do 19 k=1,nl1                                                   ! then Exchanged out an articial variable for an equality constraint. 
            if(l1(k).eq.kp)goto 2                                       ! Make sure it stays out by removing it from the l1 list.
        enddo 19
2           nl1=nl1-1
        do 21 is=k,nl1
            l1(is)=l1(is+1)
        enddo 21
        else
            kh=iposv(ip)-m1-n
            if(kh.ge.1)then                                             ! Exchanged out an m2 type constraint.
                if(l3(kh).ne.0)then                                     !If it's the rst time, correct the pivot column for the minus sign and the implicit articial ariable.
                    l3(kh)=0
                    a(m+2,kp+1)=a(m+2,kp+1)+1.
                    do 22 i=1,m+2
                        a(i,kp+1)=-a(i,kp+1)
                    enddo 22
                endif
            endif
        endif
        is=izrov(kp)                                                    ! Update lists of left- and right-hand variables.
        izrov(kp)=iposv(ip)
        iposv(ip)=is
        goto 10                                                         ! Still in phase one, go back to 10.
                                                                        ! End of phase one code for nding an initial feasible solution. Now, in phase two, optimize it.
30      call simp1(a,mp,np,0,l1,nl1,0,kp,bmax)                          ! Test the z-row for doneness.
        if(bmax.le.EPS)then                                             ! Done. Solution found. Return with the good news.
            icase=0
            return
        endif

        call simp2(a,m,n,mp,np,ip,kp)                                   ! Locate a pivot element (phase two).
        if(ip.eq.0)then                                                 ! Objective function is unbounded. Report and return.
            icase=1
            return
        endif

        call simp3(a,mp,np,m,n,ip,kp)                                   ! Exchange a left- and a right-hand variable (phase two),
        is=izrov(kp)                                                    ! update lists of left- and right-hand variables,
        izrov(kp)=iposv(ip)
        iposv(ip)=is
        goto 30                                                         ! and return for another iteration.
        END
                                                                        ! The preceding routine makes use of the following utility subroutines.

        SUBROUTINE simp1(a,mp,np,mm,ll,nll,iabf,kp,bmax)
        INTEGER iabf,kp,mm,mp,nll,np,ll(np)
        REAL bmax,a(mp,np)
                                                                        ! Determines the maximum of those elements whose index is contained in the supplied list
                                                                        ! ll, either with or without taking the absolute value, as flagged by iabf.
        INTEGER k
        REAL test
        if(nll.le.0)then No eligible columns.
            bmax=0.
        else
            kp=ll(1)
            bmax=a(mm+1,kp+1)
            do 11 k=2,nll
                if(iabf.eq.0)then
                    test=a(mm+1,ll(k)+1)-bmax
                else
                    test=abs(a(mm+1,ll(k)+1))-abs(bmax)
                endif
                if(test.gt.0.)then
                    bmax=a(mm+1,ll(k)+1)
                    kp=ll(k)
                endif
            enddo 11
        endif
        return
        END

        SUBROUTINE simp2(a,m,n,mp,np,ip,kp)
        INTEGER ip,kp,m,mp,n,np
        REAL a(mp,np),EPS
        PARAMETER (EPS=1.e-6)
                                                                        ! Locate a pivot element, taking degeneracy into account.
        INTEGER i,k
        REAL q,q0,q1,qp
        ip=0
        do 11 i=1,m
            if(a(i+1,kp+1).lt.-EPS)goto 1
        enddo 11
        return                                                          ! No possible pivots. Return with message.
1       q1=-a(i+1,1)/a(i+1,kp+1)
        ip=i
        do 13 i=ip+1,m
            if(a(i+1,kp+1).lt.-EPS)then
                q=-a(i+1,1)/a(i+1,kp+1)
                if(q.lt.q1)then
                    ip=i
                    q1=q
                else if (q.eq.q1) then                                  ! We have a degeneracy.

                    do 12 k=1,n
                        qp=-a(ip+1,k+1)/a(ip+1,kp+1)
                        q0=-a(i+1,k+1)/a(i+1,kp+1)
                        if(q0.ne.qp)goto 2
                    enddo 12
2                   if(q0.lt.qp)ip=i
                    endif
                endif
        enddo 13
        return
        END

        SUBROUTINE simp3(a,mp,np,i1,k1,ip,kp)
        INTEGER i1,ip,k1,kp,mp,np
        REAL a(mp,np)
                                                                        ! Matrix operations to exchange a left-hand and right-hand variable (see text).
        INTEGER ii,kk
        REAL piv
        piv=1./a(ip+1,kp+1)
        do 12 ii=1,i1+1
            if(ii-1.ne.ip)then
                a(ii,kp+1)=a(ii,kp+1)*piv
                do 11 kk=1,k1+1
                    if(kk-1.ne.kp)then
                        a(ii,kk)=a(ii,kk)-a(ip+1,kk)*a(ii,kp+1)
                    endif
                enddo 11
            endif
        enddo 12
        do 13 kk=1,k1+1
            if(kk-1.ne.kp)a(ip+1,kk)=-a(ip+1,kk)*piv
        enddo 13
        a(ip+1,kp+1)=piv
        return
        END