PROGRAM TransLz
!Program to transform real, solid harmonic molecular one- and two-electron integrals into complex, pure harmonic molecular integrals.
!The program will read in the original integrals from an INTDUMP file and produce a PUREHARMINTDUMP file of complex integrals.

    IMPLICIT NONE
    INTEGER :: NORB,NELEC,MS2,ORBSYM(1000),ISYM,SYML(1000),SYMLZ(1000)
    INTEGER :: iPairs,iSize,ierr,I,J,K,L,UMatInd,PairMl,A,B,MaxLz,Gap
    INTEGER , ALLOCATABLE :: SymLCounts(:,:),SymLList(:),IND(:),LzPairs(:),RefOrbLz(:,:)
    REAL*8 , ALLOCATABLE :: UMAT(:),ARR(:),TMAT(:,:)
    REAL*8 :: Z,ECore
    REAL*8 :: Norm
    LOGICAL :: UHF,FlipInt,exists,tDiatomic,LinearExists,AtomExists,tSphericalSym,tPairbyEnergy,tError
    LOGICAL , ALLOCATABLE :: FlipSign(:)
    COMPLEX*16 :: CompZ,alpha1Coeff,alpha2Coeff,beta1Coeff,beta2Coeff,gamma1Coeff,gamma2Coeff,Delta1Coeff,Delta2Coeff
    INTEGER :: alpha1,alpha2,beta1,beta2,gamma1,gamma2,delta1,delta2,SymLTot,Degen
    LOGICAL :: tNoSym
    INTEGER :: iSyms
    NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM,UHF,SYML,SYMLZ

    UHF=.false.
    tError=.false.

    OPEN(8,FILE='INTDUMP',STATUS='OLD',FORM='FORMATTED')
    READ(8,FCI)

    IF(NORB.gt.1000) STOP 'Too Many Orbitals'
    IF(UHF) STOP 'Cannot cope with qchem UHF integrals (ROHF should be ok). &
                  &Bug ajwt if needed. Remove this stop only if you know what you are doing'
    

    Norm=1.D0/sqrt(2.D0)
    
    MaxLz=0
    do i=1,NORB
        IF(abs(SYMLz(i)).gt.MaxLz) MaxLz=abs(SYMLz(i))
    enddo

    INQUIRE(FILE='ATOMIC',EXIST=AtomExists)
    INQUIRE(FILE='LINEAR',EXIST=LinearExists)
    IF(AtomExists.and.LinearExists) THEN
        WRITE(6,*) "Both ATOMIC and LINEAR files found..."
        WRITE(6,*) "Is this system a linear molecule rather than atomic system? (T/F)"
        READ(*,*) tDiatomic
    ELSEIF(AtomExists) THEN
        tDiatomic=.false.
    ELSEIF(LinearExists) THEN
        tDiatomic=.true.
    ELSE
        WRITE(6,*) "Is this system a linear molecule (T) or an atomic system (F) ? "
        READ(*,*) tDiatomic
    ENDIF


    tPairbyEnergy=.false.
!Check that Lz orbitals are written out in the correct order
    IF((.not.UHF).and.(.not.tDiatomic)) THEN
        do i=1,NORB
            
            IF(abs(SYMLz(i)).gt.SYML(i)) THEN
                STOP 'Error in Ls'
            ENDIF
        enddo

        i=1
        do while(i.le.NORB)
            !New L Block
            Degen=2*abs(SYMLz(i))+1

            do j=i+1,i+Degen-1
                IF(SYMLz(j).ne.SYMLz(j-1)+1) THEN
                    STOP 'Error with Lz ordering'
                ENDIF
            enddo
            i=i+Degen
        enddo


    ELSEIF(UHF) THEN
        tSphericalSym=.false.   !This is true if spherical symmetry is maintained (ie all states are pure L states)
        lp: do i=1,NORB
            IF(SYMLz(i).eq.0) THEN
                CYCLE
            ELSEIF(SYMLz(i).lt.0) THEN
                IF(SYMLz(i).ne.-SYMLz(i+2)) THEN
                    tSphericalSym=.true.
                    Do j=1,NORB
                        IF(SYML(j).eq.-20) THEN
                            WRITE(6,*) "******   WARNING   ******"
                            WRITE(6,*) "ROHF System, with spherical symmetry broken and unknown symmetry ordering"
                            WRITE(6,*) "Ordering system instead simply by energy. "
                            CALL FLUSH(6)
                            tPairbyEnergy=.true.
                            EXIT lp
                        ENDIF
                    ENDDO
                ENDIF

            ELSEIF(SYMLz(i).gt.0) THEN
                IF(SYMLz(i).ne.-SYMLz(i-2)) THEN
                    tSphericalSym=.true.
                    Do j=1,NORB
                        IF(SYML(j).eq.-20) THEN
                            WRITE(6,*) "******   WARNING   ******"
                            WRITE(6,*) "ROHF System, with spherical symmetry broken and unknown symmetry ordering"
                            WRITE(6,*) "Ordering system instead simply by energy. "
                            CALL FLUSH(6)
                            tPairbyEnergy=.true.
                            EXIT lp
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        enddo lp

!Check that the Lz are consistent between an alpha/beta pair
        do i=1,NORB-1,2
            IF(SYMLz(i).ne.SYMLz(i+1)) STOP 'Error with alpha/beta ordering'
        enddo

!        do while(i.le.NORB)
!            Degen=4*abs(SYMLz(i))+2
!
!            do j=i+2,i+Degen-4
!                IF(SYMLz(j).ne.SYMLz(j-2)+1) THEN
!                    STOP 'Error with Lz ordering'
!                ENDIF
!            enddo
!            do j=i+3,i+Degen-3
!                IF(SYMLz(j).ne.SYMLz(j-2)+1) THEN
!                    STOP 'Error with Lz ordering'
!                ENDIF
!            enddo
!            i=i+Degen
!        enddo
!
    ELSE
        do i=1,NORB
            IF(SYMLz(i).eq.0) THEN
                CYCLE
            ELSEIF(SYMLz(i).lt.0) THEN
                IF(SYMLz(i).ne.-SYMLz(i+1)) THEN 
                    WRITE(6,*) 'Error with Lz ordering'
                    WRITE(6,*) "Attempting to order orbitals energetically"
                    tPairbyEnergy=.true.
                    EXIT
                ENDIF
            ELSEIF(SYMLz(i).gt.0) THEN
                IF(SYMLz(i).ne.-SYMLz(i-1)) THEN
                    WRITE(6,*) 'Error with Lz ordering'
                    WRITE(6,*) "Attempting to order orbitals energetically"
                    tPairbyEnergy=.true.
                    EXIT
                ENDIF
            ENDIF
        enddo

    ENDIF


!The orbitals are written out one L state at a time, however, we cannot be sure that the Ml states are in the correct order.
!Therefore, we need to create a indexing scheme, so that each orbital can find its prospective partner.
    ALLOCATE(LzPairs(NORB),stat=ierr)
    LzPairs(:)=0

    IF(tPairByEnergy) THEN
!We didn't manage to find a consistent scheme to pair off the Lz values. Instead, pair them off by energy. Do in order of increasing energy, and decreasing energy, and hope they agree...
!Nicely, the orbitals are already ordered by fock energy, and so we can 
        IF(UHF) THEN
            do i=NORB-1,1,-2
!                WRITE(6,*) i

                IF(LzPairs(i).ne.0) THEN
                    CYCLE
                ELSEIF(LzPairs(i).gt.i) THEN
                    !We have already paired this orbital up with something higher than it.
                    CYCLE
                ELSEIF(SYMLz(i).eq.0) THEN
                    LzPairs(i)=i
                    LzPairs(i+1)=i+1
                ELSE
                    do j=i-2,1,-2
                        IF(SYMLz(j).eq.-SYMLz(i)) THEN
                            !Found pair
                            IF(LzPairs(i).ne.0) STOP 'Pair wanted already taken'
                            IF(LzPairs(j).ne.0) STOP 'Pair wanted already taken'
                            LzPairs(i)=j
                            LzPairs(j)=i
                            LzPairs(i+1)=j+1
                            LzPairs(j+1)=i+1
                            EXIT
                        ELSEIF(j.eq.1) THEN
                            STOP 'Cannot pair all orbitals'
                        ENDIF
                    enddo
                ENDIF
            enddo

        ELSE
!The orbitals *should* be ordered as -m,m, however due to the fact that the required offset to split degeneracies between orbitals can cause angular functions to overlap, we'll have to assume that they will
!pair by energy.
!There is a slight problem with this. at dissociation, atomic functions on the different atoms are degenerate, and so we will get things like -1,-1,0,0,1,1 for the atomic p functions.
!Therefore, we will just have to take the next pair if one is already taken.
!This will be the same ordering as atomic systems
            do i=1,NORB
                IF(SYMLz(i).eq.0) THEN
                    LzPairs(i)=i
                ELSEIF(LzPairs(i).gt.0) THEN !If we've already been allocated, continue
                    CYCLE
                ELSE
                    !cycle up through the orbitals until the positive value which corresponds is found.
                    do j=i+1,NORB
                        IF((SYMLz(j).eq.-SYMLz(i)).and.(LzPairs(j).eq.0)) THEN
!                            IF(LzPairs(i).ne.0) STOP 'Pair wanted already taken'
!                            IF(LzPairs(j).ne.0) STOP 'Pair wanted already taken'
                            LzPairs(i)=j
                            LzPairs(j)=i
!                            write(6,*) "Pair", i,j
                            EXIT
                        ELSEIF(j.eq.NORB) THEN
                            write(6,*) "Cannot find pair for orbital ", i
                            STOP 'Cannot pair all orbitals'
                        ENDIF
                    enddo
                ENDIF
            enddo
            do i=1,NORB
                IF(LzPairs(i).eq.0) STOP 'Orbitals not paired correctly'
            enddo

        ENDIF

!Now we want to check that this ordering is the same as if we wanted to order in increasing energy
        IF(UHF) THEN
            do i=1,NORB-1,2
!                WRITE(6,*) i
                IF(LzPairs(i).eq.0) THEN
                    STOP 'Not all orbitals paired'
                ELSEIF(LzPairs(i).lt.i) THEN
                    !We have already paired this orbital up with something lower than it.
                    CYCLE
                ELSEIF(SYMLz(i).eq.0) THEN
                    IF(LzPairs(i).ne.i) STOP '0 ml orbital not paired with itself.'
                ELSE
                    do j=i+2,NORB-1,2
                        IF(SYMLz(j).eq.-SYMLz(i)) THEN
                            !Found pair
                            IF(LzPairs(i).ne.j) STOP  &
                                   'Pairing by increasing energy does not give the same answer as decreasing energy'
                            IF(LzPairs(j).ne.i) STOP  & 
                                   'Pairing by increasing energy does not give the same answer as decreasing energy'
                            EXIT
                        ELSEIF(j.eq.1) THEN
                            STOP 'Cannot pair all orbitals in increasing energy'
                        ENDIF
                    enddo
                ENDIF
            enddo
        ELSE
!            do i=1,NORB
!                IF(LzPairs(i).eq.0) THEN
!                    STOP 'Not all orbitals paired'
!                ELSEIF(LzPairs(i).lt.i) THEN
!                    !We have already paired this orbital up with something lower than it.
!                    CYCLE
!                ELSEIF(SYMLz(i).eq.0) THEN
!                    IF(LzPairs(i).ne.i) STOP '0 ml orbital not paired with itself.'
!                ELSE
!                    do j=i+1,NORB
!                        IF(SYMLz(j).eq.-SYMLz(i)) THEN
!                            !Found pair
!                            IF(LzPairs(i).ne.j) STOP  &
!                                   'Pairing by increasing energy does not give the same answer as decreasing energy'
!                            IF(LzPairs(j).ne.i) STOP  &
!                                   'Pairing by increasing energy does not give the same answer as decreasing energy'
!                            EXIT
!                        ELSEIF(j.eq.1) THEN
!                            STOP 'Cannot pair all orbitals in increasing energy'
!                        ENDIF
!                    enddo
!                ENDIF
!            enddo
        ENDIF

    ELSE

!LzPairs(i) is the index of the orbital from the same angular momentum function, but with opposite Lz value.
        IF((.not.UHF).and.(.not.tDiatomic)) THEN
!LzPairs(i) is the index of the orbital from the same angular momentum function, but with opposite Lz value.
            i=1
            do while(i.le.NORB)
                IF(SymL(i).eq.0) THEN
                    LzPairs(i)=i
                    i=i+1
                ELSE
                    do j=i,2*SymL(i)+i
                        !Just running over block of a given L value
                        PairMl=-SymLz(j)
                        do k=i,2*SymL(i)+i
                            IF(SymLz(k).eq.PairMl) THEN
                                LzPairs(j)=k
                                EXIT
                            ENDIF
                        enddo
                    enddo
                    i=i+2*SymL(i)+1
                ENDIF
            enddo
        ELSEIF(UHF) THEN
            IF(tSphericalSym) THEN
!Ordered -2,-2,-1,-1,0,0,1,1,2,2,
!We are now assuming that Lz will always increase in an L state, Lz of the first orbital is the 'L' value.
                do i=1,NORB
                    IF(SYMLz(i).eq.0) THEN
                        LzPairs(i)=i
                    ELSE
                        Gap=2*(2*SYML(i)+1)-2-2*(2*SYML(i)-2*abs(SYMLz(i)))
                        IF(SYMLz(i).lt.0) THEN
                            LzPairs(i)=i+Gap
                        ELSE
                            LzPairs(i)=i-Gap
                        ENDIF
                    ENDIF
                enddo
            ELSE
                do i=1,NORB
                    IF(SYMLz(i).eq.0) THEN
                        LzPairs(i)=i
                    ELSEIF(SYMLz(i).lt.0) THEN
                        LzPairs(i)=i+2
                    ELSE
                        LzPairs(i)=i-2
                    ENDIF
                enddo
            ENDIF
        ELSE
            do i=1,NORB
                IF(SYMLz(i).eq.0) THEN
                    LzPairs(i)=i
                ELSEIF(SYMLz(i).lt.0) THEN
                    LzPairs(i)=i+1
                ELSE
                    LzPairs(i)=i-1
                ENDIF
            enddo
        ENDIF
    ENDIF

    do i=1,NORB
        WRITE(6,*) i,LzPairs(i)
    enddo

    iPairs=(NORB*(NORB+1))/2
    iSize=(iPairs*(iPairs+1))/2

    ALLOCATE(UMAT(0:iSize),stat=ierr)
    IF(ierr.ne.0) STOP 'Allocate error'
    ALLOCATE(ARR(NORB),stat=ierr)
    ALLOCATE(TMAT(0:NORB,0:NORB),stat=ierr)
    IF(ierr.ne.0) STOP 'Allocate error'
    UMAT(:)=0.D0
    TMAT(:,:)=0.D0
    ARR(:)=0.D0

!I now need to create a consistent x-direction. 
!This could be done by making sure that the one electron hamiltonian integral between orbitals of the same Ml are positive, though there will be symmetry problems.
!This is now done internally in qchem...
    
    DO WHILE(.true.)

        READ(8,'(1X,G20.14,4I3)',END=99) Z,I,J,K,L    
!        WRITE(6,*) I

        IF(I.ne.0.and.J.eq.0.and.K.eq.0) THEN
            !Fock energies
            Arr(I)=Z
!            WRITE(19,'(1X,G20.12,4I3)') Z,i,j,k,l
        ELSEIF(I.ne.0.and.J.ne.0.and.K.eq.0.and.L.eq.0) THEN
            !One-electron energies

            TMAT(I,J)=Z
            TMAT(J,I)=Z
!            WRITE(19,'(1X,G20.12,4I3)') Z,i,j,k,l
        ELSEIF(I.eq.0) THEN
            !Core energy
            ECore=Z
!            WRITE(19,'(1X,G20.12,4I3)') Z,i,j,k,l
        ELSE
            !2e int

            UMAT(UMatInd(I,J,K,L))=Z
!            WRITE(19,'(1X,G20.12,4I3)') Z,i,j,k,l
        ENDIF

    END DO

99  CONTINUE

!Integrals have now been read in - now to transform them to complex integrals.
    CLOSE(8)

!First, calculate whether we are a heteronuclear diatomic or homonuclear
!If heteronuclear - will have 4 (1 -> 4 in qchem) irreps - remove all symmetry info
!If homo - will have 8 (1 -> 8 in qchem) irreps - then keep the inversion symmetry - reduce to Ci point group.
!    TestSym(:)=.false.
!    do i=1,NORB
!        TestSym(ORBSYM(i))=.true.
!    enddo
!    iSyms=0
!    do i=1,8
!        if(TestSym(i)) then
!            iSyms=iSyms+1
!        endif
!    enddo
!    if((iSyms.ne.4).and.(iSyms.ne.8)) then
!        STOP 'Neither 4 or 8 irreps - neither D2h or C2v point group detected. Error.'
!    endif
    iSyms=0
    tNoSym=.false.
    do i=1,NORB
        if(ORBSYM(i).gt.iSyms) then
            iSyms=ORBSYM(i)
        elseif(ORBSYM(i).eq.0) then
            tNoSym=.true.
        endif
    enddo

    if((iSyms.le.4).or.tNoSym) then
        !Heteronuclear - Remove all sym - all calculated in Lz and now redundant, and can't keep +/- sym (unless kept in dets)
        if(tNoSym) then
            write(6,'(a)') "Symmetry irreps missing - removing all symmetry"
        else
            write(6,'(a)') "Heteronuclear molecule detected - removing all redundant point group symmetry"
        endif
        ORBSYM(:)=0
    elseif(iSyms.gt.4) then
        !Keep inversion symmetry. **** In qchem **** (i.e. not molpro), this irreps 1 -> 4 (g) and 5 -> 8 (u)
        write(6,'(a)') "Homonuclear molecule detected - reducing point group of molecule to Ci (keeping only inversion)"
        do i=1,NORB
            if(ORBSYM(i).ge.5) then
                ORBSYM(i) = 2
            else
                ORBSYM(i) = 1
            endif
        enddo
    endif



    OPEN(8,FILE='PUREHARMINTDUMP',STATUS='UNKNOWN',FORM='FORMATTED')

    WRITE(8,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',NORB,'NELEC=',NELEC,',MS2=',MS2,','
    WRITE(8,'(A9)',advance='no') 'ORBSYM='
    DO i=1,NORB
!Explicitly turn off symmetry
!        WRITE(8,'(I1,A1)',advance='no') 0,','
!Reduce symmetry from D2h to C2h
        WRITE(8,'(I1,A1)',advance='no') ORBSYM(i),','
    ENDDO
    WRITE(8,*) 
    IF(UHF) THEN
        WRITE(8,'(A7,I1,A12)') 'ISYM=',1,' UHF=.TRUE.'
    ELSE
        WRITE(8,'(A7,I1,A12)') 'ISYM=',1,' UHF=.FALSE.'
    ENDIF
    WRITE(8,'(A7)',advance='no') 'SYML='
    DO i=1,NORB
        WRITE(8,'(I3,A1)',advance='no') SYML(i),','
    ENDDO
    WRITE(8,*) 
    WRITE(8,'(A8)',advance='no') 'SYMLZ='
    DO i=1,NORB
        WRITE(8,'(I2,A1)',advance='no') SYMLZ(i),','
    ENDDO
    WRITE(8,*) 
    WRITE(8,'(A5)') '&END'

    DO i=1,NORB
        DO j=1,NORB
            A=(i*(i-1))/2+j
            DO k=1,NORB
                DO l=1,NORB
                    B=(k*(k-1))/2+l

                    IF(B.lt.A) CYCLE
                    IF((i.lt.j).and.(k.lt.l)) CYCLE
                    IF((i.gt.j).and.(k.lt.l)) CYCLE

                    !Find alpha orbitals which will combine to make i
                    !These need to be from the same L function, but have opposite Lz value
                    IF(SymLz(i).eq.0) THEN
                    !If Lz(i)=0, then alpha = i. This will include all S states and Ml=0 states of higher angular momentum functions.

                        Alpha1=i
                        Alpha2=0
                        Alpha1Coeff=DCMPLX(1.D0,0.D0)
                        Alpha2Coeff=DCMPLX(0.D0,0.D0)

                    ELSE
                    !If Lz(i).ne.0, then there will be two contributary solid harmonic
                    !Alpha1 will always be the orbital with positive ml value, and Alpha2 will have negative ml.

!                        IF(SymL(i).ne.1) THEN
                        IF(SymLz(i).lt.0) THEN  !Molecular index is negative
                            !Negative Ml - this will be alpha2
                            alpha2=i
                            alpha1=LzPairs(i)
!                            alpha1=i+2*abs(SymLz(i))    !This will be the corresponding positive Ml state
                            !Now find coefficients. Alpha needs the coefficients to be the complex conjugates...
                            IF(mod(abs(SymLz(i)),2).eq.1) THEN
                                !We have an ODD ml pair of orbitals
                                Alpha1Coeff=DCMPLX(Norm,0.D0)
                                Alpha2Coeff=DCMPLX(0.D0,Norm)
                            ELSE
                                !We have an EVEN ml pair of orbitals
                                Alpha1Coeff=DCMPLX(Norm,0.D0)
                                Alpha2Coeff=DCMPLX(0.D0,Norm)
                            ENDIF
                        ELSE
                            !Positive Ml
                            alpha1=i
                            alpha2=LzPairs(i)
!                            alpha2=i-2*abs(SymLz(i))
                            IF(mod(abs(SymLz(i)),2).eq.1) THEN
                                !We have an ODD ml pair of orbitals
                                Alpha1Coeff=DCMPLX(-Norm,0.D0)
                                Alpha2Coeff=DCMPLX(0.D0,Norm)
                            ELSE
                                !We have an EVEN ml pair of orbitals
                                Alpha1Coeff=DCMPLX(Norm,0.D0)
                                Alpha2Coeff=DCMPLX(0.D0,-Norm)
                            ENDIF
                        ENDIF
                        IF(SymLz(alpha1).ne.-SymLz(alpha2)) THEN
                            STOP 'Alphas not chosen correctly'
                        ENDIF
                    ENDIF

                    !Now we need to do the same for the betas...
                    !Find beta orbitals which will combine to make j
                    !These need to be from the same L function, but have opposite Lz value
                    IF(SymLz(j).eq.0) THEN
                    !If Lz(j)=0, then beta = j. This will include all S states and Ml=0 states of higher angular momentum functions.

                        beta1=j
                        beta2=0
                        beta1Coeff=DCMPLX(1.D0,0.D0)
                        beta2Coeff=DCMPLX(0.D0,0.D0)

                    ELSE
                    !If Lz(j).ne.0, then there will be two contributary solid harmonic
                    !beta1 will always be the orbital with positive ml value, and beta2 will have negative ml.

!                        IF(SymL(j).ne.1) THEN
                        IF(SymLz(j).lt.0) THEN
                            !Negative Ml - this will be beta2
                            beta2=j
                            beta1=LzPairs(j)
!                            beta1=j+2*abs(SymLz(j))    !This will be the corresponding positive Ml state
                            !Now find coefficients. 
                            IF(mod(abs(SymLz(j)),2).eq.1) THEN
                                !We have an ODD ml pair of orbitals
                                Beta1Coeff=DCMPLX(Norm,0.D0)
                                Beta2Coeff=DCMPLX(0.D0,-Norm)
                            ELSE
                                !We have an EVEN ml pair of orbitals
                                Beta1Coeff=DCMPLX(Norm,0.D0)
                                Beta2Coeff=DCMPLX(0.D0,-Norm)
                            ENDIF
                        ELSE
                            beta1=j
                            beta2=LzPairs(j)
!                            beta2=j-2*abs(SymLz(j))
                            IF(mod(abs(SymLz(j)),2).eq.1) THEN
                                !We have an ODD ml pair of orbitals
                                Beta1Coeff=DCMPLX(-Norm,0.D0)
                                Beta2Coeff=DCMPLX(0.D0,-Norm)
                            ELSE
                                !We have an EVEN ml pair of orbitals
                                Beta1Coeff=DCMPLX(Norm,0.D0)
                                Beta2Coeff=DCMPLX(0.D0,Norm)
                            ENDIF
                        ENDIF
                        IF(SymLz(beta1).ne.-SymLz(beta2)) THEN
                            STOP 'betas not chosen correctly'
                        ENDIF
                    ENDIF
                    
                    
                    !and the gammas...
                    !Find gamma orbitals which will combine to make k
                    !These need to be from the same L function, but have opposite Lz value
                    IF(SymLz(k).eq.0) THEN
                    !If Lz(k)=0, then gamma = k. This will include all S states and Ml=0 states of higher angular momentum functions.

                        gamma1=k
                        gamma2=0
                        gamma1Coeff=DCMPLX(1.D0,0.D0)
                        gamma2Coeff=DCMPLX(0.D0,0.D0)

                    ELSE
                    !If Lz(k).ne.0, then there will be two contributary solid harmonic
                    !gamma1 will always be the orbital with positive ml value, and gamma2 will have negative ml.

!                        IF(SymL(k).ne.1) THEN
                        IF(SymLz(k).lt.0) THEN
                            !Negative Ml - this will be gamma2
                            gamma2=k
                            gamma1=LzPairs(k)
!                            gamma1=k+2*abs(SymLz(k))    !This will be the corresponding positive Ml state
                            !Now find coefficients. Gammas will be complex conjugates
                            IF(mod(abs(SymLz(k)),2).eq.1) THEN
                                !We have an ODD ml pair of orbitals
                                Gamma1Coeff=DCMPLX(Norm,0.D0)
                                Gamma2Coeff=DCMPLX(0.D0,Norm)
                            ELSE
                                !We have an EVEN ml pair of orbitals
                                Gamma1Coeff=DCMPLX(Norm,0.D0)
                                Gamma2Coeff=DCMPLX(0.D0,Norm)
                            ENDIF
                        ELSE
                            gamma1=k
                            gamma2=LzPairs(k)
!                            gamma2=k-2*abs(SymLz(k))
                            IF(mod(abs(SymLz(k)),2).eq.1) THEN
                                !We have an ODD ml pair of orbitals
                                Gamma1Coeff=DCMPLX(-Norm,0.D0)
                                Gamma2Coeff=DCMPLX(0.D0,Norm)
                            ELSE
                                !We have an EVEN ml pair of orbitals
                                Gamma1Coeff=DCMPLX(Norm,0.D0)
                                Gamma2Coeff=DCMPLX(0.D0,-Norm)
                            ENDIF
                        ENDIF
                        IF(SymLz(gamma1).ne.-SymLz(gamma2)) THEN
                            STOP 'gammas not chosen correctly'
                        ENDIF
                    ENDIF


                    !and finally the deltas!...
                    !Find delta orbitals which will combine to make l
                    !These need to be from the same L function, but have opposite Lz value
                    IF(SymLz(l).eq.0) THEN
                    !If Lz(l)=0, then delta = l. This will include all S states and Ml=0 states of higher angular momentum functions.

                        delta1=l
                        delta2=0
                        delta1Coeff=DCMPLX(1.D0,0.D0)
                        delta2Coeff=DCMPLX(0.D0,0.D0)

                    ELSE
                    !If Lz(l).ne.0, then there will be two contributary solid harmonic
                    !delta1 will always be the orbital with positive ml value, and delta2 will have negative ml.

!                        IF(SymL(l).ne.1) THEN
                        IF(SymLz(l).lt.0) THEN
                            !Negative Ml - this will be delta2
                            delta2=l
                            delta1=LzPairs(l)
!                            delta1=l+2*abs(SymLz(l))    !This will be the corresponding positive Ml state
                            !Now find coefficients...
                            IF(mod(abs(SymLz(l)),2).eq.1) THEN
                                !We have an ODD ml pair of orbitals
                                Delta1Coeff=DCMPLX(Norm,0.D0)
                                Delta2Coeff=DCMPLX(0.D0,-Norm)
                            ELSE
                                !We have an EVEN ml pair of orbitals
                                Delta1Coeff=DCMPLX(Norm,0.D0)
                                Delta2Coeff=DCMPLX(0.D0,-Norm)
                            ENDIF
                        ELSE
                            delta1=l
                            delta2=LzPairs(l)
!                            delta2=l-2*abs(SymLz(l))
                            IF(mod(abs(SymLz(l)),2).eq.1) THEN
                                !We have an ODD ml pair of orbitals
                                Delta1Coeff=DCMPLX(-Norm,0.D0)
                                Delta2Coeff=DCMPLX(0.D0,-Norm)
                            ELSE
                                !We have an EVEN ml pair of orbitals
                                Delta1Coeff=DCMPLX(Norm,0.D0)
                                Delta2Coeff=DCMPLX(0.D0,Norm)
                            ENDIF
                        ENDIF
                        IF(SymLz(delta1).ne.-SymLz(delta2)) THEN
                            STOP 'deltas not chosen correctly'
                        ENDIF
                    ENDIF

!We now need to sum the contributions correctly!
!                    IF(i.eq.5.and.j.eq.1.and.k.eq.2.and.l.eq.5) THEN
!                        WRITE(6,*) Alpha1Coeff,Beta1Coeff,Gamma1Coeff,Delta1Coeff,UMAT(UMatInd(alpha1,beta1,gamma1,delta1)),alpha1,beta1,gamma1,delta1
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha2Coeff,Beta1Coeff,Gamma1Coeff,Delta1Coeff,UMAT(UMatInd(alpha2,beta1,gamma1,delta1)),alpha2,beta1,gamma1,delta1
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha1Coeff,Beta2Coeff,Gamma1Coeff,Delta1Coeff,UMAT(UMatInd(alpha1,beta2,gamma1,delta1)),alpha1,beta2,gamma1,delta1
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha2Coeff,Beta2Coeff,Gamma1Coeff,Delta1Coeff,UMAT(UMatInd(alpha2,beta2,gamma1,delta1)),alpha2,beta2,gamma1,delta1
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha1Coeff,Beta1Coeff,Gamma2Coeff,Delta1Coeff,UMAT(UMatInd(alpha1,beta1,gamma2,delta1)),alpha1,beta1,gamma2,delta1
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha2Coeff,Beta1Coeff,Gamma2Coeff,Delta1Coeff,UMAT(UMatInd(alpha2,beta1,gamma2,delta1)),alpha2,beta1,gamma2,delta1
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha1Coeff,Beta2Coeff,Gamma2Coeff,Delta1Coeff,UMAT(UMatInd(alpha1,beta2,gamma2,delta1)),alpha1,beta2,gamma2,delta1
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha2Coeff,Beta2Coeff,Gamma2Coeff,Delta1Coeff,UMAT(UMatInd(alpha2,beta2,gamma2,delta1)),alpha2,beta2,gamma2,delta1
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha1Coeff,Beta1Coeff,Gamma1Coeff,Delta2Coeff,UMAT(UMatInd(alpha1,beta1,gamma1,delta2)),alpha1,beta1,gamma1,delta2
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha2Coeff,Beta1Coeff,Gamma1Coeff,Delta2Coeff,UMAT(UMatInd(alpha2,beta1,gamma1,delta2)),alpha2,beta1,gamma1,delta2
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha1Coeff,Beta2Coeff,Gamma1Coeff,Delta2Coeff,UMAT(UMatInd(alpha1,beta2,gamma1,delta2)),alpha1,beta2,gamma1,delta2
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha2Coeff,Beta2Coeff,Gamma1Coeff,Delta2Coeff,UMAT(UMatInd(alpha2,beta2,gamma1,delta2)),alpha2,beta2,gamma1,delta2
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha1Coeff,Beta1Coeff,Gamma2Coeff,Delta2Coeff,UMAT(UMatInd(alpha1,beta1,gamma2,delta2)),alpha1,beta1,gamma2,delta2
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha2Coeff,Beta1Coeff,Gamma2Coeff,Delta2Coeff,UMAT(UMatInd(alpha2,beta1,gamma2,delta2)),alpha2,beta1,gamma2,delta2
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha1Coeff,Beta2Coeff,Gamma2Coeff,Delta2Coeff,UMAT(UMatInd(alpha1,beta2,gamma2,delta2)),alpha1,beta2,gamma2,delta2
!                        WRITE(6,*) "***"
!                        WRITE(6,*) Alpha2Coeff,Beta2Coeff,Gamma2Coeff,Delta2Coeff,UMAT(UMatInd(alpha2,beta2,gamma2,delta2)),alpha2,beta2,gamma2,delta2
!                        WRITE(6,*) "***"
!                    ENDIF
                    CompZ=Alpha1Coeff*Beta1Coeff*Gamma1Coeff*Delta1Coeff*UMAT(UMatInd(alpha1,beta1,gamma1,delta1))
                    CompZ=CompZ+Alpha2Coeff*Beta1Coeff*Gamma1Coeff*Delta1Coeff*UMAT(UMatInd(alpha2,beta1,gamma1,delta1))
                    CompZ=CompZ+Alpha1Coeff*Beta2Coeff*Gamma1Coeff*Delta1Coeff*UMAT(UMatInd(alpha1,beta2,gamma1,delta1))
                    CompZ=CompZ+Alpha2Coeff*Beta2Coeff*Gamma1Coeff*Delta1Coeff*UMAT(UMatInd(alpha2,beta2,gamma1,delta1))
                    CompZ=CompZ+Alpha1Coeff*Beta1Coeff*Gamma2Coeff*Delta1Coeff*UMAT(UMatInd(alpha1,beta1,gamma2,delta1))
                    CompZ=CompZ+Alpha2Coeff*Beta1Coeff*Gamma2Coeff*Delta1Coeff*UMAT(UMatInd(alpha2,beta1,gamma2,delta1))
                    CompZ=CompZ+Alpha1Coeff*Beta2Coeff*Gamma2Coeff*Delta1Coeff*UMAT(UMatInd(alpha1,beta2,gamma2,delta1))
                    CompZ=CompZ+Alpha2Coeff*Beta2Coeff*Gamma2Coeff*Delta1Coeff*UMAT(UMatInd(alpha2,beta2,gamma2,delta1))
                    CompZ=CompZ+Alpha1Coeff*Beta1Coeff*Gamma1Coeff*Delta2Coeff*UMAT(UMatInd(alpha1,beta1,gamma1,delta2))
                    CompZ=CompZ+Alpha2Coeff*Beta1Coeff*Gamma1Coeff*Delta2Coeff*UMAT(UMatInd(alpha2,beta1,gamma1,delta2))
                    CompZ=CompZ+Alpha1Coeff*Beta2Coeff*Gamma1Coeff*Delta2Coeff*UMAT(UMatInd(alpha1,beta2,gamma1,delta2))
                    CompZ=CompZ+Alpha2Coeff*Beta2Coeff*Gamma1Coeff*Delta2Coeff*UMAT(UMatInd(alpha2,beta2,gamma1,delta2))
                    CompZ=CompZ+Alpha1Coeff*Beta1Coeff*Gamma2Coeff*Delta2Coeff*UMAT(UMatInd(alpha1,beta1,gamma2,delta2))
                    CompZ=CompZ+Alpha2Coeff*Beta1Coeff*Gamma2Coeff*Delta2Coeff*UMAT(UMatInd(alpha2,beta1,gamma2,delta2))
                    CompZ=CompZ+Alpha1Coeff*Beta2Coeff*Gamma2Coeff*Delta2Coeff*UMAT(UMatInd(alpha1,beta2,gamma2,delta2))
                    CompZ=CompZ+Alpha2Coeff*Beta2Coeff*Gamma2Coeff*Delta2Coeff*UMAT(UMatInd(alpha2,beta2,gamma2,delta2))

!                    IF(ABS(CompZ).gt.1.D-09) THEN
!                        WRITE(8,'(1X,2G20.12,4I3)') CompZ,i,j,k,l
!                    ENDIF
                     IF(((SymLz(i)+SymLz(k)).ne.(SymLz(j)+SymLz(l))).and.ABS(CompZ).gt.1.D-8) THEN
                         IF(ABS(CompZ).lt.1.D-5) THEN
                             WRITE(6,*) 'Conservation of Ml error - setting to zero'
                             WRITE(6,*) i,j,k,l,CompZ
                             CompZ=DCMPLX(0.D0,0.D0)
                         ELSE
                             WRITE(6,*) i,j,k,l,CompZ
                             WRITE(6,*) 'Conservation of Ml error - stopping'
                             tError=.true.
                         ENDIF
                     ENDIF
                     IF(ABS(AIMAG(CompZ)).gt.1.D-9) THEN
                         STOP 'Complex integral found'
                     ENDIF
                     IF(ABS(DBLE(CompZ)).gt.1.D-9) THEN
                         WRITE(8,*) DBLE(CompZ),i,j,k,l
!                         WRITE(8,'(1X,G20.14,4I3)') DBLE(CompZ),i,j,k,l
                     ENDIF

                ENDDO
            ENDDO
        ENDDO
    ENDDO

!Now to do one-electron integrals.
    DO i=1,NORB
        DO j=i,NORB
            !Find alpha orbitals which will combine to make i
            !These need to be from the same L function, but have opposite Lz value
            IF(SymLz(i).eq.0) THEN
            !If Lz(i)=0, then alpha = i. This will include all S states and Ml=0 states of higher angular momentum functions.

                Alpha1=i
                Alpha2=0
                Alpha1Coeff=DCMPLX(1.D0,0.D0)
                Alpha2Coeff=DCMPLX(0.D0,0.D0)

            ELSE
            !If Lz(i).ne.0, then there will be two contributary solid harmonic
            !Alpha1 will always be the orbital with positive ml value, and Alpha2 will have negative ml.

!                        IF(SymL(i).ne.1) THEN
                IF(SymLz(i).lt.0) THEN  !Molecular index is negative
                    !Negative Ml - this will be alpha2
                    alpha2=i
                    alpha1=LzPairs(i)
!                            alpha1=i+2*abs(SymLz(i))    !This will be the corresponding positive Ml state
                    !Now find coefficients. Alpha needs the coefficients to be the complex conjugates...
                    IF(mod(abs(SymLz(i)),2).eq.1) THEN
                        !We have an ODD ml pair of orbitals
                        !Alpha1Coeff=CMPLX(-Norm,0.D0)
                        !Alpha2Coeff=CMPLX(0.D0,-Norm)
                        Alpha1Coeff=DCMPLX(Norm,0.D0)
                        Alpha2Coeff=DCMPLX(0.D0,Norm)
                    ELSE
                        !We have an EVEN ml pair of orbitals
                        !Alpha1Coeff=CMPLX(0.D0,-Norm)
                        !Alpha2Coeff=CMPLX(Norm,0.D0)
                        Alpha1Coeff=DCMPLX(Norm,0.D0)
                        Alpha2Coeff=DCMPLX(0.D0,Norm)
                    ENDIF
                ELSE
                    !Positive Ml
                    alpha1=i
                    alpha2=LzPairs(i)
!                            alpha2=i-2*abs(SymLz(i))
                    IF(mod(abs(SymLz(i)),2).eq.1) THEN
                        !We have an ODD ml pair of orbitals
                        Alpha1Coeff=DCMPLX(-Norm,0.D0)
                        Alpha2Coeff=DCMPLX(0.D0,Norm)
                    ELSE
                        !We have an EVEN ml pair of orbitals
                        Alpha1Coeff=DCMPLX(Norm,0.D0)
                        Alpha2Coeff=DCMPLX(0.D0,-Norm)
                    ENDIF
                ENDIF
                IF(SymLz(alpha1).ne.-SymLz(alpha2)) THEN
                    STOP 'Alphas not chosen correctly'
                ENDIF
            ENDIF

            !Now we need to do the same for the betas...
            !Find beta orbitals which will combine to make j
            !These need to be from the same L function, but have opposite Lz value
            IF(SymLz(j).eq.0) THEN
            !If Lz(j)=0, then beta = j. This will include all S states and Ml=0 states of higher angular momentum functions.

                beta1=j
                beta2=0
                beta1Coeff=DCMPLX(1.D0,0.D0)
                beta2Coeff=DCMPLX(0.D0,0.D0)

            ELSE
            !If Lz(j).ne.0, then there will be two contributary solid harmonic
            !beta1 will always be the orbital with positive ml value, and beta2 will have negative ml.

!                        IF(SymL(j).ne.1) THEN
                IF(SymLz(j).lt.0) THEN
                    !Negative Ml - this will be beta2
                    beta2=j
                    beta1=LzPairs(j)
!                            beta1=j+2*abs(SymLz(j))    !This will be the corresponding positive Ml state
                    !Now find coefficients. 
                    IF(mod(abs(SymLz(j)),2).eq.1) THEN
                        !We have an ODD ml pair of orbitals
                        Beta1Coeff=DCMPLX(Norm,0.D0)
                        Beta2Coeff=DCMPLX(0.D0,-Norm)
                    ELSE
                        !We have an EVEN ml pair of orbitals
                        Beta1Coeff=DCMPLX(Norm,0.D0)
                        Beta2Coeff=DCMPLX(0.D0,-Norm)
                    ENDIF
                ELSE
                    beta1=j
                    beta2=LzPairs(j)
!                            beta2=j-2*abs(SymLz(j))
                    IF(mod(abs(SymLz(j)),2).eq.1) THEN
                        !We have an ODD ml pair of orbitals
                        Beta1Coeff=DCMPLX(-Norm,0.D0)
                        Beta2Coeff=DCMPLX(0.D0,-Norm)
                    ELSE
                        !We have an EVEN ml pair of orbitals
                        Beta1Coeff=DCMPLX(Norm,0.D0)
                        Beta2Coeff=DCMPLX(0.D0,Norm)
                    ENDIF
                ENDIF
                IF(SymLz(beta1).ne.-SymLz(beta2)) THEN
                    STOP 'betas not chosen correctly'
                ENDIF
            ENDIF
                
!We now need to sum the contributions correctly!
!            IF(i.eq.j) THEN
!                WRITE(6,*) i,Alpha1Coeff,Beta1Coeff,alpha1,beta1,TMAT(alpha1,beta1)
!                WRITE(6,*) "***"
!                WRITE(6,*) Alpha2Coeff,Beta1Coeff,alpha2,beta1,TMAT(alpha2,beta1)
!                WRITE(6,*) "***"
!                WRITE(6,*) Alpha1Coeff,Beta2Coeff,alpha1,beta2,TMAT(alpha1,beta2)
!                WRITE(6,*) "***"
!                WRITE(6,*) Alpha2Coeff,Beta2Coeff,alpha2,beta2,TMAT(alpha2,beta2)
!                WRITE(6,*) "*******************"
!            ENDIF
            CompZ=Alpha1Coeff*Beta1Coeff*TMAT(alpha1,beta1)
            CompZ=CompZ+Alpha2Coeff*Beta1Coeff*TMAT(alpha2,beta1)
            CompZ=CompZ+Alpha1Coeff*Beta2Coeff*TMAT(alpha1,beta2)
            CompZ=CompZ+Alpha2Coeff*Beta2Coeff*TMAT(alpha2,beta2)

            IF(ABS(AIMAG(CompZ)).gt.1.D-9) THEN
                STOP 'Complex integral found'
            ENDIF
            IF(ABS(DBLE(CompZ)).gt.1.D-09) THEN
                WRITE(8,*) DBLE(CompZ),i,j,0,0
                !WRITE(8,'(1X,G20.14,4I3)') DBLE(CompZ),i,j,0,0
            ENDIF

        ENDDO
    ENDDO

    DO i=1,NORB
        WRITE(8,*) Arr(i),i,0,0,0
        !WRITE(8,'(1X,G20.12,4I3)') Arr(i),i,0,0,0
    ENDDO

    WRITE(8,*) ECore,0,0,0,0
    !WRITE(8,'(1X,G20.12,4I3)') ECore,0,0,0,0

    if(tError) then
        CLOSE(8,status='delete')
    else
        CLOSE(8)
    endif

    DEALLOCATE(UMAT)
    DEALLOCATE(ARR)
    DEALLOCATE(TMAT)
!    DEALLOCATE(SymLCounts)
!    DEALLOCATE(SymLList)
!    DEALLOCATE(IND)
        
END PROGRAM TransLz

FUNCTION UMatInd(I,J,K,L)
    IMPLICIT NONE
    INTEGER :: I,J,K,L,A,B,UMatInd

    IF(I.eq.0.or.J.eq.0.or.K.eq.0.or.L.eq.0) THEN
        UMatInd=0
        RETURN
    ENDIF

    IF(I.gt.J) THEN
        A=(I*(I-1))/2+J
    ELSE
        A=(J*(J-1))/2+I
    ENDIF

    IF(K.gt.L) THEN
        B=(K*(K-1))/2+L
    ELSE
        B=(L*(L-1))/2+K
    ENDIF
    IF(A.gt.B) THEN
        UMatInd=(A*(A-1))/2+B
    ELSE
        UMatInd=(B*(B-1))/2+A
    ENDIF

END FUNCTION UMatInd
