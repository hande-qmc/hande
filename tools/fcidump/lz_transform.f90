module lz_transform
   use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit, iostat_end

   implicit none

   integer, parameter :: i0 = selected_int_kind(15)
   integer, parameter :: p = selected_real_kind(15,307)
   real(p), parameter :: print_thres = 1.e-8_p ! Print integrals with real parts larger than this
   real(p), parameter :: zero_thres = 1.e-10_p ! Warn about integrals that shouldn't exist when they're larger than this
   real(p), parameter :: allow_thres = print_thres ! Set the above to zero instead of aborting if they're less than this
   real(p), parameter :: NORM = 1.0_p/sqrt(2.0_p)

   integer, parameter :: bignum = 1000

   contains

      elemental function eri_ind(i, j) result(ind)
         ! This function can be a general version for i,j,a,b, but 
         ! we can also use storage to cut down on unnecessary instructrions:
         ! i.e. by storing intermediate results ij and ab

         integer, intent(in) :: i, j
         integer :: ind

         if (i == 0 .or. j == 0) then
             ind = 0
             return
         end if

         if (i >= j) then
            ind = i*(i-1)/2 + j
         else
            ind = j*(j-1)/2 + i
         end if

      end function eri_ind

      subroutine get_lz_idx_and_coeff(lz, i, pair, conj, a1, a2, a1c, a2c)

         ! Find out which orbitals contribute to an Lz-transformed orbital,
         ! and their coefficients.
         ! In:
         !    lz: The signed Ml value.
         !    i: The (1-indexed) index of the current real orbital
         !    pair: The (1-indexed) index of the Lz-paired orbital of the current real orbital/
         !    conj: whether we should conjugate this Lz-transformed orbital: 1 for no, -1 for yes
         ! Out:
         !    a1: The (1-indexed) index of the x-like (positive Ml) contributing orbital.
         !    a2: The (1-indexed) index of the y-like (negative Ml) contributing orbital, zero if lz=0.
         !    a1c: The coefficient of a1.
         !    a2c: The coefficient of a2.

         integer, intent(in) :: lz, i, pair, conj
         integer, intent(out) :: a1, a2
         complex(p), intent(out) :: a1c, a2c

         if (lz == 0) then
            a1 = i
            a2 = 0
            a1c=cmplx(1.0_p,0.0_p)
            a2c=cmplx(0.0_p,0.0_p)
         else
            if (lz < 0) then
               a1 = pair
               a2 = i
               ! odd and even are the same coefficients
               a1c = cmplx(NORM,0.0_p)
               a2c = cmplx(0.0_p,-NORM*conj)
            else
               a1 = i
               a2 = pair
               if (mod(abs(lz), 2) == 1) then
                  ! Odd Ml
                  a1c = cmplx(-NORM,0.0_p)
                  a2c = cmplx(0.0_p,-NORM*conj)
               else
                  ! Evem Ml
                  a1c = cmplx(NORM,0.0_p)
                  a2c = cmplx(0.0_p,NORM*conj)
               end if
            end if
         end if
      end subroutine get_lz_idx_and_coeff

      subroutine write_transformed_intgrl(filename)

         character(*), intent(in) :: filename

         character(255) :: lz_filename
         integer :: ir, ierr, ios, size
         integer :: i, j, k, l, ij, kl, i1, i2, j1, j2, k1, k2, l1, l2
         integer :: NORB, NELEC, MS2, ORBSYM(bignum), ISYM, SYML(bignum), SYMLZ(bignum), PAIR(bignum)
         logical :: UHF
         real(p) :: intgrl, e_nuc
         real(p), allocatable :: tei_store(:), oei_store(:,:), mo_energies(:)
         complex(p) :: lzintgrl, i1c, i2c, j1c, j2c, k1c, k2c, l1c, l2c

         namelist /FCI/ NORB, NELEC, MS2, ORBSYM, ISYM, UHF, SYML, SYMLZ, PAIR

         write(stdout, *) 'Reading in '//trim(filename)

         open(newunit=ir, file=trim(filename), status='old', form='formatted')
         ! Read the namelist and advance to the line after the namelist
         read(ir, FCI)

         size = NORB*(NORB+1)/2
         size = size*(size+1)/2

         ! If an 0 index is given as any of i,j,k,l, zero is returned
         allocate(tei_store(0:size), source=0.0_p, stat=ierr)
         if (ierr /= 0) then
            write(stderr, *)'Cannot allocate tei_store of size ',size+1,' aborting!'
            stop
         end if

         ! If an 0 index is given as any of i,j, zero is returned
         allocate(oei_store(0:NORB, 0:NORB), source=0.0_p, stat=ierr)
         if (ierr /= 0) then
            write(stderr, *)'Cannot allocate oei_store of size ',NORB**2,' aborting!'
            stop
         end if

         allocate(mo_energies(NORB), source=0.0_p, stat=ierr)
         if (ierr /= 0) then
            write(stderr, *)'Cannot allocate mo_energies of size ',NORB,' aborting!'
            stop
         end if

         do
            read(ir, *, iostat=ios) intgrl, i, j, k, l
            ! If EOF reached
            if (ios==iostat_end) exit

            if (i /= 0 .and. j /= 0 .and. k /= 0 .and. l /= 0) then
               ! Two-electron integral
               ! 8-fold permutational symmetry, we use compression to save space
               ij = eri_ind(i, j)
               kl = eri_ind(k, l)
               tei_store(eri_ind(ij, kl)) = intgrl
            else if (i /= 0 .and. j /= 0 .and. k == 0 .and. l == 0) then
               ! One-electron integral
               oei_store(i, j) = intgrl
               oei_store(j, i) = intgrl
            else if (i /= 0 .and. j == 0 .and. k == 0 .and. l == 0) then
               ! MO energies
               mo_energies(i) = intgrl
            else if (i == 0 .and. j == 0 .and. k == 0 .and. l == 0) then
               ! Nuclear repulsion
               e_nuc = intgrl
            end if
         end do

         close(ir)

         write(stdout, *) 'Read-in done!'

         lz_filename = 'lz-'//trim(filename)

         write(stdout, *) 'Starting Lz integral transformation, storing in '//trim(lz_filename)

         write(stdout, *) '2e integrals..'

         open(newunit=ir, file=trim(lz_filename), status='replace', form='formatted')

         ! Write the namelist 
         write(ir, '(A,/,A,I0,A)') '&FCI','NORB=',NORB,','
         write(ir, '(A,I0,A,/,A,I0,A)') 'NELEC=',NELEC,',','MS2=',MS2,','
         write(ir, '(A)', advance='no') 'ORBSYM='
         do i = 1, NORB
            write(ir, '(I0,A)', advance='no') ORBSYM(i),','
         end do
         write(ir, *)
         write(ir, '(A,I0,A,/,A)') 'ISYM=',1,',','UHF=.FALSE.,'
         write(ir, '(A)',advance='no') 'SYML='
         do i = 1, NORB
            write(ir, '(I0,A)', advance='no') SYML(i),','
         end do
         write(ir, *)
         write(ir, '(A)',advance='no') 'SYMLZ='
         do i = 1, NORB
            write(ir, '(I0,A)', advance='no') SYMLZ(i),','
         end do
         write(ir, *)
         write(ir, '(A)') '&END'


         do i = 1, NORB
            ! For the given Lz, find out which two (if Lz is 0 then only one component)
            ! orbitals contribute, and then get their complex coefficients
            ! E.g. C(1,+1) \propto -1/sqrt(2) * (x+iy)
            ! And remember in Chemists' notation, (ij|kl) means i and k are complex conjugated
            call get_lz_idx_and_coeff(SYMLZ(i), i, pair(i), -1, i1, i2, i1c, i2c) ! Conjugate
            do j = 1, NORB
               ij = (i*(i-1))/2 + j
               call get_lz_idx_and_coeff(SYMLZ(j), j, pair(j), 1, j1, j2, j1c, j2c)
               do k = 1, NORB
                  call get_lz_idx_and_coeff(SYMLZ(k), k, pair(k), -1, k1, k2, k1c, k2c) ! Conjugate
                  do l = 1, NORB
                     kl = (k*(k-1))/2 + l
                     ! We have four-fold permutational symmetry here:
                     ! (ij|kl) = (ji|lk) = (kl|ij) = (lk|ji) 
                     ! instead of the usual 8-fold symmetry, as the Lz-transformed orbitals are complex.
                     ! (Note that this is also not the usual 4-fold symmetry of i<j and k<l, and 
                     ! is instead a symmetry under simultaneous transpotition)
                     ! The following checks will produce the correct indices, but does not quite take
                     ! care of simultaneous transpositions (which would be quite verbose to achieve)
                     ! so we're producing a superset of the permutationally unique indices.
                     ! E.g. the following checks will allow both (11|21) and (11|12)
                     if (kl < ij) cycle
                     if ((i < j) .and. (k < l)) cycle
                     if ((i > j) .and. (k < l)) cycle

                     call get_lz_idx_and_coeff(SYMLZ(l), l, pair(l), 1, l1, l2, l1c, l2c)
                     
                     ! Since every orbital is made up of up to 2 +- Ml components,
                     ! we're going to get up to 16 contributions to the 2e integral.
                     ! The zero integrals are taken care of by the zeroth element
                     ! of 'teint', and the fact that eri_ind returns a zero 
                     ! if any argument is zero.
                     lzintgrl =            i1c*j1c*k1c*l1c*tei_store(eri_ind(eri_ind(i1,j1),eri_ind(k1,l1)))
                     lzintgrl = lzintgrl + i2c*j1c*k1c*l1c*tei_store(eri_ind(eri_ind(i2,j1),eri_ind(k1,l1)))
                     lzintgrl = lzintgrl + i1c*j2c*k1c*l1c*tei_store(eri_ind(eri_ind(i1,j2),eri_ind(k1,l1)))
                     lzintgrl = lzintgrl + i2c*j2c*k1c*l1c*tei_store(eri_ind(eri_ind(i2,j2),eri_ind(k1,l1)))
                     lzintgrl = lzintgrl + i1c*j1c*k2c*l1c*tei_store(eri_ind(eri_ind(i1,j1),eri_ind(k2,l1)))
                     lzintgrl = lzintgrl + i2c*j1c*k2c*l1c*tei_store(eri_ind(eri_ind(i2,j1),eri_ind(k2,l1)))
                     lzintgrl = lzintgrl + i1c*j2c*k2c*l1c*tei_store(eri_ind(eri_ind(i1,j2),eri_ind(k2,l1)))
                     lzintgrl = lzintgrl + i2c*j2c*k2c*l1c*tei_store(eri_ind(eri_ind(i2,j2),eri_ind(k2,l1)))
                     lzintgrl = lzintgrl + i1c*j1c*k1c*l2c*tei_store(eri_ind(eri_ind(i1,j1),eri_ind(k1,l2)))
                     lzintgrl = lzintgrl + i2c*j1c*k1c*l2c*tei_store(eri_ind(eri_ind(i2,j1),eri_ind(k1,l2)))
                     lzintgrl = lzintgrl + i1c*j2c*k1c*l2c*tei_store(eri_ind(eri_ind(i1,j2),eri_ind(k1,l2)))
                     lzintgrl = lzintgrl + i2c*j2c*k1c*l2c*tei_store(eri_ind(eri_ind(i2,j2),eri_ind(k1,l2)))
                     lzintgrl = lzintgrl + i1c*j1c*k2c*l2c*tei_store(eri_ind(eri_ind(i1,j1),eri_ind(k2,l2)))
                     lzintgrl = lzintgrl + i2c*j1c*k2c*l2c*tei_store(eri_ind(eri_ind(i2,j1),eri_ind(k2,l2)))
                     lzintgrl = lzintgrl + i1c*j2c*k2c*l2c*tei_store(eri_ind(eri_ind(i1,j2),eri_ind(k2,l2)))
                     lzintgrl = lzintgrl + i2c*j2c*k2c*l2c*tei_store(eri_ind(eri_ind(i2,j2),eri_ind(k2,l2)))

                     ! Check for conservation of angular momentum: \iint a*(1)b(1) 1/r12 c*(2)d(2) dr1 dr2
                     ! i.e. <ac|bd>, meaning transition from ac->bd, so Ml(a)+Ml(c) must be equal to Ml(b)+Ml(d)
                     if (((symlz(i) + symlz(k)) /= (symlz(j)+symlz(l))) .and. (abs(lzintgrl) > zero_thres)) then
                        if (abs(lzintgrl) < allow_thres) then
                           write(stderr, *) i,j,k,l,' should be zero by conservation of angular momentum, ',&
                              'but has value ',lzintgrl,', which is within machine precision, setting to zero!'
                           lzintgrl = cmplx(0.0_p, 0.0_p)
                        else
                           write(stderr, *) i,j,k,l,' should be zero by conservation of angular momentum, ',&
                              'but has value ',lzintgrl,', which is outside machine precision, aborting!'
                           stop
                        end if
                     end if

                     ! Check whether the integral is strictly real
                     if (abs(aimag(lzintgrl)) > allow_thres) then
                        write(stderr, *)i,j,k,l,', like all other integrals, should have zero imaginary part, ',&
                           'but has imaginary part ',aimag(lzintgrl),', which is beyond machine precision, aborting!'
                        stop
                     end if

                     ! Print out if larger than threshold
                     if (abs(real(lzintgrl, kind=p)) > print_thres) then
                        write(ir, '(1X, ES28.20, 4I4)') real(lzintgrl,kind=p), i, j, k, l
                     end if
                  end do
               end do
            end do
         end do

         write(stdout, *) '1e integrals...'
         do i = 1, NORB
            call get_lz_idx_and_coeff(SYMLZ(i), i, pair(i), -1, i1, i2, i1c, i2c) ! Conjugate
            do j = 1, NORB
               call get_lz_idx_and_coeff(SYMLZ(j), j, pair(j), 1, j1, j2, j1c, j2c)

               lzintgrl =          + i1c*j1c*oei_store(i1,j1)
               lzintgrl = lzintgrl + i2c*j1c*oei_store(i2,j1)
               lzintgrl = lzintgrl + i1c*j2c*oei_store(i1,j2)
               lzintgrl = lzintgrl + i2c*j2c*oei_store(i2,j2)

               ! Check for conservation of angular momentum
               if ((symlz(i) /= symlz(j)) .and. (abs(lzintgrl) > zero_thres)) then
                  if (abs(lzintgrl) < allow_thres) then
                     write(stderr, *) i,j,' should be zero by conservation of angular momentum, ',&
                        'but has value ',lzintgrl,', which is within machine precision, setting to zero!'
                     lzintgrl = cmplx(0.0_p, 0.0_p)
                  else
                     write(stderr, *) i,j,' should be zero by conservation of angular momentum, ',&
                        'but has value ',lzintgrl,', which is outside machine precision, aborting!'
                     stop
                  end if
               end if

               ! Check whether the integral is strictly real
               if (abs(aimag(lzintgrl)) > allow_thres) then
                  write(stderr, *) i,j,', like all other integrals, should have zero imaginary part, ',&
                     'but has imaginary part ',aimag(lzintgrl),', which is beyond machine precision, aborting!'
                  stop
               end if

               ! Print out if larger than threshold
               if (abs(real(lzintgrl, kind=p)) > print_thres) then
                  write(ir, '(1X, ES28.20, 4I4)') real(lzintgrl,kind=p), i, j, 0, 0
               end if
            end do
         end do

         do i = 1, NORB
            write(ir, '(1X, ES28.20, 4I4)') mo_energies(i), i, 0, 0, 0
         end do

         write(ir, '(1X, ES28.20, 4I4)') e_nuc, 0, 0, 0, 0

         close(ir)

         write(stdout, *) 'Fortran transform done!'
      end subroutine write_transformed_intgrl

end module lz_transform

program lz_fcidump

   use lz_transform

   implicit none

   character(255) :: filename

   if (command_argument_count() /= 1) then
      write(stderr, *)'Please provide only the filename as the input argument!'
      stop
   end if

   call get_command_argument(1, filename)

   call write_transformed_intgrl(trim(filename))

end program lz_fcidump
