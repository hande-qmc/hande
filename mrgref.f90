Module m_mrgref
Integer, Parameter :: kdp = selected_real_kind(15)
public :: mrgref
private :: kdp
private :: R_mrgref, I_mrgref, D_mrgref
interface mrgref
  module procedure d_mrgref, r_mrgref, i_mrgref
end interface mrgref
contains

Subroutine D_mrgref (XVALT, IRNGT)
!   Ranks array XVALT into index array IRNGT, using merge-sort
! __________________________________________________________
!   This version is not optimized for performance, and is thus
!   not as difficult to read as some other ones.
!   Michel Olagnon - April 2000
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (In) :: XVALT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
!
      Integer, Dimension (:), Allocatable :: JWRKT
      Integer :: LMTNA, LMTNC
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      If (NVAL <= 0) Then
         Return
      End If
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XVALT(IIND-1) < XVALT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      Allocate (JWRKT(1:NVAL))
      LMTNC = 2
      LMTNA = 2
!
!  Iteration. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
         IWRK = 0
!
!   Loop on merges of A and B into C
!
         Do
            IINDA = IWRKF
            IWRKD = IWRKF + 1
            IWRKF = IINDA + LMTNC
            JINDA = IINDA + LMTNA
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDB = JINDA
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B (no need to do anything)
!
            If (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(JINDA+1))) Then
               IWRK = IWRKF
               Cycle
            End If
!
!  One steps in the C subset, that we create in the final rank array
!
            Do
               If (IWRK >= IWRKF) Then
!
!  Make a copy of the rank array for next iteration
!
                  IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
                  Exit
               End If
!
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (IINDA < JINDA) Then
                  If (IINDB < IWRKF) Then
                     If (XVALT(IRNGT(IINDA+1)) > XVALT(IRNGT(IINDB+1))) &
                    & Then
                        IINDB = IINDB + 1
                        JWRKT (IWRK) = IRNGT (IINDB)
                     Else
                        IINDA = IINDA + 1
                        JWRKT (IWRK) = IRNGT (IINDA)
                     End If
                  Else
!
!  Only A still with unprocessed values
!
                     IINDA = IINDA + 1
                     JWRKT (IWRK) = IRNGT (IINDA)
                  End If
               Else
!
!  Only B still with unprocessed values
!
                  IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
                  IWRK = IWRKF
                  Exit
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!  Clean up
!
      Deallocate (JWRKT)
      Return
!
End Subroutine D_mrgref

Subroutine R_mrgref (XVALT, IRNGT)
!   Ranks array XVALT into index array IRNGT, using merge-sort
! __________________________________________________________
!   This version is not optimized for performance, and is thus
!   not as difficult to read as some other ones.
!   Michel Olagnon - April 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XVALT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
!
      Integer, Dimension (:), Allocatable :: JWRKT
      Integer :: LMTNA, LMTNC
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      If (NVAL <= 0) Then
         Return
      End If
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XVALT(IIND-1) < XVALT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      Allocate (JWRKT(1:NVAL))
      LMTNC = 2
      LMTNA = 2
!
!  Iteration. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
         IWRK = 0
!
!   Loop on merges of A and B into C
!
         Do
            IINDA = IWRKF
            IWRKD = IWRKF + 1
            IWRKF = IINDA + LMTNC
            JINDA = IINDA + LMTNA
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDB = JINDA
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B (no need to do anything)
!
            If (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(JINDA+1))) Then
               IWRK = IWRKF
               Cycle
            End If
!
!  One steps in the C subset, that we create in the final rank array
!
            Do
               If (IWRK >= IWRKF) Then
!
!  Make a copy of the rank array for next iteration
!
                  IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
                  Exit
               End If
!
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (IINDA < JINDA) Then
                  If (IINDB < IWRKF) Then
                     If (XVALT(IRNGT(IINDA+1)) > XVALT(IRNGT(IINDB+1))) &
                    & Then
                        IINDB = IINDB + 1
                        JWRKT (IWRK) = IRNGT (IINDB)
                     Else
                        IINDA = IINDA + 1
                        JWRKT (IWRK) = IRNGT (IINDA)
                     End If
                  Else
!
!  Only A still with unprocessed values
!
                     IINDA = IINDA + 1
                     JWRKT (IWRK) = IRNGT (IINDA)
                  End If
               Else
!
!  Only B still with unprocessed values
!
                  IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
                  IWRK = IWRKF
                  Exit
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!  Clean up
!
      Deallocate (JWRKT)
      Return
!
End Subroutine R_mrgref
Subroutine I_mrgref (XVALT, IRNGT)
!   Ranks array XVALT into index array IRNGT, using merge-sort
! __________________________________________________________
!   This version is not optimized for performance, and is thus
!   not as difficult to read as some other ones.
!   Michel Olagnon - April 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XVALT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
!
      Integer, Dimension (:), Allocatable :: JWRKT
      Integer :: LMTNA, LMTNC
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      If (NVAL <= 0) Then
         Return
      End If
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XVALT(IIND-1) < XVALT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      Allocate (JWRKT(1:NVAL))
      LMTNC = 2
      LMTNA = 2
!
!  Iteration. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
         IWRK = 0
!
!   Loop on merges of A and B into C
!
         Do
            IINDA = IWRKF
            IWRKD = IWRKF + 1
            IWRKF = IINDA + LMTNC
            JINDA = IINDA + LMTNA
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDB = JINDA
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B (no need to do anything)
!
            If (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(JINDA+1))) Then
               IWRK = IWRKF
               Cycle
            End If
!
!  One steps in the C subset, that we create in the final rank array
!
            Do
               If (IWRK >= IWRKF) Then
!
!  Make a copy of the rank array for next iteration
!
                  IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
                  Exit
               End If
!
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (IINDA < JINDA) Then
                  If (IINDB < IWRKF) Then
                     If (XVALT(IRNGT(IINDA+1)) > XVALT(IRNGT(IINDB+1))) &
                    & Then
                        IINDB = IINDB + 1
                        JWRKT (IWRK) = IRNGT (IINDB)
                     Else
                        IINDA = IINDA + 1
                        JWRKT (IWRK) = IRNGT (IINDA)
                     End If
                  Else
!
!  Only A still with unprocessed values
!
                     IINDA = IINDA + 1
                     JWRKT (IWRK) = IRNGT (IINDA)
                  End If
               Else
!
!  Only B still with unprocessed values
!
                  IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
                  IWRK = IWRKF
                  Exit
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!  Clean up
!
      Deallocate (JWRKT)
      Return
!
End Subroutine I_mrgref
end module m_mrgref
