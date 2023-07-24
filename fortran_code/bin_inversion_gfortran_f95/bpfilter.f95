!CC FILE: bpfilter.f
MODULE BPFILTER


   USE CONSTANTS, ONLY : PI, TWOPI
   IMPLICIT NONE
!   REAL, PARAMETER, PRIVATE :: PI = 3.14159265, TWOPI = 2.0 * PI


contains


   subroutine BandPassFilter(data1, NP, DT, n_but, i_pass, lowcorner, highcorner)
   implicit none
!c Apply a simple band pass filter to data
   real*4 :: data1(*)
   real*4 :: lowcorner, highcorner, dt
   integer :: np, n_but, i_pass
!Cf2py intent(out) data1
!Cf2py intent(in) data1,NP,DT,n_but,i_pass,lowcorner,highcorner
   if (lowcorner.gt.0.0) then
      call IIRFILT( DATA1, NP, 'BUTTER  ',n_but, 'BP      ', lowcorner, highcorner, DT, i_pass)
   else
      call IIRFILT( DATA1, NP, 'BUTTER  ',n_but, 'LP      ', lowcorner, highcorner, DT, i_pass)
   endif
   end subroutine BandPassFilter

!  Copyright 1990  Regents of the University of California
!
!
!  Author:  Dave Harris
!
!           Lawrence Livermore National Laboratory
!           L-205
!           P.O. Box 808
!           Livermore, CA  94550
!           USA
!
!           (415) 423-0617
!                                                               APPLY
!  Subroutine to apply an iir filter to a DATA1 sequence.
!    The filter is assumed to be stored as second order sections.
!    Filtering is in-place.
!    Zero-phase (forward and reverse) is an option.
!
!  Input Arguments:
!  ----------------
!
!    DATA1                           Array containing DATA1
!
!    NSAMPS                         Number of DATA1 samples
!
!    ZP                             Logical variable, true for
!                                     zero phase filtering, false
!                                     for single pass filtering
!
!    SN                             Numerator polynomials for second
!                                     order sections.
!
!    SD                             Denominator polynomials for second
!                                     order sections.
!
!    NSECTS                         Number of second-order sections
!
!  Output Arguments:
!  -----------------
!
!    DATA1                          DATA1 array (same as input)
!
!
   SUBROUTINE APPLY( DATA1, NSAMPS, ZP, SN, SD, NSECTS )
   IMPLICIT NONE
   REAL*4 :: SN(30), SD(30)
   REAL*4 :: OUTPUT,DATA1(*)
   LOGICAL :: ZP
   INTEGER :: NSAMPS, NSECTS, JPTR, J, I
   REAL*4 :: A1, A2, B0, B1, B2, X1, X2, Y1, Y2
   
   JPTR = 1
   DO J = 1, NSECTS

      X1 = 0.0
      X2 = 0.0
      Y1 = 0.0
      Y2 = 0.0
      B0 = SN(JPTR)
      B1 = SN(JPTR + 1)
      B2 = SN(JPTR + 2)
      A1 = SD(JPTR + 1)
      A2 = SD(JPTR + 2)

      DO I = 1, NSAMPS

         OUTPUT = B0 * DATA1(I) + B1 * X1 + B2 * X2
         OUTPUT = OUTPUT - ( A1 * Y1 + A2 * Y2 )
         Y2 = Y1
         Y1 = OUTPUT
         X2 = X1
         X1 = DATA1(I)
         DATA1(I) = OUTPUT

      ENDDO

      JPTR = JPTR + 3

   ENDDO

   IF (   ZP ) THEN

      JPTR = 1
      DO J = 1, NSECTS

         X1 = 0.0
         X2 = 0.0
         Y1 = 0.0
         Y2 = 0.0
         B0 = SN(JPTR)
         B1 = SN(JPTR + 1)
         B2 = SN(JPTR + 2)
         A1 = SD(JPTR + 1)
         A2 = SD(JPTR + 2)

         DO I = NSAMPS, 1, -1

            OUTPUT = B0 * DATA1(I) + B1 * X1 + B2 * X2
            OUTPUT = OUTPUT - ( A1 * Y1 + A2 * Y2 )
            Y2 = Y1
            Y1 = OUTPUT
            X2 = X1
            X1 = DATA1(I)
            DATA1(I) = OUTPUT

         ENDDO

         JPTR = JPTR + 3

      ENDDO

   END IF

   END SUBROUTINE APPLY


! BEROOTS -- SUBROUTINE TO RETURN BESSEL POLES FOR
!   NORMALIZED LOWPASS FILTER
!
! LAST MODIFIED:  April 15, 1992. Changed P and RTYPE to adjustable
!                 array by using an "*" rather than a "1".
!
!  OUTPUT ARGUMENTS:
!  -----------------
!      P              COMPLEX ARRAY CONTAINING POLES
!                       CONTAINS ONLY ONE FROM EACH
!                       COMPLEX CONJUGATE PAIR, AND
!                       ALL REAL POLES
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION
!                       TYPE:
!                         (SP)  SINGLE REAL POLE
!                         (CP)  COMPLEX CONJUGATE POLE PAIR
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS
!
!      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS
!
!  INPUT ARGUMENTS:
!  ----------------
!
!      IORD           DESIRED FILTER ORDER
!
!
   SUBROUTINE BEROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )
!
   IMPLICIT NONE
   COMPLEX :: P(*)
   INTEGER :: NSECTS, IORD
   CHARACTER(LEN=3) RTYPE(*)
   REAL*4 :: DCVALUE
!       !-- LUDAN TEST IORD
!
   IF (   IORD .EQ. 1 ) THEN

      P(1) = CMPLX( -1.0, 0.0 )
      RTYPE(1) = 'SP'

   ELSE IF (  IORD .EQ. 2 ) THEN

      P(1) = CMPLX( -1.1016013,  0.6360098 )
      RTYPE(1) = 'CP'

   ELSE IF (  IORD .EQ. 3 ) THEN

      P(1) = CMPLX( -1.0474091, 0.9992645 )
      RTYPE(1) = 'CP'
      P(2) = CMPLX( -1.3226758, 0.0 )
      RTYPE(2) = 'SP'

   ELSE IF (  IORD .EQ. 4 ) THEN

      P(1) = CMPLX( -0.9952088,  1.2571058 )
      RTYPE(1) = 'CP'
      P(2) = CMPLX( -1.3700679, 0.4102497 )
      RTYPE(2) = 'CP'

   ELSE IF (  IORD .EQ. 5 ) THEN

      P(1) = CMPLX( -0.9576766,  1.4711244 )
      RTYPE(1) = 'CP'
      P(2) = CMPLX( -1.3808774,  0.7179096 )
      RTYPE(2) = 'CP'
      P(3) = CMPLX( -1.5023160, 0.0 )
      RTYPE(3) = 'SP'

   ELSE IF (  IORD .EQ. 6 ) THEN

      P(1) = CMPLX( -0.9306565,  1.6618633 )
      RTYPE(1) = 'CP'
      P(2) = CMPLX( -1.3818581,  0.9714719 )
      RTYPE(2) = 'CP'
      P(3) = CMPLX( -1.5714904,  0.3208964 )
      RTYPE(3) = 'CP'

   ELSE IF (  IORD .EQ. 7 ) THEN

      P(1) = CMPLX( -0.9098678,  1.8364514 )
      RTYPE(1) = 'CP'
      P(2) = CMPLX( -1.3789032,  1.1915667 )
      RTYPE(2) = 'CP'
      P(3) = CMPLX( -1.6120388,  0.5892445 )
      RTYPE(3) = 'CP'
      P(4) = CMPLX( -1.6843682, 0.0 )
      RTYPE(4) = 'SP'

   ELSE IF (  IORD .EQ. 8 ) THEN

      P(1) = CMPLX( -0.8928710,  1.9983286 )
      RTYPE(1) = 'CP'
      P(2) = CMPLX( -1.3738431,  1.3883585 )
      RTYPE(2) = 'CP'
      P(3) = CMPLX( -1.6369417,  0.8227968 )
      RTYPE(3) = 'CP'
      P(4) = CMPLX( -1.7574108,  0.2728679 )
      RTYPE(4) = 'CP'

   END IF

   NSECTS = IORD - IORD/2
!      LUDAN TEST NSECTS


   DCVALUE = 1.0
!
!  DONE
!
   END SUBROUTINE BEROOTS
!
!
!  Transforms an analog filter to a digital filter via the bilinear transformati
!    Assumes both are stored as second order sections.  The transformation is
!    done in-place.
!
!  Input Arguments:
!  ----------------
!
!    SN                   Array containing numerator polynomial coefficients for
!                           second order sections.  Packed head-to-tail.
!
!    SD                   Array containing denominator polynomial coefficients f
!                           second order sections.  Packed head-to-tail.
!
!    NSECTS               Number of second order sections.
!
!
   SUBROUTINE BILIN2( SN, SD, NSECTS )

   IMPLICIT NONE
   REAL*4 :: SN(30), SD(30)
   INTEGER :: IPTR, NSECTS, I
   REAL*4 :: A0, A1, A2, SCALE1

   IPTR = 1
   DO I = 1, NSECTS

      A0 = SD(IPTR)
      A1 = SD(IPTR + 1)
      A2 = SD(IPTR + 2)

      SCALE1 = A2 + A1 + A0
      SD(IPTR)   = 1.
      SD(IPTR + 1) = (2.0 * (A0 - A2)) / SCALE1
      SD(IPTR + 2) = (A2 - A1 + A0) / SCALE1

      A0 = SN(IPTR)
      A1 = SN(IPTR + 1)
      A2 = SN(IPTR + 2)

      SN(IPTR)   = (A2 + A1 + A0) / SCALE1
      SN(IPTR + 1) = (2.0 * (A0 - A2)) / SCALE1
      SN(IPTR + 2) = (A2 - A1 + A0) / SCALE1

      IPTR = IPTR + 3

   ENDDO

   END SUBROUTINE BILIN2
!
!
! BUROOTS -- SUBROUTINE TO COMPUTE BUTTERWORTH POLES FOR
!   NORMALIZED LOWPASS FILTER
!
! LAST MODIFIED:  SEPTEMBER 7, 1990
!
!  OUTPUT ARGUMENTS:
!  -----------------
!      P              COMPLEX ARRAY CONTAINING POLES
!                       CONTAINS ONLY ONE FROM EACH
!                       COMPLEX CONJUGATE PAIR, AND
!                       ALL REAL POLES
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION
!                       TYPE:
!                         (SP)  SINGLE REAL POLE
!                         (CP)  COMPLEX CONJUGATE POLE PAIR
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS
!
!      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS
!
!  INPUT ARGUMENTS:
!  ----------------
!
!      IORD           DESIRED FILTER ORDER
!
!
   SUBROUTINE BUROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )

   IMPLICIT NONE
   COMPLEX :: P(10)
   REAL*4 :: DCVALUE
   INTEGER :: IORD, NSECTS
   CHARACTER(LEN=3) :: RTYPE(10)
   INTEGER :: HALF, K
   REAL*4 :: ANGLE


   HALF = IORD / 2
!
! TEST FOR ODD ORDER, AND ADD POLE AT -1
!
   NSECTS = 0
   IF (    2*HALF .LT. IORD ) THEN
      P(1) = CMPLX( -1., 0. )
      RTYPE(1) = 'SP'
      NSECTS = 1
   END IF

   DO K = 1, HALF
      ANGLE = PI * ( .5 + FLOAT(2 * K - 1) / FLOAT(2 * IORD) )
      NSECTS = NSECTS + 1
      P(NSECTS) = CMPLX( COS(ANGLE), SIN(ANGLE) )
      RTYPE(NSECTS) = 'CP'
   ENDDO

   DCVALUE = 1.0

   END SUBROUTINE BUROOTS



!
!  Subroutine to alter the cutoff of a filter.  Assumes that the
!    filter is structured as second order sections.  Changes
!    the cutoffs of normalized lowpass and highpass filters through
!    a simple polynomial transformation.
!
!  Input Arguments:
!  ----------------
!
!    F                       New cutoff frequency
!
!  Input/Output Arguments:
!  -----------------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!    NSECTS                  Number of second order sectionsects
!
!
   SUBROUTINE CUTOFFS( SN, SD, NSECTS, F )

   IMPLICIT NONE
   INTEGER :: NSECTS
   REAL*4 :: SN(30), SD(30), F
   INTEGER :: IPTR, I
   REAL*4 :: SCALE1

   SCALE1 = 2.*3.14159265*F

   IPTR = 1
   DO I = 1, NSECTS

      SN( IPTR + 1 ) = SN( IPTR + 1 ) / SCALE1
      SN( IPTR + 2 ) = SN( IPTR + 2 ) / (SCALE1*SCALE1)
      SD( IPTR + 1 ) = SD( IPTR + 1 ) / SCALE1
      SD( IPTR + 2 ) = SD( IPTR + 2 ) / (SCALE1*SCALE1)
      IPTR = IPTR + 3

   ENDDO

   END SUBROUTINE CUTOFFS
!
!
!  Subroutine to design IIR digital filters from analog prototypes.
!
!  Input Arguments:
!  ----------------
!
!    IORD                Filter order (10 MAXIMUM)
!
!    TYPE                Character*2 variable containing filter type
!                          LOWPASS (LP)
!                          HIGHPASS (HP)
!                          BANDPASS (BP)
!                          BANDREJECT (BR)
!
!    APROTO              Character*2 variable designating analog prototype
!                          Butterworth (BU)
!                          Bessel (BE)
!
!    FL                  Low-frequency cutoff
!
!    FH                  High-frequency cutoff
!
!    TS                  Sampling interval (in seconds)
!
!  Output Arguments:
!  -----------------
!
!    SN                  Array containing numerator coefficients of
!                        second-order sections packed head-to-tail.
!
!    SD                  Array containing denominator coefficients
!                        of second-order sections packed head-to-tail.
!
!    NSECTS              Number of second-order sections.
!
!
   SUBROUTINE DESIGN( IORD, FTYPE, APROTO, FL, FH, TS, SN, SD, NSECTS )
 
   IMPLICIT NONE
   COMPLEX :: P(10), Z(10)
   CHARACTER(LEN=2) :: FTYPE, APROTO
   CHARACTER(LEN=3) :: STYPE(10)
   REAL*4 :: SN(30), SD(30), FL, FH, TS
   INTEGER :: IORD, NSECTS
   REAL*4 :: FLW, FHW, DCVALUE
!
!  Analog prototype selection
!
   IF (     APROTO .EQ. 'BU' ) THEN

      CALL BUROOTS( P, STYPE, DCVALUE, NSECTS, IORD )

   ELSE IF (    APROTO .EQ. 'BE' ) THEN

      CALL BEROOTS( P, STYPE, DCVALUE, NSECTS, IORD )

   END IF
!
!  Analog mapping selection
!
   IF (     FTYPE .EQ. 'BP' ) THEN
         
      FLW = WARP( FL*TS/2., 2. )
      FHW = WARP( FH*TS/2., 2. )
      CALL LPTBP( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )

   ELSE IF (   FTYPE .EQ. 'BR' ) THEN

      FLW = WARP( FL*TS/2., 2. )
      FHW = WARP( FH*TS/2., 2. )
      CALL LPTBR( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )

   ELSE IF (   FTYPE .EQ. 'LP' ) THEN

      FHW = WARP( FH*TS/2., 2. )
      CALL LP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )
      CALL CUTOFFS( SN, SD, NSECTS, FHW )

   ELSE IF (   FTYPE .EQ. 'HP' ) THEN

      FLW = WARP( FL*TS/2., 2. )
      CALL LPTHP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )
      CALL CUTOFFS( SN, SD, NSECTS, FLW )

   END IF
!
!  Bilinear analog to digital transformation
!
   CALL BILIN2( SN, SD, NSECTS )

   END SUBROUTINE DESIGN








!  ARGUMENTS:
!  ----------
!
!    DATA1           REAL ARRAY CONTAINING SEQUENCE TO BE FILTERED
!                     ORIGINAL DATA1 DESTROYED, REPLACED BY FILTERED DATA1
!
!    NSAMPS         NUMBER OF SAMPLES IN DATA1
!
!
!    APROTO         CHARACTER*8 VARIABLE, CONTAINS TYPE OF ANALOG
!                     PROTOTYPE FILTER
!                     '(BU)TTER  ' -- BUTTERWORTH FILTER
!                     '(BE)SSEL  ' -- BESSEL FILTER
!
!    IORD           ORDER (#POLES) OF ANALOG PROTOTYPE
!                   NOT TO EXCEED 10 IN THIS CONFIGURATION.  4 - 5
!                   SHOULD BE AMPLE.
!
!    TYPE           CHARACTER*8 VARIABLE CONTAINING FILTER TYPE
!                     'LP' -- LOW PASS
!                     'HP' -- HIGH PASS
!                     'BP' -- BAND PASS
!                     'BR' -- BAND REJECT
!
!    FLO            LOW FREQUENCY CUTOFF OF FILTER (HERTZ)
!                   IGNORED IF TYPE = 'LP'
!
!    FHI            HIGH FREQUENCY CUTOFF OF FILTER (HERTZ)
!                   IGNORED IF TYPE = 'HP'
!
!    TS             SAMPLING INTERVAL (SECONDS)
!
!    PASSES           INTEGER VARIABLE CONTAINING THE NUMBER OF PASSES
!                   1 -- FORWARD FILTERING ONLY
!                   2 -- FORWARD AND REVERSE (I.E. ZERO PHASE) FILTERING
!
!
!  SUBPROGRAMS REFERENCED:  BILIN2, BUROOTS, WARP, CUTOFFS, LPTHP, LPTBP,
!    LP, LPTBR, BEROOTS, DESIGN, APPLY
!
   SUBROUTINE IIRFILT( DATA1, NSAMPS, APROTO, IORD, FTYPE, FLO, FHI, TS, PASSES )

   IMPLICIT NONE
   REAL*4 :: DATA1(*)
   CHARACTER(LEN=8) :: FTYPE, APROTO
   INTEGER :: NSAMPS, PASSES, IORD
   REAL*4 ::  FLO, FHI, TS, SN(30), SD(30)
   LOGICAL :: ZP
   INTEGER :: NSECTS
!
!  Filter designed
   CALL DESIGN( IORD, FTYPE(1:2), APROTO(1:2), FLO, FHI, TS, SN, SD, NSECTS )
!
!  Filter DATA1
!
   IF (   PASSES .EQ. 1 ) THEN
      ZP = .FALSE.
   ELSE
      ZP = .TRUE.
   END IF
   CALL APPLY( DATA1, NSAMPS, ZP, SN, SD, NSECTS )

   END SUBROUTINE IIRFILT
!
!  Copyright 1990  Regents of the University of California
!
!
!
!  Subroutine to generate second order section parameterization
!    from an pole-zero description for lowpass filters.
!
!  Input Arguments:
!  ----------------
!
!    P                       Array containing poles
!
!    Z                       Array containing zeros
!
!    RTYPE                   Character array containing root type information
!                              (SP)  Single real pole or
!                              (CP)  Complex conjugate pole pair
!                              (CPZ) Complex conjugate pole and zero pairs
!
!    DCVALUE                 Zero-frequency value of prototype filter
!
!    NSECTS                  Number of second-order sections
!
!  Output Arguments:
!  -----------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!
   SUBROUTINE LP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )

   IMPLICIT NONE
   COMPLEX :: P(*), Z(*)
   CHARACTER(LEN=3) :: RTYPE(*)
   REAL*4 :: SN(*), SD(*), DCVALUE
   REAL*4 :: SCALE1
   INTEGER :: IPTR, I, NSECTS

   IPTR = 1
   DO I = 1, NSECTS

      IF (   RTYPE(I) .EQ. 'CPZ' ) THEN

         SCALE1 = REAL( P(I) * CONJG( P(I) ) ) &
     &   / REAL( Z(I) * CONJG( Z(I) ) )
         SN( IPTR )     = REAL( Z(I) * CONJG( Z(I) ) ) * SCALE1
         SN( IPTR + 1 ) = -2. * REAL( Z(I) ) * SCALE1
         SN( IPTR + 2 ) = 1. * SCALE1
         SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )
         SD( IPTR + 1 ) = -2. * REAL( P(I) )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3

      ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN

         SCALE1 = REAL( P(I) * CONJG( P(I) ) )
         SN( IPTR )     = SCALE1
         SN( IPTR + 1 ) = 0.
         SN( IPTR + 2 ) = 0.
         SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )
         SD( IPTR + 1 ) = -2. * REAL( P(I) )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3

      ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN

         SCALE1 = -REAL( P(I) )
         SN( IPTR )     = SCALE1
         SN( IPTR + 1 ) = 0.
         SN( IPTR + 2 ) = 0.
         SD( IPTR )     = -REAL( P(I) )
         SD( IPTR + 1 ) = 1.
         SD( IPTR + 2 ) = 0.
         IPTR = IPTR + 3

      END IF

   ENDDO

   SN(1) = DCVALUE * SN(1)
   SN(2) = DCVALUE * SN(2)
   SN(3) = DCVALUE * SN(3)

   END SUBROUTINE LP


!
!  Subroutine to convert an prototype lowpass filter to a bandpass filter via
!    the analog polynomial transformation.  The lowpass filter is
!    described in terms of its poles and zeros (as input to this routine).
!    The output consists of the parameters for second order sections.
!
!  Input Arguments:
!  ----------------
!
!    P                       Array containing poles
!
!    Z                       Array containing zeros
!
!    RTYPE                   Character array containing type information
!                              (SP) single real pole  or
!                              (CP) complex conjugate pole pair  or
!                              (CPZ) complex conjugate pole/zero pairs
!
!    DCVALUE                 Zero frequency value of filter
!
!    NSECTS                  Number of second-order sections upon input
!
!    FL                      Low-frequency cutoff
!
!    FH                      High-frequency cutoff
!
!  Output Arguments:
!  -----------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!    NSECTS                  Number of second order sections upon output
!                              This subroutine doubles the number of
!                              sections.
!
!
   SUBROUTINE LPTBP( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )

   IMPLICIT NONE
   INTEGER :: NSECTS
   COMPLEX :: P(*), Z(*), CTEMP, P1, P2, Z1, Z2, S, H
   CHARACTER(LEN=3) :: RTYPE(*)
   REAL*4 :: SN(*), SD(*), DCVALUE, FL, FH
   REAL*4 :: A, B, SCALE1
   INTEGER :: IPTR, I, N

   A = TWOPI*TWOPI*FL*FH
   B = TWOPI*( FH - FL )

   N = NSECTS
   NSECTS = 0
   IPTR = 1
   DO I = 1, N

      IF (    RTYPE(I) .EQ. 'CPZ' ) THEN

         CTEMP = ( B*Z(I) )**2 - 4.*A
         CTEMP = CSQRT( CTEMP )
         Z1 = 0.5*( B*Z(I) + CTEMP )
         Z2 = 0.5*( B*Z(I) - CTEMP )
         CTEMP = ( B*P(I) )**2 - 4.*A
         CTEMP = CSQRT( CTEMP )
         P1 = 0.5*( B*P(I) + CTEMP )
         P2 = 0.5*( B*P(I) - CTEMP )
         SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )
         SN( IPTR + 1 ) = -2. * REAL( Z1 )
         SN( IPTR + 2 ) = 1.
         SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
         SD( IPTR + 1 ) = -2. * REAL( P1 )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3
         SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )
         SN( IPTR + 1 ) = -2. * REAL( Z2 )
         SN( IPTR + 2 ) = 1.
         SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
         SD( IPTR + 1 ) = -2. * REAL( P2 )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3

         NSECTS = NSECTS + 2

      ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN

         CTEMP = ( B*P(I) )**2 - 4.*A
         CTEMP = CSQRT( CTEMP )
         P1 = 0.5*( B*P(I) + CTEMP )
         P2 = 0.5*( B*P(I) - CTEMP )
         SN( IPTR )     = 0.
         SN( IPTR + 1 ) = B
         SN( IPTR + 2 ) = 0.
         SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
         SD( IPTR + 1 ) = -2. * REAL( P1 )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3
         SN( IPTR )     = 0.
         SN( IPTR + 1 ) = B
         SN( IPTR + 2 ) = 0.
         SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
         SD( IPTR + 1 ) = -2. * REAL( P2 )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3

         NSECTS = NSECTS + 2

      ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN

         SN( IPTR )     = 0.
         SN( IPTR + 1 ) = B
         SN( IPTR + 2 ) = 0.
         SD( IPTR )     = A
         SD( IPTR + 1 ) = -B*REAL( P(I) )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3

         NSECTS = NSECTS + 1

      END IF

   ENDDO
!
!  Scaling - use the fact that the bandpass filter amplitude at sqrt( omega_l *
!            equals the amplitude of the lowpass prototype at d.c.
!
   S = CMPLX( 0., SQRT(A) )
   H = CMPLX( 1., 0. )

   IPTR = 1
   DO I = 1, NSECTS
      H = H * ( ( SN(IPTR+2)*S + SN(IPTR+1) )*S + SN(IPTR) ) &
     &      / ( ( SD(IPTR+2)*S + SD(IPTR+1) )*S + SD(IPTR) )
      IPTR = IPTR + 3
   ENDDO
   SCALE1 = DCVALUE / SQRT( REAL( H ) * CONJG( H ) )

   SN(1) = SN(1) * SCALE1
   SN(2) = SN(2) * SCALE1
   SN(3) = SN(3) * SCALE1

   END SUBROUTINE LPTBP
!
!
!  Subroutine to convert a lowpass filter to a band reject filter
!    via an analog polynomial transformation.  The lowpass filter is
!    described in terms of its poles and zeros (as input to this routine).
!    The output consists of the parameters for second order sections.
!
!  Input Arguments:
!  ----------------
!
!    P                       Array containing poles
!
!    Z                       Array containing zeros
!
!    RTYPE                   Character array containing type information
!                              (SP)  single real pole or
!                              (CP)  complex conjugate pole pair
!                              (CPZ) complex conjugate pole/zero pairs
!
!    DCVALUE                 Zero-frequency value of prototype filter
!
!    NSECTS                  Number of second-order sections
!                              prior to transformation
!
!    FL                      Low-frequency cutoff
!
!    FH                      High-frequency cutoff
!
!  Output Arguments:
!  -----------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!    NSECTS                  Number of second order sections following
!                              transformation.  The number is doubled.
!
!
   SUBROUTINE LPTBR( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )

   IMPLICIT NONE
   COMPLEX :: P(*), Z(*), CINV, CTEMP, P1, P2, Z1, Z2
   CHARACTER(LEN=3) :: RTYPE(*)
   REAL*4 :: SN(*), SD(*), FL, FH
   REAL*4 :: H, A, B, SCALE1, DCVALUE
   INTEGER :: IPTR, I, NSECTS, N

   A = TWOPI*TWOPI*FL*FH
   B = TWOPI*( FH - FL )

   N = NSECTS
   NSECTS = 0
   IPTR = 1
   DO I = 1, N

      IF (    RTYPE(I) .EQ. 'CPZ' ) THEN

         CINV = 1./Z(I)
         CTEMP = ( B*CINV )**2 - 4.*A
         CTEMP = CSQRT( CTEMP )
         Z1 = 0.5*( B*CINV + CTEMP )
         Z2 = 0.5*( B*CINV - CTEMP )
         CINV = 1./P(I)
         CTEMP = ( B*CINV )**2 - 4.*A
         CTEMP = CSQRT( CTEMP )
         P1 = 0.5*( B*CINV + CTEMP )
         P2 = 0.5*( B*CINV - CTEMP )
         SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )
         SN( IPTR + 1 ) = -2. * REAL( Z1 )
         SN( IPTR + 2 ) = 1.
         SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
         SD( IPTR + 1 ) = -2. * REAL( P1 )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3
         SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )
         SN( IPTR + 1 ) = -2. * REAL( Z2 )
         SN( IPTR + 2 ) = 1.
         SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
         SD( IPTR + 1 ) = -2. * REAL( P2 )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3

         NSECTS = NSECTS + 2

      ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN

         CINV = 1./P(I)
         CTEMP = ( B*CINV )**2 - 4.*A
         CTEMP = CSQRT( CTEMP )
         P1 = 0.5*( B*CINV + CTEMP )
         P2 = 0.5*( B*CINV - CTEMP )
         SN( IPTR )     = A
         SN( IPTR + 1 ) = 0.
         SN( IPTR + 2 ) = 1.
         SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
         SD( IPTR + 1 ) = -2. * REAL( P1 )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3
         SN( IPTR )     = A
         SN( IPTR + 1 ) = 0.
         SN( IPTR + 2 ) = 1.
         SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
         SD( IPTR + 1 ) = -2. * REAL( P2 )
         SD( IPTR + 2 ) = 1.
         IPTR = IPTR + 3

         NSECTS = NSECTS + 2

      ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN

         SN( IPTR )     = A
         SN( IPTR + 1 ) = 0.
         SN( IPTR + 2 ) = 1.
         SD( IPTR )     = -A*REAL( P(I) )
         SD( IPTR + 1 ) = B
         SD( IPTR + 2 ) = -REAL( P(I) )
         IPTR = IPTR + 3

         NSECTS = NSECTS + 1

      END IF

   ENDDO
!
!  Scaling - use the fact that the bandreject filter amplitude  at d.c.
!            equals the lowpass prototype amplitude at d.c.
!
   H = 1.0

   IPTR = 1
   DO I = 1, NSECTS
      H = H * SN(IPTR) / SD(IPTR)
      IPTR = IPTR + 3
   ENDDO
   SCALE1 = DCVALUE / ABS(H)
   SN(1) = SN(1) * SCALE1
   SN(2) = SN(2) * SCALE1
   SN(3) = SN(3) * SCALE1

   END SUBROUTINE LPTBR
!
!
!  Subroutine to convert a lowpass filter to a highpass filter via
!    an analog polynomial transformation.  The lowpass filter is
!    described in terms of its poles and zeroes (as input to this routine).
!    The output consists of the parameters for second order sections.
!
!  Input Arguments:
!  ----------------
!
!    P                       Array containing poles
!
!    Z                       Array containing zeroes
!
!    RTYPE                   Character array containing root type information
!                              (SP) single real pole or
!                              (CP)  complex conjugate pair
!                              (CPZ) complex pole/zero pairs
!
!    DCVALUE                 Zero-frequency value of prototype filter
!
!    NSECTS                  Number of second-order sections
!
!  Output Arguments:
!  -----------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!
   SUBROUTINE LPTHP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )

   IMPLICIT NONE
   COMPLEX :: P(*), Z(*)
   CHARACTER(LEN=3) :: RTYPE(*)
   REAL*4 :: SN(*), SD(*), DCVALUE, SCALE1
   INTEGER :: NSECTS
   INTEGER :: IPTR, I

   IPTR = 1
   DO I = 1, NSECTS

      IF (     RTYPE(I) .EQ. 'CPZ' ) THEN

         SCALE1 = REAL( P(I) * CONJG( P(I) ) ) &
     &            / REAL( Z(I) * CONJG( Z(I) ) )
         SN( IPTR )     = 1.  *  SCALE1
         SN( IPTR + 1 ) = -2. * REAL( Z(I) )  *  SCALE1
         SN( IPTR + 2 ) = REAL( Z(I) * CONJG( Z(I) ) )  *  SCALE1
         SD( IPTR )     = 1.
         SD( IPTR + 1 ) = -2. * REAL( P(I) )
         SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )
         IPTR = IPTR + 3

      ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN

         SCALE1 = REAL( P(I) * CONJG( P(I) ) )
         SN( IPTR )     = 0.
         SN( IPTR + 1 ) = 0.
         SN( IPTR + 2 ) = SCALE1
         SD( IPTR )     = 1.
         SD( IPTR + 1 ) = -2. * REAL( P(I) )
         SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )
         IPTR = IPTR + 3

      ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN

         SCALE1 = -REAL( P(I) )
         SN( IPTR )     = 0.
         SN( IPTR + 1 ) = SCALE1
         SN( IPTR + 2 ) = 0.
         SD( IPTR )     = 1.
         SD( IPTR + 1 ) = -REAL( P(I) )
         SD( IPTR + 2 ) = 0.
         IPTR = IPTR + 3

      END IF

   ENDDO

   SN(1) = SN(1) * DCVALUE
   SN(2) = SN(2) * DCVALUE
   SN(3) = SN(3) * DCVALUE

   END SUBROUTINE LPTHP
!
!
! WARP -- FUNCTION, APPLIES TANGENT FREQUENCY WARPING TO COMPENSATE
!         FOR BILINEAR ANALOG -> DIGITAL TRANSFORMATION
!
! ARGUMENTS:
! ----------
!
!      F       ORIGINAL DESIGN FREQUENCY SPECIFICATION (HERTZ)
!      TS      SAMPLING INTERVAL (SECONDS)
!
!  LAST MODIFIED:  SEPTEMBER 20, 1990
!
   FUNCTION WARP( F , TS ) RESULT(ANSWER)

   IMPLICIT NONE
   REAL*4 :: ANSWER, F, TS, ANGLE
   ANGLE = TWOPI*F*TS/2.
   ANSWER = 2.*TAN(ANGLE)/TS
   ANSWER = ANSWER/TWOPI

   END FUNCTION WARP


END MODULE BPFILTER
