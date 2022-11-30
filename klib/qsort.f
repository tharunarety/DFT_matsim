      SUBROUTINE QSORT(A,IA,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(*),IA(*)
C ---------------------------------------------------------------------
C *************  RISOE COMPUTER LIBRARY          MAY 1980  ************
C ***********  PETER KIRKEGAARD,  RISOE NATIONAL LABORATORY  ***********
C ----------------------------------------------------------------------
C PURPOSE
C     SORTS THE ELEMENTS OF A NUMERICAL ARRAY IN ASCENDING ORDER.
C ----------------------------------------------------------------------
C METHOD
C     THIS ROUTINE IS BASED ON VAN EMDEN'S IMPROVED QUICKSORT ALGORITHM
C     WITH THE ADDITIONAL FEATURE THAT THE ORIGINAL INDICES OF THE SOR-
C     TED ELEMENTS ARE SAVED.
C ----------------------------------------------------------------------
C DESCRIPTION OF PARAMETERS       (IN=INPUT,OUT=OUTPUT,I/O=INPUT/OUTPUT)
C  A      I/O  ARRAY TO BE SORTED. ON EXIT THE ELEMENTS OF A WILL BE IN
C              ASCENDING ORDER.
C  IA     OUT  ARRAY KEEPING TRACK OF THE ORIGINAL INDICES. ON EXIT
C              IA(I) IS THE OLD INDEX OF THE ELEMENT APPEARING AT INDEX
C              I AFTER SORTING.
C  N      IN   NUMBER OF ELEMENTS TO BE SORTED.
C ----------------------------------------------------------------------
C REFERENCES
C     (1) M.H.VAN EMDEN, ALGORITHM 402, INCREASING THE EFFICIENCY OF
C         QUICKSORT, C.ACM 13(1970)563-567.
C ----------------------------------------------------------------------
C     HLS: 18-JAN-87
C ----------------------------------------------------------------------
      DIMENSION ILT(20),IUT(20)                                         
C
      DO 1 I=1,N
    1 IA(I)=I
      IF(N.EQ.1)GO TO 19
      IL1=1
      IU1=N
      ISTACK=0
    2 ISTACK=ISTACK+1
      IL=IL1
      IU=IU1
    3 IP=IL
      IQ=IU
      X=A(IP)
      Z=A(IQ)
      I=0
      J=IQ-IP-1
      IF(X.LE.Z)GO TO 4
      XSAV=X
      X=Z
      Z=XSAV
      A(IP)=X
      A(IQ)=Z
      ISAV=IA(IP)
      IA(IP)=IA(IQ)
      IA(IQ)=ISAV
    4 IF(IU-IL.LE.1)GO TO 5
      XX=X
      IX=IP
      ZZ=Z
      IZ=IQ
    6 IP=IP+1
      IF(IP.GE.IQ)GO TO 7
      X=A(IP)
      IF(X-XX)6,8,8
    7 IP=IQ-1
      GO TO 9
    8 IQ=IQ-1
      IF(IP.GE.IQ)GO TO 10
      Z=A(IQ)
      IF(Z-ZZ)11,11,8
   10 IQ=IP
      IP=IP-1
      Z=X
      X=A(IP)
   11 IF(X.LE.Z)GO TO 12
      XSAV=X
      X=Z
      Z=XSAV
      A(IP)=X
      A(IQ)=Z
      ISAV=IA(IP)
      IA(IP)=IA(IQ)
      IA(IQ)=ISAV
   12 IF(X.LE.XX)GO TO 13
      XX=X
      I=I+1
      IX=IP
   13 IF(Z.GE.ZZ)GO TO 6
      ZZ=Z
      I=I+1
      IZ=IQ
      GO TO 6
    9 IF(IP.EQ.IX.OR.X.EQ.XX)GO TO 14
      A(IP)=XX
      A(IX)=X
      ISAV=IA(IX)
      IA(IX)=IA(IP)
      IA(IP)=ISAV
   14 IF(IQ.EQ.IZ.OR.Z.EQ.ZZ)GO TO 15
      A(IQ)=ZZ
      A(IZ)=Z
      ISAV=IA(IZ)
      IA(IZ)=IA(IQ)
      IA(IQ)=ISAV
   15 IF(IU-IQ.LE.IP-IL)GO TO 16
      IL1=IL
      IU1=IP-1
      IL=IQ+1
      GO TO 17
   16 IU1=IU
      IL1=IQ+1
      IU=IP-1
   17 IF(I.EQ.J)GO TO 5
      IF(IU1.LE.IL1)GO TO 18
      ILT(ISTACK)=IL
      IUT(ISTACK)=IU
      GO TO 2
   18 IF(IU.GT.IL)GO TO 3
    5 ISTACK=ISTACK-1
      IF(ISTACK)19,19,20
   20 IL=ILT(ISTACK)
      IU=IUT(ISTACK)
      IF(IU-IL)5,5,3
   19 RETURN
C
      END
