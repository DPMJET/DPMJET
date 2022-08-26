
      DOUBLE PRECISION FUNCTION DT_CSTRMA(p1, p2)

      IMPLICIT NONE
      DOUBLE PRECISION am1, am2, p1, p2, ptoch, ech
      DIMENSION p1(4), p2(4)

      ptoch = SQRT((p1(1)+p2(1))**2+(p1(2)+p2(2))
     &              **2+(p1(3)+p2(3))**2)
      ech = p1(4) + p2(4)
      DT_CSTRMA = (ech+ptoch)*(ech-ptoch)
      RETURN
      END