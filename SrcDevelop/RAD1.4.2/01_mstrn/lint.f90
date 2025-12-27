SUBROUTINE LINT(K,X,Y,N,A,B)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::K,N
  REAL(8),INTENT(IN)::X(K),Y(K),A(N)
  REAL(8),INTENT(OUT)::B(N)
  INTEGER::I,J,MIN,MAX
  REAL(8)::WGT
!--exec
  MIN=COUNT(A(:)<X(1))
  MAX=COUNT(A(:)>X(K))
  DO I=1,N
     IF(I<=MIN) THEN
        WGT=(X(1)-A(I))/(X(2)-X(1))
        B(I)=Y(1)-(Y(2)-Y(1))*WGT
     ELSEIF(I<=N-MAX) THEN
        DO J=1,K-1
           IF(X(J)<=A(I) .and. X(J+1)>A(I)) then
              WGT=(A(I)-X(J))/(X(J+1)-X(J))
              B(I)=Y(J)*(1-WGT)+Y(J+1)*WGT
              EXIT
           ENDIF
        ENDDO
     ELSE
        WGT=(A(I)-X(K))/(X(K)-X(K-1))
        B(I)=Y(K)+(Y(K)-Y(K-1))*WGT
     ENDIF
  ENDDO
END SUBROUTINE LINT

