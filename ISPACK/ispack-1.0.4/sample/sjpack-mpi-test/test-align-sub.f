************************************************************************
* ISPACK FORTRAN SUBROUTINE LIBRARY FOR SCIENTIFIC COMPUTING
* Copyright (C) 1998--2015 Keiichi Ishioka <ishioka@gfd-dennou.org>
*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
* 02110-1301 USA.
************************************************************************
      SUBROUTINE SUB(T,Q,G,WG,W)

      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mpif.h'
      PARAMETER(MM=170,JM=256,IM=512)            
      PARAMETER(NM=MM,NN=MM,IPOW=0)
      DIMENSION S((2*NN+1-MM)*MM+NN+1),G(IM,JM)
      DIMENSION SD((2*NN+1-MM)*MM+NN+1)
      DIMENSION IT(4),T(IM*6)      
      DIMENSION P(JM/2*(MM+4)),Q(JM/2*7),R((MM+1)*(2*NM-MM-1)+1)
      DIMENSION WS((2*NN+1-MM)*MM+NN+1),WG((IM+2)*JM),W((JM+1)*IM)      

      DO L=1,(2*NN+1-MM)*MM+NN+1
        CALL SJL2NM(MM,L,N,M)
        S(L)=1D0/((ABS(M)+1D0)*(N+1D0))
      END DO
      CALL BSCOPY((2*NN+1-MM)*MM+NN+1,S,SD)

      CALL MPI_INIT(IERR)
      CALL SJPINI(MM,NM,JM,JC,IM,P,R,IT,T)
      CALL SJTS2G(MM,NM,NN,IM,JC,S,G,IT,T,P,Q,R,WS,WG,W,IPOW)
      CALL SJPG2S(MM,NM,NN,IM,JC,S,G,IT,T,P,Q,R,WS,WG,W,IPOW)
      CALL MPI_FINALIZE(IERR)      

      EPS=1D-13
      
      INDEX=0
      DO L=1,(2*NN+1-MM)*MM+NN+1
        BUNBO=ABS(S(L))+ABS(SD(L))+1
        BUNSHI=ABS(S(L)-SD(L))
        IF(BUNSHI/BUNBO.GT.EPS) THEN
          INDEX=INDEX+1
          WRITE(6,'(I6,2F25.17)') L,SD(L),S(L)
        END IF
      END DO

      print *,'test:'
      IF(INDEX.EQ.0) THEN
        print *,'** OK'
      ELSE
        print *,'** Fail'
      END IF

      END
