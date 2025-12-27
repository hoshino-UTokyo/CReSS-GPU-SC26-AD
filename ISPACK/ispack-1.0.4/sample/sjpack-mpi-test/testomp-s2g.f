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
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mpif.h'
*      PARAMETER(MM=682,JM=1024,IM=2048)
*      PARAMETER(MM=341,JM=512,IM=1024)
      PARAMETER(MM=170,JM=256,IM=512)            
      PARAMETER(NM=MM,NN=MM,IPOW=0)
      PARAMETER(NT=32)
      DIMENSION S((2*NN+1-MM)*MM+NN+1),G(IM,JM)
      DIMENSION SD((2*NN+1-MM)*MM+NN+1),GD(IM,JM)
      DIMENSION IT(4),T(IM*6)      
      DIMENSION P(JM/2*(MM+4)),Q(JM/2*7*NT),R((MM+1)*(2*NM-MM-1)+1)
      DIMENSION PP(JM/2*(MM+4))
      DIMENSION WS((2*NN+1-MM)*MM+NN+1+2*(NM+1)*NT)
      DIMENSION WG((IM+2)*JM),W((JM+1)*IM)      

      CALL SJINIT(MM,NM,JM,IM,P,R,IT,T)

      DO L=1,(2*NN+1-MM)*MM+NN+1
        CALL SJL2NM(MM,L,N,M)
        S(L)=1D0/((ABS(M)+1D0)*(N+1D0))
      END DO

      CALL SJTSOG(MM,NM,NN,IM,JM,S,G,IT,T,P,Q,R,WS,WG,W,IPOW)

      CALL MPI_INIT(IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NP,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IPROC,IERR)

      CALL LJACHK(IA)
      JA=JM/(2**(IA+1))
      JH=JM/2
      
      JPH=((JA-1)/NP+1)*(2**IA)
      JS=JPH*IPROC+1

      CALL SJPINI(MM,NM,JM,JC,IM,PP,R,IT,T)
      CALL BSSET0((JM+1)*IM,W)
      CALL SJTS2G(MM,NM,NN,IM,JC,S,GD,IT,T,PP,Q,R,WS,WG,W,IPOW)

      EPS=1D-15
      INDEX=0
      DO J=1,JC/2
        DO I=1,IM
          DIFF=ABS(G(I,JM/2-(JS+J-1)+1)-GD(I,JC/2-J+1))
          BUNBO=ABS(G(I,JM/2-(JS+J-1)+1))+ABS(GD(I,JC/2-J+1))+1
          IF(DIFF/BUNBO.GT.EPS) THEN
            print *,DIFF/BUNBO
            INDEX=INDEX+1
          END IF
          DIFF=ABS(G(I,JM/2+(JS+J-1))-GD(I,JC/2+J))
          BUNBO=ABS(G(I,JM/2+(JS+J-1)))+ABS(GD(I,JC/2+J))+1
          IF(DIFF/BUNBO.GT.EPS) THEN
            print *,DIFF/BUNBO
            INDEX=INDEX+1
          END IF
        END DO
      END DO

      print *,'test:'
      IF(INDEX.EQ.0) THEN
        print *,'** OK'
      ELSE
        print *,'** Fail'
      END IF

      CALL MPI_FINALIZE(IERR)      

      END
