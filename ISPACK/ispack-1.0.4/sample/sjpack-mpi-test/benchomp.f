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
      PARAMETER(MM=170,JM=256,IM=512)            
      PARAMETER(NM=MM,NN=MM,IPOW=0)
      PARAMETER(NTRY=8,NTR=100)
      DIMENSION S((2*NN+1-MM)*MM+NN+1),G(IM,JM)
      DIMENSION SD((2*NN+1-MM)*MM+NN+1)
      DIMENSION IT(4),T(IM*6)      
      PARAMETER(NT=32)
      DIMENSION P(JM/2*(MM+4)),Q(JM/2*7*NT),R((MM+1)*(2*NM-MM-1)+1)
      DIMENSION WS((2*NN+1-MM)*MM+NN+1+2*(NM+1)*NT)
      DIMENSION WG((IM+2)*JM),W((JM+1)*IM)      
!$    INTEGER omp_get_max_threads

      RC=1D0*5*IM*LOG(1D0*IM)/LOG(2D0)*0.5D0*JM+1D0*(MM+1)*(MM+1)*JM
      ! 逆変換/正変換の1回あたりの演算数概算

      DO L=1,(2*NN+1-MM)*MM+NN+1
        CALL SJL2NM(MM,L,N,M)
        S(L)=1D0/((ABS(M)+1D0)*(N+1D0))
      END DO

      CALL MPI_INIT(IERR)

      CALL SJPINI(MM,NM,JM,JC,IM,P,R,IT,T)
      DO ITRY=1,NTRY
        CALL APTIME(TIM0)
        DO ITR=1,NTR      
          CALL SJTSOG(MM,NM,NN,IM,JC,S,G,IT,T,P,Q,R,WS,WG,W,IPOW)
          CALL SJPGOS(MM,NM,NN,IM,JC,SD,G,IT,T,P,Q,R,WS,WG,W,IPOW)
        END DO
        CALL APTIME(TIM1)
        TIM=TIM1-TIM0
        IF(ITRY.EQ.1) THEN
          TIMMIN=TIM
        ELSE IF(TIM.LT.TIMMIN) THEN
          TIMMIN=TIM
        END IF
      END DO

!$      MAXTD=omp_get_max_threads()
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IPROC,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NP,IERR)
      IF(IPROC.EQ.0) THEN
        GFLOPS=2*RC*NTR/TIMMIN/1D9
        WRITE(6,'(A,I3,A,I3,A,F6.3,A)') 
     &      'SJPACK-MPI(',NP,' processes, ',MAXTD,' threads): ',
     &      GFLOPS,' GFlops'
      END IF

      CALL MPI_FINALIZE(IERR)      

      END
