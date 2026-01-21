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
************************************************************************
*     TRANSFORM GRID TO SPECTRA (using MPI)                   2011/09/08
*      (SJTGOSを利用する Hybrid 並列)
************************************************************************
      SUBROUTINE SJPGOS(MM,NM,NN,IM,JM,S,G,IT,T,P,Q,R,WS,WG,W,IPOW)

      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mpif.h'
      DIMENSION S((2*NN+1-MM)*MM+NN+1)
      DIMENSION WS((2*NN+1-MM)*MM+NN+1)
      DIMENSION P(JM/2,MM+4),Q(JM/2,7),R(*)
      DIMENSION W(*),G(0:IM-1,JM)
      DIMENSION WG(0:IM+1,JM)
      DIMENSION IT(2,2),T(IM*3,2)

      IF(JM.EQ.0) THEN
        CALL BSSET0((2*NN+1-MM)*MM+NN+1,S)
      ELSE
        CALL SJTGOS(MM,NM,NN,IM,JM,S,G,IT,T,P,Q,R,WS,WG,W,IPOW)
      END IF

      CALL MPI_ALLREDUCE(S,WS,(2*NN+1-MM)*MM+NN+1,
     &    MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

      CALL BSCOPY((2*NN+1-MM)*MM+NN+1,WS,S)

      END
