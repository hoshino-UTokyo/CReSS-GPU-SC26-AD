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
*     INITIALIZATION OF SJPACK-MPI                            2011/08/25
*-----------------------------------------------------------------------
      SUBROUTINE SJPINI(MM,NM,JM,JC,IM,P,R,IT,T)

      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mpif.h'      
      DIMENSION P(JM/2,MM+4),R((MM+1)*(2*NM-MM-1)+1),IT(2,2),T(IM*3,2)

      CALL LJACHK(IA)

      IF(IA.EQ.0) THEN
        IF(MOD(JM,2).NE.0) CALL BSDMSG('E','SJPINT','JM must be even.')
      ELSE IF(IA.EQ.1) THEN
        IF(MOD(JM,4).NE.0) CALL BSDMSG('E','SJPINT',
     &      'JM must be a multiple of 4 to use SSE.')
      ELSE IF(IA.EQ.2) THEN
        IF(MOD(JM,8).NE.0) CALL BSDMSG('E','SJPINT',
     &      'JM must be a multiple of 8 to use AVX.')
      END IF

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IPROC,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NP,IERR)

      JA=JM/(2**(IA+1))
      JH=JM/2
      
      JPH=((JA-1)/NP+1)*(2**IA)
      JS=JPH*IPROC+1
      JE=MIN(JPH*(IPROC+1),JH)
      IF(JE.GE.JS) THEN
        JCH=JE-JS+1
        JC=JCH*2
      ELSE
        JC=0
        JS=1
        JE=1
      END IF

      IF(JC.GT.0) THEN
        CALL LJPINI(MM,NM,JM,JC,JS,P,R)
        CALL FJRINI(IM,1,1,IT(1,1),T(1,1))
        CALL FJRINI(IM,2,2,IT(1,2),T(1,2))
      END IF

      END
