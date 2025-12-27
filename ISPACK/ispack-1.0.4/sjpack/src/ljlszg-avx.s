########################################################################
# ISPACK FORTRAN SUBROUTINE LIBRARY FOR SCIENTIFIC COMPUTING
# Copyright (C) 1998--2015 Keiichi Ishioka <ishioka@gfd-dennou.org>
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301 USA.
########################################################################
.text
.globl ljlszg_
ljlszg_:
       movl   (%rdi), %edi  # : JH が rdi に
       vbroadcastsd (%rsi), %ymm0 # S を ymm0 の4箇所に
       vbroadcastsd (%rdx), %ymm1 # R を ymm1 の4箇所に
       movq  8(%rsp), %r10  # : W のベースアドレス
       
       shlq $3,%rdi # JH*8 が rdi に
       
       movq $0,%rsi
       subq %rdi,%rsi
       
       addq %rdi,%rcx
       addq %rdi,%r8
       addq %rdi,%r9       
       addq %rdi,%r10
       
.align 16
.L0:
       vmulpd (%rcx,%rsi),%ymm1,%ymm4 # Y*R
       vmovapd (%r8,%rsi), %ymm2 # QA	
       vmulpd %ymm2,%ymm4,%ymm4 # R*Y*QA
       vmulpd %ymm0,%ymm2,%ymm2 # S*QA
       vaddpd (%r9,%rsi),%ymm4,%ymm4 # 更新された QB が ymm4 に
       vaddpd (%r10,%rsi),%ymm2,%ymm2 # 更新された W が ymm2 に
       vmovapd %ymm4,(%r9,%rsi)
       vmovapd %ymm2,(%r10,%rsi)	
       addq $32,%rsi
       jnz .L0
       
#------------------------------------

       ret
       
