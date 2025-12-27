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
.globl ljnszg_
ljnszg_:
       movl (%rdi), %edi  # : JH が rdi に
       vbroadcastsd  (%rsi), %ymm0 # S1 を ymm0 の4箇所に
       vbroadcastsd  (%rdx), %ymm2 # S2 を ymm2 の4箇所に
       vbroadcastsd  (%rcx), %ymm4 # R を ymm4 の4箇所に	

       movq  8(%rsp), %r10  # : QB のベースアドレス       
       movq 16(%rsp), %r11  # : W1 のベースアドレス
       movq 24(%rsp), %rsi  # : W2 のベースアドレス                     
       
      # Y : r8, QA: r9, QB: r10, W1R: r11, W2R: rsi

       shlq $3,%rdi # JH*8 が rdi に
       movq $0,%rdx
       subq %rdi,%rdx
       
       subq %rdx,%r8	
       subq %rdx,%r9
       subq %rdx,%r10
       subq %rdx,%r11
       subq %rdx,%rsi

.align 16
.L0:
       vmulpd (%r8,%rdx),%ymm4,%ymm5 # Y*R		
       vmovapd (%r9,%rdx), %ymm6 # QA
       vmulpd %ymm6,%ymm5,%ymm5 # R*Y*QA	
       vaddpd (%r10,%rdx), %ymm5,%ymm5 # 更新された QB が ymm5 に
       
       vmulpd %ymm0,%ymm6,%ymm7 # S1*QA
       vaddpd (%r11,%rdx),%ymm7,%ymm7 # 更新された W1 が ymm7 に
       vmulpd %ymm2,%ymm6,%ymm9 # S2*QA
       vaddpd (%rsi,%rdx),%ymm9,%ymm9 # 更新された W2 が ymm9 に
       
       vmovaps %ymm5,(%r10,%rdx) # 更新された QB をストア	
       vmovaps %ymm7,(%r11,%rdx)	
       vmovaps %ymm9,(%rsi,%rdx)	

       addq $32,%rdx
       jnz .L0
       
       ret
       
