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
.globl ljnswg_
ljnswg_:
       movl (%rdi), %edi  # : JH が rdi に
       vbroadcastsd  (%rsi), %ymm0 # S1R を ymm0 の4箇所に
       vbroadcastsd 8(%rsi), %ymm1 # S1I を ymm1 の4箇所に
       vbroadcastsd  (%rdx), %ymm2 # S2R を ymm2 の4箇所に
       vbroadcastsd 8(%rdx), %ymm3 # S2I を ymm3 の4箇所に
       vbroadcastsd  (%rcx), %ymm4 # R を ymm4 の4箇所に	

       movq  8(%rsp), %r10  # : QB のベースアドレス       
       movq 16(%rsp), %r11  # : W1R のベースアドレス
       movq 24(%rsp), %rsi  # : W2R のベースアドレス                     
       
      # Y : r8, QA: r9, QB: r10, W1R: r11, W2R: rsi
      # W1I: rcx, W2I: rdi

       shlq $3,%rdi # JH*8 が rdi に
       movq $0,%rdx
       subq %rdi,%rdx
       
       movq %r11,%rcx
       addq %rdi,%rcx
       addq %rsi,%rdi

       subq %rdx,%r8	
       subq %rdx,%r9
       subq %rdx,%r10
       subq %rdx,%r11
       subq %rdx,%rcx
       subq %rdx,%rsi
       subq %rdx,%rdi	

.align 16
.L0:
       vmulpd (%r8,%rdx),%ymm4,%ymm5 # Y*R		
       vmovapd (%r9,%rdx), %ymm6 # QA
       vmulpd %ymm6,%ymm5,%ymm5 # R*Y*QA	
       vaddpd (%r10,%rdx), %ymm5,%ymm5 # 更新された QB が ymm5 に
       
       vmulpd %ymm0,%ymm6,%ymm7 # S1R*QA
       vaddpd (%r11,%rdx),%ymm7,%ymm7 # 更新された W1R が ymm7 に
       vmulpd %ymm1,%ymm6,%ymm8 # S1I*QA
       vaddpd (%rcx,%rdx),%ymm8,%ymm8 # 更新された W1I が ymm8 に
       vmulpd %ymm2,%ymm6,%ymm9 # S2R*QA
       vaddpd (%rsi,%rdx),%ymm9,%ymm9 # 更新された W2R が ymm9 に
       vmulpd %ymm3,%ymm6,%ymm10 # S2I*QA
       vaddpd (%rdi,%rdx),%ymm10,%ymm10 # 更新された W2I が ymm10 に
       
       vmovaps %ymm5,(%r10,%rdx) # 更新された QB をストア	
       vmovaps %ymm7,(%r11,%rdx)	
       vmovaps %ymm8,(%rcx,%rdx)	
       vmovaps %ymm9,(%rsi,%rdx)	
       vmovaps %ymm10,(%rdi,%rdx)

       addq $32,%rdx
       jnz .L0
       
       ret
       
