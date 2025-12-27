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
.globl ljngzs_
ljngzs_:
       vzeroall 
       vbroadcastsd  (%rcx), %ymm4 # R を ymm4 の4箇所に
       movl (%rdi), %ecx  # : JH が rcx に	

       shlq $3,%rcx # JH*8 が rcx に
       movq $0,%rdi
       subq %rcx,%rdi
       movq %rdi,%rax

       movq  8(%rsp), %r10  # : QB のベースアドレス       
       movq 16(%rsp), %r11  # : W1 のベースアドレス
       movq 24(%rsp), %rcx  # : W2 のベースアドレス                     
       
      # Y : r8, QA: r9, QB: r10, W1R: r11, W2R: rcx

       subq %rdi,%r8	
       subq %rdi,%r9
       subq %rdi,%r10
       subq %rdi,%r11
       subq %rdi,%rcx

.align 16
.L0:
       vmulpd (%r8,%rdi),%ymm4,%ymm5 # Y*R		
       vmovapd (%r9,%rdi), %ymm6 # QA
       vmulpd %ymm6,%ymm5,%ymm5 # R*Y*QA	
       vaddpd (%r10,%rdi), %ymm5,%ymm5 # 更新された QB が ymm5 に
       vmulpd (%r11,%rdi),%ymm6,%ymm7 # W1*QA が ymm7 に
       vmulpd (%rcx,%rdi),%ymm6,%ymm9 # W2*QA が ymm9 に
       vaddpd %ymm7,%ymm0,%ymm0
       vaddpd %ymm9,%ymm2,%ymm2
       vmovaps %ymm5,(%r10,%rdi) # 更新された QB をストア	
       addq $32,%rdi
       jnz .L0

       addq %rax,%r8		
       vmovapd (%r8),%ymm5 # Y の先頭を一旦 ymm5 に退避
       vmovapd %ymm0,(%r8) # ymm0 を Y の先頭に
       fldl (%r8)
       faddl 8(%r8)
       faddl 16(%r8)
       faddl 24(%r8)			
       fstpl (%rsi) # S1 の総和を (%rsi) に

       vmovapd %ymm2,(%r8) # ymm2 を Y の先頭に
       fldl (%r8)
       faddl 8(%r8)
       faddl 16(%r8)
       faddl 24(%r8)			
       fstpl (%rdx) # S2 の総和を (%rdx) に

       vmovapd %ymm5,(%r8) # Y の先頭を復元
       
       ret
       
