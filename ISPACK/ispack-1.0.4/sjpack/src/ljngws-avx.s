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
.globl ljngws_
ljngws_:
       pushq %r12
       vzeroall 
       vbroadcastsd  (%rcx), %ymm4 # R を ymm4 の4箇所に
       movl (%rdi), %ecx  # : JH が rcx に	

       shlq $3,%rcx # JH*8 が rcx に
       movq $0,%rdi
       subq %rcx,%rdi

       movq 16(%rsp), %r10  # : QB のベースアドレス       
       movq 24(%rsp), %r11  # : W1R のベースアドレス
       movq 32(%rsp), %r12  # : W2R のベースアドレス                     
       
      # Y : r8, QA: r9, QB: r10, W1R: r11, W2R: r12, W1I: rax, W2I: rcx

       subq %rdi,%r8	
       subq %rdi,%r9
       subq %rdi,%r10
       subq %rdi,%r11
       subq %rdi,%r12

       movq %r11,%rax
       movq %r12,%rcx
       subq %rdi,%rax
       subq %rdi,%rcx

.align 16
.L0:
       vmulpd (%r8,%rdi),%ymm4,%ymm5 # Y*R		
       vmovapd (%r9,%rdi), %ymm6 # QA
       vmulpd %ymm6,%ymm5,%ymm5 # R*Y*QA	
       vaddpd (%r10,%rdi), %ymm5,%ymm5 # 更新された QB が ymm5 に
       vmulpd (%r11,%rdi),%ymm6,%ymm7 # W1R*QA が ymm7 に
       vmulpd (%rax,%rdi),%ymm6,%ymm8 # W1I*QA が ymm8 に	
       vmulpd (%r12,%rdi),%ymm6,%ymm9 # W2R*QA が ymm9 に
       vmulpd (%rcx,%rdi),%ymm6,%ymm10 # W2I*QA が ymm10 に	
       vaddpd %ymm7,%ymm0,%ymm0
       vaddpd %ymm8,%ymm1,%ymm1
       vaddpd %ymm9,%ymm2,%ymm2
       vaddpd %ymm10,%ymm3,%ymm3
       vmovaps %ymm5,(%r10,%rdi) # 更新された QB をストア	
       addq $32,%rdi
       jnz .L0

       vmovapd (%r11),%ymm5 # W1I の先頭を一旦 ymm5 に退避
       vmovapd %ymm0,(%r11) # ymm0 を W1I の先頭に
       fldl (%r11)
       faddl 8(%r11)
       faddl 16(%r11)
       faddl 24(%r11)			
       fstpl (%rsi) # S1R の総和を (%rsi) に
       vmovapd %ymm1,(%r11) # ymm0 を W1I の先頭に
       fldl (%r11)
       faddl 8(%r11)
       faddl 16(%r11)
       faddl 24(%r11)			
       fstpl 8(%rsi) # S1I の総和を 8(%rsi) に
       vmovapd %ymm2,(%r11) # ymm2 を W1I の先頭に
       fldl (%r11)
       faddl 8(%r11)
       faddl 16(%r11)
       faddl 24(%r11)			
       fstpl (%rdx) # S2R の総和を (%rdx) に
       vmovapd %ymm3,(%r11) # ymm3 を W1I の先頭に
       fldl (%r11)
       faddl 8(%r11)
       faddl 16(%r11)
       faddl 24(%r11)			
       fstpl 8(%rdx) # S2I の総和を 8(%rdx) に
       vmovapd %ymm5,(%r11) # W1I の先頭を復元

       popq %r12
       
       ret
       
