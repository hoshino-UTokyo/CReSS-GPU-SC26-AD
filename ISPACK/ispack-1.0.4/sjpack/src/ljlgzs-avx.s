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
.globl ljlgzs_
ljlgzs_:
       vzeroall	# 全 ymm を 0 にセット
       movl   (%rdi), %edi  # : JH が rdi に
       vbroadcastsd (%rdx), %ymm1 # R を ymm1 の4箇所に	
       movq  8(%rsp), %r10  # : W のベースアドレス
       
       shlq $3,%rdi # JH*8 が rdi に
       xorq %rdx,%rdx
       subq %rdi,%rdx
       
       addq %rdi,%rcx
       addq %rdi,%r8       
       addq %rdi,%r9
       addq %rdi,%r10
       
.align 16
.L0:
       vmulpd (%rcx,%rdx),%ymm1,%ymm3 # Y*R
       vmovapd (%r8,%rdx),%ymm2 # QA		
       vmulpd %ymm2,%ymm3,%ymm3 # R*Y*QA
       vaddpd (%r9,%rdx),%ymm3,%ymm3 # 更新された QB が ymm3 に
       vmovapd %ymm3,(%r9,%rdx)
       vmulpd (%r10,%rdx),%ymm2,%ymm2 # W*QA	
       vaddpd %ymm2,%ymm0,%ymm0 # S=S+W*QA
       addq $32,%rdx
       jnz .L0

       subq %rdi,%r10
       vmovapd (%r10),%ymm2 # W の先頭を一旦 ymm2 に退避
       vmovapd %ymm0,(%r10) # ymm0 を一旦 W の先頭に退避
       fldl (%r10)
       faddl 8(%r10)
       faddl 16(%r10)
       faddl 24(%r10)			
       fstpl (%rsi) # S の総和を rsi に
       vmovapd %ymm2,(%r10) # W の先頭を復元

       ret
       
