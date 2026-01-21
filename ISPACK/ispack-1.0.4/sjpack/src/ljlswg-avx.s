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
.globl ljlswg_
ljlswg_:
       movl   (%rdi), %edi  # : JH が rdi に
       vbroadcastsd (%rsi), %ymm0 # SR を ymm0 の4箇所に
       vbroadcastsd 8(%rsi), %ymm5 # SI を ymm5 の4箇所に
       vbroadcastsd (%rdx), %ymm1 # R を ymm1 の4箇所に
       movq  8(%rsp), %r10  # : WR のベースアドレス
       
       shlq $3,%rdi # JH*8 が rdi に
       xorq %rdx,%rdx
       subq %rdi,%rdx
       
       movq %r10,%r11
       addq %rdi,%r11  # : WI のベースアドレス
       
       subq %rdx,%r8
       subq %rdx,%r9
       subq %rdx,%r10
       subq %rdx,%r11                     
       subq %rdx,%rcx
       
.align 16
.L0:
       vmulpd (%rcx,%rdx),%ymm1,%ymm4 # Y*R	
       vmovapd (%r8,%rdx), %ymm2 # QA
       vmulpd %ymm2,%ymm4,%ymm4 # R*Y*QA
       vmulpd %ymm0,%ymm2,%ymm3 # SR*QA
       vmulpd %ymm5,%ymm2,%ymm2 # SI*QA
       vaddpd (%r9,%rdx),%ymm4,%ymm4 # 更新された QB が ymm4 に
       vaddpd (%r10,%rdx),%ymm3,%ymm3 # 更新された WR が ymm3 に
       vaddpd (%r11,%rdx),%ymm2,%ymm2 # 更新された WI が ymm2 に	
       vmovapd %ymm4,(%r9,%rdx) # QB をセーブ
       vmovapd %ymm3,(%r10,%rdx) # WR をセーブ
       vmovapd %ymm2,(%r11,%rdx) # WI をセーブ
       addq $32,%rdx
       jnz .L0
       
       ret
       
