	.file	"FatropLinearAlgebraBlasfeo.cpp"
	.text
	.globl	_ZN6fatrop6max_elEiiP12blasfeo_dmatii
	.type	_ZN6fatrop6max_elEiiP12blasfeo_dmatii, @function
_ZN6fatrop6max_elEiiP12blasfeo_dmatii:
.LFB2018:
	.cfi_startproc
	endbr64
	pushq	%r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	pushq	%r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	pushq	%r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	pushq	%rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	pushq	%rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	movl	%ecx, %r12d
	cmpl	%esi, %r8d
	jge	.L8
	movl	%esi, %ebx
	movq	%rdx, %rbp
	movl	%r8d, %r13d
	movl	%ecx, %esi
	pxor	%xmm1, %xmm1
	movq	.LC1(%rip), %xmm2
	jmp	.L7
.L4:
	addl	$1, %r9d
	cmpl	%r9d, %edi
	je	.L3
.L6:
	movl	%r9d, %eax
	andl	$-4, %eax
	imull	%r11d, %eax
	addl	%r10d, %eax
	movl	%r9d, %r14d
	andl	$3, %r14d
	addl	%r14d, %eax
	cltq
	movsd	(%rdx,%rax,8), %xmm0
	andpd	%xmm2, %xmm0
	comisd	%xmm1, %xmm0
	jb	.L4
	movl	%ecx, %r13d
	movl	%r9d, %esi
	movapd	%xmm0, %xmm1
	jmp	.L4
.L3:
	addl	$1, %r8d
	cmpl	%r8d, %ebx
	je	.L2
.L7:
	cmpl	%edi, %r12d
	jge	.L3
	movq	8(%rbp), %rdx
	movl	36(%rbp), %r11d
	leal	0(,%r8,4), %r10d
	movl	%r12d, %r9d
	movl	%r8d, %ecx
	jmp	.L6
.L8:
	movl	%r8d, %r13d
	movl	%ecx, %esi
.L2:
	salq	$32, %r13
	movl	%esi, %eax
	orq	%r13, %rax
	popq	%rbx
	.cfi_def_cfa_offset 40
	popq	%rbp
	.cfi_def_cfa_offset 32
	popq	%r12
	.cfi_def_cfa_offset 24
	popq	%r13
	.cfi_def_cfa_offset 16
	popq	%r14
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE2018:
	.size	_ZN6fatrop6max_elEiiP12blasfeo_dmatii, .-_ZN6fatrop6max_elEiiP12blasfeo_dmatii
	.globl	_ZN6fatrop7LU_FACTEiiiRiP12blasfeo_dmatPKiS4_
	.type	_ZN6fatrop7LU_FACTEiiiRiP12blasfeo_dmatPKiS4_, @function
_ZN6fatrop7LU_FACTEiiiRiP12blasfeo_dmatPKiS4_:
.LFB2019:
	.cfi_startproc
	endbr64
	pushq	%r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	pushq	%r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	pushq	%r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	pushq	%rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	pushq	%rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	cmpl	%edx, %edi
	movl	%edx, %r14d
	cmovle	%edi, %r14d
	testl	%r14d, %r14d
	jle	.L13
	movl	%edi, %ebp
	movl	%edx, %r12d
	movq	%r8, %r13
	movl	$0, %ebx
.L15:
	movl	%ebx, %r8d
	movl	%ebx, %ecx
	movq	%r13, %rdx
	movl	%r12d, %esi
	movl	%ebp, %edi
	call	_ZN6fatrop6max_elEiiP12blasfeo_dmatii
	addl	$1, %ebx
	cmpl	%ebx, %r14d
	jne	.L15
.L13:
	popq	%rbx
	.cfi_def_cfa_offset 40
	popq	%rbp
	.cfi_def_cfa_offset 32
	popq	%r12
	.cfi_def_cfa_offset 24
	popq	%r13
	.cfi_def_cfa_offset 16
	popq	%r14
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE2019:
	.size	_ZN6fatrop7LU_FACTEiiiRiP12blasfeo_dmatPKiS4_, .-_ZN6fatrop7LU_FACTEiiiRiP12blasfeo_dmatPKiS4_
	.type	_GLOBAL__sub_I__ZN6fatrop6max_elEiiP12blasfeo_dmatii, @function
_GLOBAL__sub_I__ZN6fatrop6max_elEiiP12blasfeo_dmatii:
.LFB2576:
	.cfi_startproc
	endbr64
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	leaq	_ZStL8__ioinit(%rip), %rdi
	call	_ZNSt8ios_base4InitC1Ev@PLT
	leaq	__dso_handle(%rip), %rdx
	leaq	_ZStL8__ioinit(%rip), %rsi
	movq	_ZNSt8ios_base4InitD1Ev@GOTPCREL(%rip), %rdi
	call	__cxa_atexit@PLT
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE2576:
	.size	_GLOBAL__sub_I__ZN6fatrop6max_elEiiP12blasfeo_dmatii, .-_GLOBAL__sub_I__ZN6fatrop6max_elEiiP12blasfeo_dmatii
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I__ZN6fatrop6max_elEiiP12blasfeo_dmatii
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC1:
	.long	4294967295
	.long	2147483647
	.long	0
	.long	0
	.hidden	__dso_handle
	.ident	"GCC: (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	 1f - 0f
	.long	 4f - 1f
	.long	 5
0:
	.string	 "GNU"
1:
	.align 8
	.long	 0xc0000002
	.long	 3f - 2f
2:
	.long	 0x3
3:
	.align 8
4:
