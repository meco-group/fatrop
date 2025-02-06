__asm__(" .text\n\t"
" .p2align 4,,15\n\t"
" .def inner_kernel_dgemv_add_n_4_lib4; .scl 2; .type 32; .endef; inner_kernel_dgemv_add_n_4_lib4:\n\t"
"\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jle 2f\n\t"
"\n\t"
" cmpl $ 4, %r10d\n\t"
" jl 0f\n\t"
"\n\t"
"\n\t"
" .p2align 3\n\t"
"1:\n\t"
"\n\t"
" movddup 0(%r12), %xmm12\n\t"
" movapd 0(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm0\n\t"
" movapd 16(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm1\n\t"
" subl $ 4, %r10d\n\t"
"\n\t"
" movddup 8(%r12), %xmm12\n\t"
" movapd 32(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm2\n\t"
" movapd 48(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm3\n\t"
"\n\t"
" movddup 16(%r12), %xmm12\n\t"
" movapd 64(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm0\n\t"
" movapd 80(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm1\n\t"
"\n\t"
" movddup 24(%r12), %xmm12\n\t"
" movapd 96(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm2\n\t"
" movapd 112(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm3\n\t"
"\n\t"
" addq $ 128, %r11\n\t"
" addq $ 32, %r12\n\t"
"\n\t"
" cmpl $ 3, %r10d\n\t"
"\n\t"
" jg 1b\n\t"
"\n\t"
"\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jle 2f\n\t"
"\n\t"
"0:\n\t"
"\n\t"
" movddup 0(%r12), %xmm12\n\t"
" movapd 0(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm0\n\t"
" movapd 16(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm1\n\t"
"\n\t"
" addq $ 32, %r11\n\t"
" addq $ 8, %r12\n\t"
"\n\t"
" subl $ 1, %r10d\n\t"
" cmpl $ 0, %r10d\n\t"
"\n\t"
" jg 0b\n\t"
"\n\t"
"2:\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .def inner_kernel_dgemv_add_t_4_lib4; .scl 2; .type 32; .endef; inner_kernel_dgemv_add_t_4_lib4:\n\t"
"\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jle 2f\n\t"
"\n\t"
" cmpl $ 4, %r10d\n\t"
" jl 0f\n\t"
"\n\t"
"\n\t"
" .p2align 3\n\t"
"1:\n\t"
"\n\t"
" movupd 0(%r13), %xmm12\n\t"
"\n\t"
" movapd 0(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm0\n\t"
" subl $ 4, %r10d\n\t"
"\n\t"
" movapd 32(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm1\n\t"
"\n\t"
" movapd 64(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm2\n\t"
"\n\t"
" movapd 96(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm3\n\t"
"\n\t"
" movupd 16(%r13), %xmm12\n\t"
"\n\t"
" movapd 16(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm0\n\t"
"\n\t"
" movapd 48(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm1\n\t"
"\n\t"
" movapd 80(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm2\n\t"
"\n\t"
" movapd 112(%r11), %xmm8\n\t"
" mulpd %xmm12, %xmm8\n\t"
" addpd %xmm8, %xmm3\n\t"
"\n\t"
" addq %r12, %r11\n\t"
" addq $ 32, %r13\n\t"
"\n\t"
" cmpl $ 3, %r10d\n\t"
" jg 1b\n\t"
"\n\t"
"\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jle 2f\n\t"
"\n\t"
"0:\n\t"
"\n\t"
" movsd 0(%r13), %xmm12\n\t"
"\n\t"
" movsd 0(%r11), %xmm8\n\t"
" mulsd %xmm12, %xmm8\n\t"
" addsd %xmm8, %xmm0\n\t"
" subl $ 1, %r10d\n\t"
"\n\t"
" movsd 32(%r11), %xmm8\n\t"
" mulsd %xmm12, %xmm8\n\t"
" addsd %xmm8, %xmm1\n\t"
"\n\t"
" movsd 64(%r11), %xmm8\n\t"
" mulsd %xmm12, %xmm8\n\t"
" addsd %xmm8, %xmm2\n\t"
"\n\t"
" movsd 96(%r11), %xmm8\n\t"
" mulsd %xmm12, %xmm8\n\t"
" addsd %xmm8, %xmm3\n\t"
"\n\t"
" addq $ 8, %r11\n\t"
" addq $ 8, %r13\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jg 0b\n\t"
"\n\t"
"\n\t"
"2:\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .def inner_kernel_dgemv_add_nt_4_lib4; .scl 2; .type 32; .endef; inner_kernel_dgemv_add_nt_4_lib4:\n\t"
"\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jle 2f\n\t"
"\n\t"
" cmpl $ 4, %r10d\n\t"
" jl 0f\n\t"
"\n\t"
"\n\t"
" .p2align 3\n\t"
"1:\n\t"
"\n\t"
" movupd 0(%r13), %xmm9\n\t"
" movupd 16(%r13), %xmm10\n\t"
" movupd 0(%r14), %xmm11\n\t"
" movupd 16(%r14), %xmm12\n\t"
"\n\t"
" subl $ 4, %r10d\n\t"
"\n\t"
" movapd 0(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm9, %xmm14\n\t"
" addpd %xmm14, %xmm0\n\t"
" mulpd %xmm4, %xmm15\n\t"
" addpd %xmm15, %xmm11\n\t"
"\n\t"
" movapd 16(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm10, %xmm14\n\t"
" addpd %xmm14, %xmm0\n\t"
" mulpd %xmm4, %xmm15\n\t"
" addpd %xmm15, %xmm12\n\t"
"\n\t"
" movapd 32(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm9, %xmm14\n\t"
" addpd %xmm14, %xmm1\n\t"
" mulpd %xmm5, %xmm15\n\t"
" addpd %xmm15, %xmm11\n\t"
"\n\t"
" movapd 48(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm10, %xmm14\n\t"
" addpd %xmm14, %xmm1\n\t"
" mulpd %xmm5, %xmm15\n\t"
" addpd %xmm15, %xmm12\n\t"
"\n\t"
" movapd 64(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm9, %xmm14\n\t"
" addpd %xmm14, %xmm2\n\t"
" mulpd %xmm6, %xmm15\n\t"
" addpd %xmm15, %xmm11\n\t"
"\n\t"
" movapd 80(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm10, %xmm14\n\t"
" addpd %xmm14, %xmm2\n\t"
" mulpd %xmm6, %xmm15\n\t"
" addpd %xmm15, %xmm12\n\t"
"\n\t"
" movapd 96(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm9, %xmm14\n\t"
" addpd %xmm14, %xmm3\n\t"
" mulpd %xmm7, %xmm15\n\t"
" addpd %xmm15, %xmm11\n\t"
"\n\t"
" movapd 112(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm10, %xmm14\n\t"
" addpd %xmm14, %xmm3\n\t"
" mulpd %xmm7, %xmm15\n\t"
" addpd %xmm15, %xmm12\n\t"
"\n\t"
" movupd %xmm11, 0(%r14)\n\t"
" movupd %xmm12, 16(%r14)\n\t"
"\n\t"
" addq %r12, %r11\n\t"
" addq $ 32, %r13\n\t"
" addq $ 32, %r14\n\t"
"\n\t"
" cmpl $ 3, %r10d\n\t"
" jg 1b\n\t"
"\n\t"
"\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jle 2f\n\t"
"\n\t"
"0:\n\t"
"\n\t"
" movsd 0(%r13), %xmm9\n\t"
" movsd 0(%r14), %xmm11\n\t"
"\n\t"
" subl $ 1, %r10d\n\t"
"\n\t"
" movsd 0(%r11), %xmm14\n\t"
" movsd %xmm14, %xmm15\n\t"
" mulsd %xmm9, %xmm14\n\t"
" addsd %xmm14, %xmm0\n\t"
" mulsd %xmm4, %xmm15\n\t"
" addsd %xmm15, %xmm11\n\t"
"\n\t"
" movsd 32(%r11), %xmm14\n\t"
" movsd %xmm14, %xmm15\n\t"
" mulsd %xmm9, %xmm14\n\t"
" addsd %xmm14, %xmm1\n\t"
" mulsd %xmm5, %xmm15\n\t"
" addsd %xmm15, %xmm11\n\t"
"\n\t"
" movsd 64(%r11), %xmm14\n\t"
" movsd %xmm14, %xmm15\n\t"
" mulsd %xmm9, %xmm14\n\t"
" addsd %xmm14, %xmm2\n\t"
" mulsd %xmm6, %xmm15\n\t"
" addsd %xmm15, %xmm11\n\t"
"\n\t"
" movsd 96(%r11), %xmm14\n\t"
" movsd %xmm14, %xmm15\n\t"
" mulsd %xmm9, %xmm14\n\t"
" addsd %xmm14, %xmm3\n\t"
" mulsd %xmm7, %xmm15\n\t"
" addsd %xmm15, %xmm11\n\t"
"\n\t"
" movsd %xmm11, 0(%r14)\n\t"
"\n\t"
" addq $ 8, %r11\n\t"
" addq $ 8, %r13\n\t"
" addq $ 8, %r14\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jg 0b\n\t"
"\n\t"
"2:\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .def inner_edge_dgemv_add_t_4_lib4; .scl 2; .type 32; .endef; inner_edge_dgemv_add_t_4_lib4:\n\t"
"\n\t"
"\n\t"
" cmpl $ 0, %r14d\n\t"
" jle 2f\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jle 2f\n\t"
"\n\t"
" movl $ 4, %r15d\n\t"
" subl %r14d, %r15d\n\t"
" cmpl %r10d, %r15d\n\t"
" cmovgl %r10d, %r15d\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"1:\n\t"
" movsd 0(%r13), %xmm12\n\t"
"\n\t"
" movsd 0(%r11), %xmm8\n\t"
" mulsd %xmm12, %xmm8\n\t"
" addsd %xmm8, %xmm0\n\t"
" subl $ 1, %r10d\n\t"
"\n\t"
" movsd 32(%r11), %xmm8\n\t"
" mulsd %xmm12, %xmm8\n\t"
" addsd %xmm8, %xmm1\n\t"
"\n\t"
" movsd 64(%r11), %xmm8\n\t"
" mulsd %xmm12, %xmm8\n\t"
" addsd %xmm8, %xmm2\n\t"
"\n\t"
" movsd 96(%r11), %xmm8\n\t"
" mulsd %xmm12, %xmm8\n\t"
" addsd %xmm8, %xmm3\n\t"
"\n\t"
" subl $ 1, %r10d\n\t"
" subl $ 1, %r15d\n\t"
" addq $ 8, %r11\n\t"
" addq $ 8, %r13\n\t"
"\n\t"
" cmpl $ 0, %r15d\n\t"
" jg 1b\n\t"
"\n\t"
" cmpl $ 0, %r10d\n\t"
" jle 2f\n\t"
"\n\t"
" addq %r12, %r11\n\t"
" subq $ 32, %r11\n\t"
"\n\t"
"2:\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .def inner_edge_dsymv_add_nt_4_lib4; .scl 2; .type 32; .endef; inner_edge_dsymv_add_nt_4_lib4:\n\t"
"\n\t"
"\n\t"
" xorpd %xmm13, %xmm13\n\t"
"\n\t"
" movupd 0(%r13), %xmm9\n\t"
" movupd 16(%r13), %xmm10\n\t"
" movupd 0(%r14), %xmm11\n\t"
" movupd 16(%r14), %xmm12\n\t"
"\n\t"
"\n\t"
" movapd 0(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm9, %xmm14\n\t"
" addpd %xmm14, %xmm0\n\t"
" movsd %xmm13, %xmm15\n\t"
" mulpd %xmm4, %xmm15\n\t"
" addpd %xmm15, %xmm11\n\t"
"\n\t"
" movapd 16(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm10, %xmm14\n\t"
" addpd %xmm14, %xmm0\n\t"
" mulpd %xmm4, %xmm15\n\t"
" addpd %xmm15, %xmm12\n\t"
"\n\t"
"\n\t"
" movapd 32(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" movsd %xmm13, %xmm14\n\t"
" mulpd %xmm9, %xmm14\n\t"
" addpd %xmm14, %xmm1\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movapd 48(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm10, %xmm14\n\t"
" addpd %xmm14, %xmm1\n\t"
" mulpd %xmm5, %xmm15\n\t"
" addpd %xmm15, %xmm12\n\t"
" movapd 80(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" mulpd %xmm10, %xmm14\n\t"
" addpd %xmm14, %xmm2\n\t"
" movsd %xmm13, %xmm15\n\t"
" mulpd %xmm6, %xmm15\n\t"
" addpd %xmm15, %xmm12\n\t"
" movapd 112(%r11), %xmm14\n\t"
" movapd %xmm14, %xmm15\n\t"
" movsd %xmm13, %xmm14\n\t"
" mulpd %xmm10, %xmm14\n\t"
" addpd %xmm14, %xmm3\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movupd %xmm11, 0(%r14)\n\t"
" movupd %xmm12, 16(%r14)\n\t"
"\n\t"
" addq %r12, %r11\n\t"
" addq $ 32, %r13\n\t"
" addq $ 32, %r14\n\t"
"\n\t"
" subq $ 4, %r10\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .macro INNER_BLEND_N_SCALE_AB_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" addpd %xmm2, %xmm0\n\t"
" addpd %xmm3, %xmm1\n\t"
"\n\t"
"\n\t"
" movddup 0(%r10), %xmm15\n\t"
" mulpd %xmm15, %xmm0\n\t"
" mulpd %xmm15, %xmm1\n\t"
"\n\t"
"\n\t"
" movddup 0(%r11), %xmm15\n\t"
"\n\t"
" xorpd %xmm14, %xmm14\n\t"
" ucomisd %xmm14, %xmm15\n\t"
" je 0f\n\t"
"\n\t"
" movupd 0(%r12), %xmm14\n\t"
" mulpd %xmm15, %xmm14\n\t"
" addpd %xmm14, %xmm0\n\t"
" movupd 16(%r12), %xmm14\n\t"
" mulpd %xmm15, %xmm14\n\t"
" addpd %xmm14, %xmm1\n\t"
"\n\t"
"0:\n\t"
"\n\t"
"\n\t"
" .endm\n\t"
" .macro INNER_BLEND_T_SCALE_AB_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" haddpd %xmm1, %xmm0\n\t"
" haddpd %xmm3, %xmm2\n\t"
" movapd %xmm2, %xmm1\n\t"
"\n\t"
"\n\t"
" movddup 0(%r10), %xmm15\n\t"
" mulpd %xmm15, %xmm0\n\t"
" mulpd %xmm15, %xmm1\n\t"
"\n\t"
"\n\t"
" movddup 0(%r11), %xmm15\n\t"
"\n\t"
" xorpd %xmm14, %xmm14\n\t"
" ucomisd %xmm14, %xmm15\n\t"
" je 0f\n\t"
"\n\t"
" movupd 0(%r12), %xmm14\n\t"
" mulpd %xmm15, %xmm14\n\t"
" addpd %xmm14, %xmm0\n\t"
" movupd 16(%r12), %xmm14\n\t"
" mulpd %xmm15, %xmm14\n\t"
" addpd %xmm14, %xmm1\n\t"
"\n\t"
"0:\n\t"
"\n\t"
"\n\t"
" .endm\n\t"
" .macro INNER_BLEND_T_SCALE_A1_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" haddpd %xmm1, %xmm0\n\t"
" haddpd %xmm3, %xmm2\n\t"
" movapd %xmm2, %xmm1\n\t"
"\n\t"
"\n\t"
" movddup 0(%r10), %xmm15\n\t"
" mulpd %xmm15, %xmm0\n\t"
" mulpd %xmm15, %xmm1\n\t"
"\n\t"
"\n\t"
" movupd 0(%r11), %xmm14\n\t"
" addpd %xmm14, %xmm0\n\t"
" movupd 16(%r11), %xmm14\n\t"
" addpd %xmm14, %xmm1\n\t"
"\n\t"
"\n\t"
" .endm\n\t"
" .macro INNER_STORE_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movupd %xmm0, 0(%r10)\n\t"
" movupd %xmm1, 16(%r10)\n\t"
"\n\t"
"\n\t"
" .endm\n\t"
" .macro INNER_STORE_4_VS_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" cmpl $ 0, %r11d\n\t"
" jle 0f\n\t"
"\n\t"
" movsd %xmm0, 0(%r10)\n\t"
"\n\t"
" cmpl $ 1, %r11d\n\t"
" jle 0f\n\t"
"\n\t"
" movhpd %xmm0, 8(%r10)\n\t"
"\n\t"
" cmpl $ 2, %r11d\n\t"
" jle 0f\n\t"
"\n\t"
" movsd %xmm1, 16(%r10)\n\t"
"\n\t"
" cmpl $ 3, %r11d\n\t"
" jle 0f\n\t"
"\n\t"
" movhpd %xmm1, 24(%r10)\n\t"
"\n\t"
"0:\n\t"
"\n\t"
"\n\t"
" .endm\n\t"
" .p2align 4,,15\n\t"
" .globl kernel_dgemv_n_4_lib4; .def kernel_dgemv_n_4_lib4; .scl 2; .type 32; .endef; kernel_dgemv_n_4_lib4:\n\t"
"\n\t"
" subq $256, %rsp; movq %rbx, (%rsp); movq %rbp, 8(%rsp); movq %r12, 16(%rsp); movq %r13, 24(%rsp); movq %r14, 32(%rsp); movq %r15, 40(%rsp); movq %rdi, 48(%rsp); movq %rsi, 56(%rsp); movups %xmm6, 64(%rsp); movups %xmm7, 80(%rsp); movups %xmm8, 96(%rsp); movups %xmm9, 112(%rsp); movups %xmm10, 128(%rsp); movups %xmm11, 144(%rsp); movups %xmm12, 160(%rsp); movups %xmm13, 176(%rsp); movups %xmm14, 192(%rsp); movups %xmm15, 208(%rsp);\n\t"
"\n\t"
"\n\t"
"\n\t"
" xorpd %xmm0, %xmm0; movapd %xmm0, %xmm1; movapd %xmm0, %xmm2; movapd %xmm0, %xmm3\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rcx, %r10\n\t"
" movq %r8, %r11\n\t"
" movq %r9, %r12\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_kernel_dgemv_add_n_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rdx, %r10\n\t"
" movq 256 + 40(%rsp), %r11\n\t"
" movq 256 + 48(%rsp), %r12\n\t"
"\n\t"
"\n\t"
" INNER_BLEND_N_SCALE_AB_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq 256 + 56(%rsp), %r10\n\t"
"\n\t"
"\n\t"
" INNER_STORE_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq (%rsp), %rbx; movq 8(%rsp), %rbp; movq 16(%rsp), %r12; movq 24(%rsp), %r13; movq 32(%rsp), %r14; movq 40(%rsp), %r15; movq 48(%rsp), %rdi; movq 56(%rsp), %rsi; movups 64(%rsp), %xmm6; movups 80(%rsp), %xmm7; movups 96(%rsp), %xmm8; movups 112(%rsp), %xmm9; movups 128(%rsp), %xmm10; movups 144(%rsp), %xmm11; movups 160(%rsp), %xmm12; movups 176(%rsp), %xmm13; movups 192(%rsp), %xmm14; movups 208(%rsp), %xmm15; addq $256, %rsp;\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .globl kernel_dgemv_n_4_vs_lib4; .def kernel_dgemv_n_4_vs_lib4; .scl 2; .type 32; .endef; kernel_dgemv_n_4_vs_lib4:\n\t"
"\n\t"
" subq $256, %rsp; movq %rbx, (%rsp); movq %rbp, 8(%rsp); movq %r12, 16(%rsp); movq %r13, 24(%rsp); movq %r14, 32(%rsp); movq %r15, 40(%rsp); movq %rdi, 48(%rsp); movq %rsi, 56(%rsp); movups %xmm6, 64(%rsp); movups %xmm7, 80(%rsp); movups %xmm8, 96(%rsp); movups %xmm9, 112(%rsp); movups %xmm10, 128(%rsp); movups %xmm11, 144(%rsp); movups %xmm12, 160(%rsp); movups %xmm13, 176(%rsp); movups %xmm14, 192(%rsp); movups %xmm15, 208(%rsp);\n\t"
"\n\t"
"\n\t"
"\n\t"
" xorpd %xmm0, %xmm0; movapd %xmm0, %xmm1; movapd %xmm0, %xmm2; movapd %xmm0, %xmm3\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rcx, %r10\n\t"
" movq %r8, %r11\n\t"
" movq %r9, %r12\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_kernel_dgemv_add_n_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rdx, %r10\n\t"
" movq 256 + 40(%rsp), %r11\n\t"
" movq 256 + 48(%rsp), %r12\n\t"
"\n\t"
"\n\t"
" INNER_BLEND_N_SCALE_AB_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq 256 + 56(%rsp), %r10\n\t"
" movq 256 + 64(%rsp), %r11\n\t"
"\n\t"
"\n\t"
" INNER_STORE_4_VS_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq (%rsp), %rbx; movq 8(%rsp), %rbp; movq 16(%rsp), %r12; movq 24(%rsp), %r13; movq 32(%rsp), %r14; movq 40(%rsp), %r15; movq 48(%rsp), %rdi; movq 56(%rsp), %rsi; movups 64(%rsp), %xmm6; movups 80(%rsp), %xmm7; movups 96(%rsp), %xmm8; movups 112(%rsp), %xmm9; movups 128(%rsp), %xmm10; movups 144(%rsp), %xmm11; movups 160(%rsp), %xmm12; movups 176(%rsp), %xmm13; movups 192(%rsp), %xmm14; movups 208(%rsp), %xmm15; addq $256, %rsp;\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .globl kernel_dgemv_t_4_lib4; .def kernel_dgemv_t_4_lib4; .scl 2; .type 32; .endef; kernel_dgemv_t_4_lib4:\n\t"
"\n\t"
" subq $256, %rsp; movq %rbx, (%rsp); movq %rbp, 8(%rsp); movq %r12, 16(%rsp); movq %r13, 24(%rsp); movq %r14, 32(%rsp); movq %r15, 40(%rsp); movq %rdi, 48(%rsp); movq %rsi, 56(%rsp); movups %xmm6, 64(%rsp); movups %xmm7, 80(%rsp); movups %xmm8, 96(%rsp); movups %xmm9, 112(%rsp); movups %xmm10, 128(%rsp); movups %xmm11, 144(%rsp); movups %xmm12, 160(%rsp); movups %xmm13, 176(%rsp); movups %xmm14, 192(%rsp); movups %xmm15, 208(%rsp);\n\t"
"\n\t"
"\n\t"
"\n\t"
" xorpd %xmm0, %xmm0; movapd %xmm0, %xmm1; movapd %xmm0, %xmm2; movapd %xmm0, %xmm3\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rcx, %r10\n\t"
" movq %r9, %r11\n\t"
" movq 256 + 40(%rsp), %r12\n\t"
" sall $ 5, %r12d\n\t"
"\n\t"
" movq 256 + 48(%rsp), %r13\n\t"
" movq %r8, %r14\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_edge_dgemv_add_t_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_kernel_dgemv_add_t_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rdx, %r10\n\t"
" movq 256 + 56(%rsp), %r11\n\t"
" movq 256 + 64(%rsp), %r12\n\t"
"\n\t"
"\n\t"
" INNER_BLEND_T_SCALE_AB_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq 256 + 72(%rsp), %r10\n\t"
"\n\t"
"\n\t"
" INNER_STORE_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq (%rsp), %rbx; movq 8(%rsp), %rbp; movq 16(%rsp), %r12; movq 24(%rsp), %r13; movq 32(%rsp), %r14; movq 40(%rsp), %r15; movq 48(%rsp), %rdi; movq 56(%rsp), %rsi; movups 64(%rsp), %xmm6; movups 80(%rsp), %xmm7; movups 96(%rsp), %xmm8; movups 112(%rsp), %xmm9; movups 128(%rsp), %xmm10; movups 144(%rsp), %xmm11; movups 160(%rsp), %xmm12; movups 176(%rsp), %xmm13; movups 192(%rsp), %xmm14; movups 208(%rsp), %xmm15; addq $256, %rsp;\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .globl kernel_dgemv_t_4_vs_lib4; .def kernel_dgemv_t_4_vs_lib4; .scl 2; .type 32; .endef; kernel_dgemv_t_4_vs_lib4:\n\t"
"\n\t"
" subq $256, %rsp; movq %rbx, (%rsp); movq %rbp, 8(%rsp); movq %r12, 16(%rsp); movq %r13, 24(%rsp); movq %r14, 32(%rsp); movq %r15, 40(%rsp); movq %rdi, 48(%rsp); movq %rsi, 56(%rsp); movups %xmm6, 64(%rsp); movups %xmm7, 80(%rsp); movups %xmm8, 96(%rsp); movups %xmm9, 112(%rsp); movups %xmm10, 128(%rsp); movups %xmm11, 144(%rsp); movups %xmm12, 160(%rsp); movups %xmm13, 176(%rsp); movups %xmm14, 192(%rsp); movups %xmm15, 208(%rsp);\n\t"
"\n\t"
"\n\t"
"\n\t"
" xorpd %xmm0, %xmm0; movapd %xmm0, %xmm1; movapd %xmm0, %xmm2; movapd %xmm0, %xmm3\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rcx, %r10\n\t"
" movq %r9, %r11\n\t"
" movq 256 + 40(%rsp), %r12\n\t"
" sall $ 5, %r12d\n\t"
"\n\t"
" movq 256 + 48(%rsp), %r13\n\t"
" movq %r8, %r14\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_edge_dgemv_add_t_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_kernel_dgemv_add_t_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rdx, %r10\n\t"
" movq 256 + 56(%rsp), %r11\n\t"
" movq 256 + 64(%rsp), %r12\n\t"
"\n\t"
"\n\t"
" INNER_BLEND_T_SCALE_AB_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq 256 + 72(%rsp), %r10\n\t"
" movq 256 + 80(%rsp), %r11\n\t"
"\n\t"
"\n\t"
" INNER_STORE_4_VS_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq (%rsp), %rbx; movq 8(%rsp), %rbp; movq 16(%rsp), %r12; movq 24(%rsp), %r13; movq 32(%rsp), %r14; movq 40(%rsp), %r15; movq 48(%rsp), %rdi; movq 56(%rsp), %rsi; movups 64(%rsp), %xmm6; movups 80(%rsp), %xmm7; movups 96(%rsp), %xmm8; movups 112(%rsp), %xmm9; movups 128(%rsp), %xmm10; movups 144(%rsp), %xmm11; movups 160(%rsp), %xmm12; movups 176(%rsp), %xmm13; movups 192(%rsp), %xmm14; movups 208(%rsp), %xmm15; addq $256, %rsp;\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .globl kernel_dgemv_nt_4_lib4; .def kernel_dgemv_nt_4_lib4; .scl 2; .type 32; .endef; kernel_dgemv_nt_4_lib4:\n\t"
"\n\t"
" subq $256, %rsp; movq %rbx, (%rsp); movq %rbp, 8(%rsp); movq %r12, 16(%rsp); movq %r13, 24(%rsp); movq %r14, 32(%rsp); movq %r15, 40(%rsp); movq %rdi, 48(%rsp); movq %rsi, 56(%rsp); movups %xmm6, 64(%rsp); movups %xmm7, 80(%rsp); movups %xmm8, 96(%rsp); movups %xmm9, 112(%rsp); movups %xmm10, 128(%rsp); movups %xmm11, 144(%rsp); movups %xmm12, 160(%rsp); movups %xmm13, 176(%rsp); movups %xmm14, 192(%rsp); movups %xmm15, 208(%rsp);\n\t"
"\n\t"
"\n\t"
"\n\t"
" xorpd %xmm0, %xmm0; movapd %xmm0, %xmm1; movapd %xmm0, %xmm2; movapd %xmm0, %xmm3\n\t"
"\n\t"
"\n\t"
" movq %rdx, %r10\n\t"
" movddup 0(%r10), %xmm15\n\t"
"\n\t"
" movq 256 + 48(%rsp), %r10\n\t"
"\n\t"
" movddup 0(%r10), %xmm4\n\t"
" mulpd %xmm15, %xmm4\n\t"
" movddup 8(%r10), %xmm5\n\t"
" mulpd %xmm15, %xmm5\n\t"
" movddup 16(%r10), %xmm6\n\t"
" mulpd %xmm15, %xmm6\n\t"
" movddup 24(%r10), %xmm7\n\t"
" mulpd %xmm15, %xmm7\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rcx, %r10\n\t"
" movq %r9, %r11\n\t"
" movq 256 + 40(%rsp), %r12\n\t"
" sall $ 5, %r12d\n\t"
"\n\t"
" movq 256 + 56(%rsp), %r13\n\t"
" movq 256 + 80(%rsp), %r14\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_kernel_dgemv_add_nt_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %r8, %r10\n\t"
" movq 256 + 64(%rsp), %r11\n\t"
" movq 256 + 72(%rsp), %r12\n\t"
"\n\t"
"\n\t"
" INNER_BLEND_T_SCALE_AB_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq 256 + 88(%rsp), %r10\n\t"
"\n\t"
"\n\t"
" INNER_STORE_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq (%rsp), %rbx; movq 8(%rsp), %rbp; movq 16(%rsp), %r12; movq 24(%rsp), %r13; movq 32(%rsp), %r14; movq 40(%rsp), %r15; movq 48(%rsp), %rdi; movq 56(%rsp), %rsi; movups 64(%rsp), %xmm6; movups 80(%rsp), %xmm7; movups 96(%rsp), %xmm8; movups 112(%rsp), %xmm9; movups 128(%rsp), %xmm10; movups 144(%rsp), %xmm11; movups 160(%rsp), %xmm12; movups 176(%rsp), %xmm13; movups 192(%rsp), %xmm14; movups 208(%rsp), %xmm15; addq $256, %rsp;\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .globl kernel_dgemv_nt_4_vs_lib4; .def kernel_dgemv_nt_4_vs_lib4; .scl 2; .type 32; .endef; kernel_dgemv_nt_4_vs_lib4:\n\t"
"\n\t"
" subq $256, %rsp; movq %rbx, (%rsp); movq %rbp, 8(%rsp); movq %r12, 16(%rsp); movq %r13, 24(%rsp); movq %r14, 32(%rsp); movq %r15, 40(%rsp); movq %rdi, 48(%rsp); movq %rsi, 56(%rsp); movups %xmm6, 64(%rsp); movups %xmm7, 80(%rsp); movups %xmm8, 96(%rsp); movups %xmm9, 112(%rsp); movups %xmm10, 128(%rsp); movups %xmm11, 144(%rsp); movups %xmm12, 160(%rsp); movups %xmm13, 176(%rsp); movups %xmm14, 192(%rsp); movups %xmm15, 208(%rsp);\n\t"
"\n\t"
"\n\t"
"\n\t"
" xorpd %xmm0, %xmm0\n\t"
" movapd %xmm0, %xmm1\n\t"
" movapd %xmm0, %xmm2\n\t"
" movapd %xmm0, %xmm3\n\t"
"\n\t"
" movapd %xmm0, %xmm4\n\t"
" movapd %xmm0, %xmm5\n\t"
" movapd %xmm0, %xmm6\n\t"
" movapd %xmm0, %xmm7\n\t"
"\n\t"
"\n\t"
" movq %rdx, %r10\n\t"
" movddup 0(%r10), %xmm15\n\t"
"\n\t"
" movq 256 + 48(%rsp), %r10\n\t"
" movq 256 + 96(%rsp), %r11\n\t"
"\n\t"
" movddup 0(%r10), %xmm4\n\t"
" mulpd %xmm15, %xmm4\n\t"
" cmpl $ 2, %r11d\n\t"
" jl 0f\n\t"
" movddup 8(%r10), %xmm5\n\t"
" mulpd %xmm15, %xmm5\n\t"
" cmpl $ 3, %r11d\n\t"
" jl 0f\n\t"
" movddup 16(%r10), %xmm6\n\t"
" mulpd %xmm15, %xmm6\n\t"
" je 0f\n\t"
" movddup 24(%r10), %xmm7\n\t"
" mulpd %xmm15, %xmm7\n\t"
"0:\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rcx, %r10\n\t"
" movq %r9, %r11\n\t"
" movq 256 + 40(%rsp), %r12\n\t"
" sall $ 5, %r12d\n\t"
"\n\t"
" movq 256 + 56(%rsp), %r13\n\t"
" movq 256 + 80(%rsp), %r14\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_kernel_dgemv_add_nt_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %r8, %r10\n\t"
" movq 256 + 64(%rsp), %r11\n\t"
" movq 256 + 72(%rsp), %r12\n\t"
"\n\t"
"\n\t"
" INNER_BLEND_T_SCALE_AB_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq 256 + 88(%rsp), %r10\n\t"
" movq 256 + 96(%rsp), %r11\n\t"
"\n\t"
"\n\t"
" INNER_STORE_4_VS_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq (%rsp), %rbx; movq 8(%rsp), %rbp; movq 16(%rsp), %r12; movq 24(%rsp), %r13; movq 32(%rsp), %r14; movq 40(%rsp), %r15; movq 48(%rsp), %rdi; movq 56(%rsp), %rsi; movups 64(%rsp), %xmm6; movups 80(%rsp), %xmm7; movups 96(%rsp), %xmm8; movups 112(%rsp), %xmm9; movups 128(%rsp), %xmm10; movups 144(%rsp), %xmm11; movups 160(%rsp), %xmm12; movups 176(%rsp), %xmm13; movups 192(%rsp), %xmm14; movups 208(%rsp), %xmm15; addq $256, %rsp;\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t"
" .p2align 4,,15\n\t"
" .globl kernel_dsymv_l_4_lib4; .def kernel_dsymv_l_4_lib4; .scl 2; .type 32; .endef; kernel_dsymv_l_4_lib4:\n\t"
"\n\t"
" subq $256, %rsp; movq %rbx, (%rsp); movq %rbp, 8(%rsp); movq %r12, 16(%rsp); movq %r13, 24(%rsp); movq %r14, 32(%rsp); movq %r15, 40(%rsp); movq %rdi, 48(%rsp); movq %rsi, 56(%rsp); movups %xmm6, 64(%rsp); movups %xmm7, 80(%rsp); movups %xmm8, 96(%rsp); movups %xmm9, 112(%rsp); movups %xmm10, 128(%rsp); movups %xmm11, 144(%rsp); movups %xmm12, 160(%rsp); movups %xmm13, 176(%rsp); movups %xmm14, 192(%rsp); movups %xmm15, 208(%rsp);\n\t"
"\n\t"
"\n\t"
"\n\t"
" xorpd %xmm0, %xmm0; movapd %xmm0, %xmm1; movapd %xmm0, %xmm2; movapd %xmm0, %xmm3\n\t"
"\n\t"
"\n\t"
" movq %rdx, %r10\n\t"
" movddup 0(%r10), %xmm15\n\t"
"\n\t"
" movq 256 + 40(%rsp), %r10\n\t"
"\n\t"
" movddup 0(%r10), %xmm4\n\t"
" mulpd %xmm15, %xmm4\n\t"
" movddup 8(%r10), %xmm5\n\t"
" mulpd %xmm15, %xmm5\n\t"
" movddup 16(%r10), %xmm6\n\t"
" mulpd %xmm15, %xmm6\n\t"
" movddup 24(%r10), %xmm7\n\t"
" mulpd %xmm15, %xmm7\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rcx, %r10\n\t"
" movq %r8, %r11\n\t"
" movq %r9, %r12\n\t"
" sall $ 5, %r12d\n\t"
" movq 256 + 40(%rsp), %r13\n\t"
" movq 256 + 48(%rsp), %r14\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_edge_dsymv_add_nt_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" call inner_kernel_dgemv_add_nt_4_lib4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq %rdx, %r10\n\t"
" movq 256 + 48(%rsp), %r11\n\t"
"\n\t"
"\n\t"
" INNER_BLEND_T_SCALE_A1_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq 256 + 48(%rsp), %r10\n\t"
"\n\t"
"\n\t"
" INNER_STORE_4_LIB4\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
"\n\t"
" movq (%rsp), %rbx; movq 8(%rsp), %rbp; movq 16(%rsp), %r12; movq 24(%rsp), %r13; movq 32(%rsp), %r14; movq 40(%rsp), %r15; movq 48(%rsp), %rdi; movq 56(%rsp), %rsi; movups 64(%rsp), %xmm6; movups 80(%rsp), %xmm7; movups 96(%rsp), %xmm8; movups 112(%rsp), %xmm9; movups 128(%rsp), %xmm10; movups 144(%rsp), %xmm11; movups 160(%rsp), %xmm12; movups 176(%rsp), %xmm13; movups 192(%rsp), %xmm14; movups 208(%rsp), %xmm15; addq $256, %rsp;\n\t"
"\n\t"
" ret\n\t"
"\n\t"
"\n\t");