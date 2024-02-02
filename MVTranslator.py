import re
class TranslationPattern:
    def __init__(self, fr, to):
        # remove all spaces from fr
        fr = re.sub(r"\s+", "", fr)
        # to = re.sub(r"[ \t]+", "", to)
        self.fr = fr 
        self.to = to
        self.args = re.findall(r"\$\w+\$", self.fr)
        # set up the prototype
        self.prototype = []
        fr_left = fr
        for arg in self.args:
            res = re.split(re.escape(arg), fr_left)
            assert(len(res) == 2)
            self.prototype.append(res[0])
            fr_left = res[1]
        self.prototype.append(fr_left)
        pass
    def get_args(self, line):
        # remove all spaces from line
        # line = re.sub(r"\s+", "", line)
        # find characters between prototype[i] and prototype[i+1]
        res = []
        for i in range(len(self.prototype) - 1):
            pat = re.escape(self.prototype[i]) + r"(.*?)" + re.escape(self.prototype[i+1])
            res.append(re.search(pat, line).groups()[0])
            line = ''.join([self.prototype[i+1]] + re.split(pat, line,1)[2:])
        return res
    def translate_line(self, line):
        line = line.group()
        # remove all spaces from line
        line = re.sub(r"\s+", "", line)
        args_to = self.get_args(line)
        res = self.to
        for i in range(len(self.args)):
            # print(self.args[i])
            res = re.sub(re.escape(self.args[i]), args_to[i], res)
        return res
    def get_regex(self):
        # construct regex pattern from self.prototype
        pat =''.join([re.escape(proto) + r"(.*?)" for proto in self.prototype][:-1] + [re.escape(self.prototype[-1])])
        return pat

class Translator:
    def __init__(self):
        self.patterns = []
        pass
    def add_pattern(self, pattern):
        self.patterns.append(pattern)

    def translate(self, input):
        # remove all spaces from input
        # res = input
        res = re.sub(r"[ \t]+", "", input)
        for pattern in self.patterns:
            # get regex of pattern
            regex = pattern.get_regex()
            # create list of regex elements
            res = re.sub(regex, pattern.translate_line, res)
        return res

if __name__ == "__main__":
    transl1 = Translator()
    text = """

    fatrop_int K = OCP->K;
    // make variables local for efficiency
    OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    OCPMACRO(MAT *, Ggt_ineq, _p);
    SOLVERMACRO(MAT *, Ppt, _p);
    SOLVERMACRO(MAT *, AL, _p);
    SOLVERMACRO(MAT *, RSQrqt_tilde, _p);
    SOLVERMACRO(MAT *, Ggt_tilde, _p);
    SOLVERMACRO(MAT *, Llt, _p);
    SOLVERMACRO(MAT *, Llt_shift, _p);
    SOLVERMACRO(MAT *, LlIt, _p);
    SOLVERMACRO(MAT *, Ggt_ineq_temp, _p);
    SOLVERMACRO(VEC *, ux, _p);
    SOLVERMACRO(VEC *, lam, _p);
    SOLVERMACRO(VEC *, delta_s, _p);
    SOLVERMACRO(VEC *, sigma_total, _p);
    SOLVERMACRO(VEC *, gradb_total, _p);
    OCPMACRO(fatrop_int *, nu, _p);
    OCPMACRO(fatrop_int *, nx, _p);
    OCPMACRO(fatrop_int *, ng, _p);
    OCPMACRO(fatrop_int *, ng_ineq, _p);

    """
    fr = "$A$($dims$ + 1, $arg1$, $N$),"
    to = "FatropMemoryVecBF v_$A$($dims$, $N$);"
    transl1.add_pattern(TranslationPattern(fr, to))
    fr = "$A$(vector<int>($arg1$), $dim$, 1),"
    to = "FatropMemoryVecBF v_$A$($dim$, 1);"
    transl1.add_pattern(TranslationPattern(fr, to))
    print(transl1.translate(text))

    transl = Translator()
    text = """
    MAT *RSQrq_hat_curr_p;
    double delta_cmin1 = 1 / inertia_correction_c;
    fatrop_int *offs_ineq_p = (fatrop_int *)OCP->aux.ineq_offs.data();
    fatrop_int *offs_g_ineq_p = (fatrop_int *)OCP->aux.g_ineq_offs.data();

    /////////////// recursion ///////////////

    // last stage
    {
        const fatrop_int nx = nx_p[K - 1];
        const fatrop_int nu = nu_p[K - 1]; // this should be zero but is included here in case of misuse
        const fatrop_int ng = ng_p[K - 1];
        const fatrop_int ng_ineq = ng_ineq_p[K - 1];
        const fatrop_int offs_ineq_k = offs_ineq_p[K - 1];
        // const fatrop_int offs_g_ineq_k = offs_g_ineq_p[K - 1];
        // Pp_Km1 <- Qq_Km1
        GECP(nx + 1, nx, RSQrqt_p + (K - 1), nu, nu, Ppt_p + K - 1, 0, 0);
        DIARE(nx, inertia_correction_w, Ppt_p + K - 1, 0, 0);
        //// inequalities
        if (ng_ineq > 0)
        {
            GECP(nx + 1, ng_ineq, Ggt_ineq_p + K - 1, nu, 0, Ggt_ineq_temp_p, 0, 0);
            for (fatrop_int i = 0; i < ng_ineq; i++)
            {
                double scaling_factor = VECEL(sigma_total_p, offs_ineq_k + i) + inertia_correction_w;
                double grad_barrier = VECEL(gradb_total_p, offs_ineq_k + i);
                COLSC(nx + 1, scaling_factor, Ggt_ineq_temp_p, 0, i);
                MATEL(Ggt_ineq_temp_p, nx, i) = grad_barrier + (scaling_factor)*MATEL(Ggt_ineq_p + K - 1, nu + nx, i);
            }
            // add the penalty
            SYRK_LN_MN(nx + 1, nx, ng_ineq, 1.0, Ggt_ineq_temp_p, 0, 0, Ggt_ineq_p + K - 1, nu, 0, 1.0, Ppt_p + K - 1, 0, 0, Ppt_p + K - 1, 0, 0);
            // TRTR_L(nx, Ppt_p + K - 1, 0, 0, Ppt_p + K - 1, 0, 0);
        }

        GECP(nx + 1, ng, Ggt_p + K - 1, nu, 0, Ggt_tilde_p + K - 1, 0, 0); // needless operation because feature not implemented yet
        SYRK_LN_MN(nx + 1, nx, ng, delta_cmin1, Ggt_tilde_p + K - 1, 0, 0, Ggt_tilde_p + K - 1, 0, 0, 1.0, Ppt_p + K - 1, 0, 0, Ppt_p + K - 1, 0, 0);
        TRTR_L(nx, Ppt_p + K - 1, 0, 0, Ppt_p + K - 1, 0, 0);
    }
    for (fatrop_int k = K - 2; k >= 0; --k)
    {
        const fatrop_int nu = nu_p[k];
        const fatrop_int nx = nx_p[k];
        const fatrop_int nxp1 = nx_p[k + 1];
        const fatrop_int ng = ng_p[k];
        const fatrop_int ng_ineq = ng_ineq_p[k];
        const fatrop_int offs_ineq_k = offs_ineq_p[k];
        // const fatrop_int offs_g_ineq_k = offs_g_ineq_p[k];
        //////// SUBSDYN
        {
            // AL <- [BAb]^T_k P_kp1
            GEMM_NT(nu + nx + 1, nxp1, nxp1, 1.0, BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
            // AL[-1,:] <- AL[-1,:] + p_kp1^T
            GEAD(1, nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, AL_p, nx + nu, 0);
            // RSQrqt_stripe <- AL[BA] + RSQrqt
            SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, 1.0, AL_p, 0, 0, BAbt_p + k, 0, 0, 1.0, RSQrqt_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
            // RSQrqt + 1/d_c At@A
            DIARE(nu + nx, inertia_correction_w, RSQrqt_tilde_p + k, 0, 0);
            SYRK_LN_MN(nu + nx + 1, nu + nx, ng, delta_cmin1, Ggt_p + k, 0, 0, Ggt_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);

            //// inequalities
            if (ng_ineq > 0)
            {
                GECP(nu + nx + 1, ng_ineq, Ggt_ineq_p + k, 0, 0, Ggt_ineq_temp_p, 0, 0);
                for (fatrop_int i = 0; i < ng_ineq; i++)
                {
                    double scaling_factor = VECEL(sigma_total_p, offs_ineq_k + i) + inertia_correction_w;
                    double grad_barrier = VECEL(gradb_total_p, offs_ineq_k + i);
                    COLSC(nu + nx, scaling_factor, Ggt_ineq_temp_p, 0, i);
                    MATEL(Ggt_ineq_temp_p, nu + nx, i) = grad_barrier + (scaling_factor)*MATEL(Ggt_ineq_p + k, nu + nx, i);
                }
                // add the penalty
                SYRK_LN_MN(nu + nx + 1, nu + nx, ng_ineq, 1.0, Ggt_ineq_temp_p, 0, 0, Ggt_ineq_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
            }
        }
        //////// TRANSFORM_AND_SUBSEQ
        {
            RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
        }
        //////// SCHUR
        {
            // DLlt_k = [chol(R_hatk); Llk@chol(R_hatk)^-T]
            POTRF_L_MN(nu + nx + 1, nu, RSQrq_hat_curr_p, 0, 0, Llt_p + k, 0, 0);
            if (!check_reg(nu, Llt_p + k, 0, 0))
                return 1;
            // Pp_k = Qq_hatk - L_k^T @ Ll_k
            // SYRK_LN_MN(nx+1, nx, nu-rank_k, -1.0,Llt_p+k, nu-rank_k,0, Llt_p+k, nu-rank_k,0, 1.0, RSQrq_hat_curr_p, nu-rank_k, nu-rank_k,Pp+k,0,0); // feature not implmented yet
            GECP(nx + 1, nu, Llt_p + k, nu, 0, Llt_shift_p, 0, 0); // needless operation because feature not implemented yet
            SYRK_LN_MN(nx + 1, nx, nu, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, RSQrq_hat_curr_p, nu, nu, Ppt_p + k, 0, 0);
        }
        TRTR_L(nx, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
    }
    //////// FIRST_STAGE
    {
        const fatrop_int nx = nx_p[0];
        {
            POTRF_L_MN(nx + 1, nx, Ppt_p, 0, 0, LlIt_p, 0, 0);
            if (!check_reg(nx, LlIt_p, 0, 0))
                return 2;
        }
    }
    ////// FORWARD_SUBSTITUTION:
    // first stage
    {
        const fatrop_int nx = nx_p[0];
        const fatrop_int nu = nu_p[0];
        // calculate xIb
        ROWEX(nx, -1.0, LlIt_p, nx, 0, ux_p, nu);
        // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
        TRSV_LTN(nx, LlIt_p, 0, 0, ux_p, nu, ux_p, nu);
    }
    fatrop_int *offs_ux = (fatrop_int *)OCP->aux.ux_offs.data();
    fatrop_int *offs_g = (fatrop_int *)OCP->aux.g_offs.data();
    fatrop_int *offs_dyn_eq = (fatrop_int *)OCP->aux.dyn_eq_offs.data();
    // other stages
    // for (fatrop_int k = 0; k < K - 1; k++)
    // fatrop_int dyn_eqs_ofs = offs_g[K - 1] + ng_p[K - 1]; // this value is incremented at end of recursion
    for (fatrop_int k = 0; k < K - 1; k++)
    {
        const fatrop_int nx = nx_p[k];
        const fatrop_int nu = nu_p[k];
        const fatrop_int nxp1 = nx_p[k + 1];
        const fatrop_int nup1 = nu_p[k + 1];
        const fatrop_int offsp1 = offs_ux[k + 1];
        const fatrop_int offs = offs_ux[k];
        const fatrop_int offs_dyn_eq_k = offs_dyn_eq[k];
        /// calculate ukb_tilde
        // -Lkxk - lk
        ROWEX(nu, -1.0, Llt_p + k, nu + nx, 0, ux_p, offs);
        // assume aliasing of last two eliments is allowed
        GEMV_T(nx, nu, -1.0, Llt_p + k, nu, 0, ux_p, offs + nu, 1.0, ux_p, offs, ux_p, offs);
        TRSV_LTN(nu, Llt_p + k, 0, 0, ux_p, offs, ux_p, offs);
        // calculate xkp1
        ROWEX(nxp1, 1.0, BAbt_p + k, nu + nx, 0, ux_p, offsp1 + nup1);
        GEMV_T(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, ux_p, offs, 1.0, ux_p, offsp1 + nup1, ux_p, offsp1 + nup1);
        // calculate lam_dyn xp1
        ROWEX(nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, lam_p, offs_dyn_eq_k);
        GEMV_T(nxp1, nxp1, 1.0, Ppt_p + (k + 1), 0, 0, ux_p, offsp1 + nup1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
        // // calculate lam_eq xk
        // ROWEX(ng, -1.0, Ggt_p + k, 0, 0, lam_p, offs_g_k);
        // GEMV_T(nu + nx, ng, delta_cmin1, Ggt_p + k, 0, 0, ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
    }
    // calculate lam_eq xk
    for (fatrop_int k = 0; k < K; k++)
    {
        const fatrop_int nx = nx_p[k];
        const fatrop_int nu = nu_p[k];
        const fatrop_int ng = ng_p[k];
        const fatrop_int offs = offs_ux[k];
        const fatrop_int offs_g_k = offs_g[k];
        ROWEX(ng, delta_cmin1, Ggt_p + k, nu + nx, 0, lam_p, offs_g_k);
        GEMV_T(nu + nx, ng, delta_cmin1, Ggt_p + k, 0, 0, ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
    }

    for (fatrop_int k = 0; k < K; k++)
    {
        const fatrop_int nx = nx_p[k];
        const fatrop_int nu = nu_p[k];
        const fatrop_int ng_ineq = ng_ineq_p[k];
        const fatrop_int offs = offs_ux[k];
        const fatrop_int offs_g_ineq_k = offs_g_ineq_p[k];
        const fatrop_int offs_ineq_k = offs_ineq_p[k];
        if (ng_ineq > 0)
        {
            // calculate delta_s
            ROWEX(ng_ineq, 1.0, Ggt_ineq_p + k, nu + nx, 0, delta_s_p, offs_ineq_k);
            // GEMV_T(nu + nx, ng_ineq, 1.0, Ggt_ineq_p + k, 0, 0, ux_p, offs, 1.0, delta_s_p, offs_ineq_k, delta_s_p, offs_ineq_k);
            GEMV_T(nu + nx, ng_ineq, 1.0, Ggt_ineq_p + k, 0, 0, ux_p, offs, 1.0, delta_s_p, offs_ineq_k, delta_s_p, offs_ineq_k);
            // calculate lamineq
            for (fatrop_int i = 0; i < ng_ineq; i++)
            {
                double scaling_factor = VECEL(sigma_total_p, offs_ineq_k + i) + inertia_correction_w;
                double grad_barrier = VECEL(gradb_total_p, offs_ineq_k + i);
                double ds = VECEL(delta_s_p, offs_ineq_k + i);
                VECEL(lam_p, offs_g_ineq_k + i) = scaling_factor * ds + grad_barrier;
            }
        }
    }

"""

    fr = "GEMM_NT( $arg1$ + 1, $dim1$,   $dim2$,   $alpha$,   $B$, 0, 0, $A$, 0, 0, $beta$, $C$ , 0, $offsC$, $D$, 0, $offsD$);"
    to = "GEMV_N($dim1$, $dim2$, $alpha$, $A$, 0,0, v_$B$, 0, $beta$, v_$C$, $offsC$, v_$D$, $offsD$);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "COLSC($arg1$+ 1, $scaling_factor$, $A$, 0, $i$);"
    to = "VECEL(v_$A$, $i$) *= $scaling_factor$;"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "MATEL($A$, $arg2$, $i$) = $c$ + $scaling_factor$*MATEL($B$, $arg3$, $arg4$);"
    to = "VECEL(v_$A$, $i$) = $c$ + $scaling_factor$*t_VECEL(v_$B$, $arg4$);" 
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "SYRK_LN_MN( $arg1$ + 1, $dim1$,   $dim2$,   $alpha$,   $B$, 0, 0, $A$, $argx$, 0, $beta$, $C$ , $arg_x$, $offs_C$, $D$, 0, 0);"
    to = "GEMV_N($dim1$, $dim2$, $alpha$, $A$, 0,0, v_$B$, 0,0, $beta$, v_$C$, $offs_C$, v_$D$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GECP($arg1$ + 1, $dim$, $A$, $arg2$, $offsA$, $B$, 0, $offsB$);"
    to = "VECCP($dim$, v_$A$, $offsA$, v_$B$, $offsB$);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GETR($arg1$+ 1, $dim$, $A$, $arg2$, 0, $B$, 0, 0);"
    to = "VECCP($dim$, v_$A$, 0, v_$B$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GEAD(1, $dim$, $alpha$, $A$, $arg1$, 0, $B$, $arg2$, 0);"
    to = "AXPY($dim$, $alpha$, v_$A$, 0, v_$B$, 0, v_$B$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GEADTR(1, $dim$, 1.0, $A$, 0, $offsA$, $B$, $arg1$, $offsB$);"
    to = "AXPY($dim$, 1.0, v_$A$, 0, v_$B$, $offsB$, v_$B$, $offsB$);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "LU_FACT_transposed($gamma$, $arg2$ + 1, $arg3$, $rank$ , $G$, $Pl$, $Pr$);"
    ### Ggt_stripe [U1_T \ L1T, L2T; U2T, etaT]
    to = """
    $Pl$ -> PV($rank$, v_$G$, 0);
    // L1^-1 g_stipe[:rho]
    TRSV_UTU($rank$, $G$, 0,0, v_$G$, 0, v_$G$, 0);
    // -L2 L1^-1 g_stripe[:rho] + g_stripe[rho:]
    GEMV_T($rank$ , $gamma$-$rank$, -1.0, $G$, 0, $rank$, v_$G$, 0, 1.0, v_$G$, $rank$, v_$G$, $rank$)
    """
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GETR( $arg1$  + 1, $dim$, $A$, nu, rank_k, $B$, 0, 0);"
    to = "VECCP($dim$, v_$A$, rank_k, v_$B$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "TRSM_RLNN($arg1$ + 1, $dim$, $alpha$, $A$, 0, 0, $B$, $arg2$, 0, $C$, 0, 0);"
    to = "VECPSC($dim$, $alpha$, v_$B$,0, v_$C$,0); TRSV_LTN($dim$, $A$, 0, 0, v_$C$, 0, v_$C$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "$P$->PM($rank$, $A$);"
    to = "$P$->PV($rank$, v_$A$);" 
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "$P$->MPt($rank$, $A$);"
    to = "$P$->VPt($rank$, v_$A$);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "POTRF_L_MN($arg1$ + 1, $dim$, $A$, 0, 0, $B$, 0, 0);"
    to = "TRSV_LNN($dim$, $B$, 0, 0, v_$B$, 0, v_$B$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GEMM_NT($dim1$, $arg1$ + 1, $dim2$, $alpha$, $A$, 0, 0, $B$, $arg2$, 0, $beta$, $C$, 0, 0, $D$, 0, 0);"
    to = "GEMV_N($dim1$, $dim2$, $alpha$, $A$, 0, 0, v_$B$, 0, $beta$, v_$C$, 0, v_$D$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GETR($dim$, $arg1$ + 1, $A$, 0, 0, $B$, 0, 0);"
    to = "VECCP($dim$, v_$A$, 0, v_$B$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "ROWEX($dim$, 1.0, $A$, $arg1$, $offsA$, $b$, $offs_b$);"
    to = "VECCP($dim$, v_$A$, $offsA$, $b$, $offs_b$);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "ROWEX($dim$, $alpha$, $A$, $arg1$, $offsA$, $b$, $offs_b$);"
    to = "VECCPSC($dim$, $alpha$, v_$A$, $offsA$, $b$, $offs_b$);"
    transl.add_pattern(TranslationPattern(fr, to))

    print(transl.translate(text))