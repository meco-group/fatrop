import re
class TranslationPattern:
    def __init__(self, fr, to):
        # remove all spaces from fr
        fr = re.sub(r"\s+", "", fr)
        to = re.sub(r"[ \t]+", "", to)
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
        print(pat)
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
    Ppt(dims.nx + 1, dims.nx, dims.K),
    Hh(dims.nx, dims.nx + 1, dims.K),
    AL(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nx)), 1),
    RSQrqt_tilde(dims.nu + dims.nx + 1, dims.nx + dims.nu, dims.K),
    Ggt_stripe(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nx + dims.nu)), 1),
    Ggt_tilde(dims.nu + dims.nx + 1, dims.nx + dims.nu, dims.K),
    GgLt(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nu + dims.nx)), 1),
    RSQrqt_hat(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nx + dims.nu)), 1),
    Llt(dims.nu + dims.nx + 1, dims.nu, dims.K),
    Llt_shift(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nu)), 1),
    GgIt_tilde(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
    GgLIt(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
    HhIt(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
    PpIt_hat(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
    LlIt(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
    Ggt_ineq_temp(vector<int>(1, max(dims.nu + dims.nx) + 1), vector<int>(1, max(dims.ng_ineq)), 1),
    """
    fr = "$A$($dims$ + 1, $arg1$, $N$)"
    to = "v_$A$($dims$, $N$)"
    transl1.add_pattern(TranslationPattern(fr, to))
    fr = "$A$(vector<int>($arg1$), $dim$, 1)"
    to = "v_$A$($dim$, 1)"
    transl1.add_pattern(TranslationPattern(fr, to))
    print(transl1.translate(text))

    transl = Translator()
    text = """
    GECP(nx + 1, ng_ineq, Ggt_ineq_p + K - 1, nu, 0, Ggt_ineq_temp_p, 0, 0);
    GEMM_NT(nu + nx + 1, nxp1,   nxp1,   1.0,    BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
    //COLSC(nx + 1, scaling_factor, Ggt_ineq_temp_p, 0, i);
    MATEL(Ggt_ineq_temp_p, nx, i) = grad_barrier + (scaling_factor)*MATEL(Ggt_ineq_p + K - 1, nu + nx, i);
    SYRK_LN_MN(nu + nx + 1, nu + nx, ng_ineq, 1.0, Ggt_ineq_temp_p, 0, 0, Ggt_ineq_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
    GETR(nx + 1, ng, Ggt_p + (K - 1), nu, 0, Hh_p + (K - 1), 0, 0);
    GECP(nu + nx + 1, ng, Ggt_p + k, 0, 0, Ggt_stripe_p, 0, 0);


    GEMM_NT(nu + nx + 1, nxp1, nxp1, 1.0, BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
    GEAD(1, nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, AL_p, nx + nu, 0);
    SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, 1.0, AL_p, 0, 0, BAbt_p + k, 0, 0, 1.0, RSQrqt_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
    GECP(nu + nx + 1, ng_ineq, Ggt_ineq_p + k, 0, 0, Ggt_ineq_temp_p, 0, 0);
    //COLSC(nu + nx, scaling_factor, Ggt_ineq_temp_p, 0, i);
    MATEL(Ggt_ineq_temp_p, nu + nx, i) = grad_barrier + (scaling_factor)*MATEL(Ggt_ineq_p + k, nu + nx, i);
    SYRK_LN_MN(nu + nx + 1, nu + nx, ng_ineq, 1.0, Ggt_ineq_temp_p, 0, 0, Ggt_ineq_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
    GECP(nu + nx + 1, ng, Ggt_p + k, 0, 0, Ggt_stripe_p, 0, 0);
    GEMM_NT(nu + nx + 1, Hp1_size, nxp1, 1.0, BAbt_p + k, 0, 0, Hh_p + (k + 1), 0, 0, 0.0, Ggt_stripe_p, 0, ng, Ggt_stripe_p, 0, ng);
    GEADTR(1, Hp1_size, 1.0, Hh_p + (k + 1), 0, nxp1, Ggt_stripe_p, nu + nx, ng);
    LU_FACT_transposed(gamma_k, nu + nx + 1, nu, rank_k, Ggt_stripe_p, Pl_p + k, Pr_p + k);
    GETR(nx + 1, gamma_k - rank_k, Ggt_stripe_p, nu, rank_k, Hh_p + k, 0, 0);
    TRSM_RLNN(nu - rank_k + nx + 0, rank_k, -1.0, Ggt_stripe_p, 0, 0, Ggt_stripe_p, rank_k, 0, Ggt_tilde_p + k, 0, 0);
    """

    fr = "GEMM_NT( $arg1$ + 1, $dim1$,   $dim2$,   $alpha$,   $B$, 0, 0, $A$, 0, 0, $beta$, $C$ , 0, $offsC$, $D$, 0, $offsD$);"
    to = "GEMV_T($dim1$, $dim2$, $alpha$, $A$, 0,0, v_$B$, 0,0, $beta$, v_$C$, $offsC$, v_$D$, $offsD$);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "COLSC($arg1$+ 1, $scaling_factor$, $A$, 0, $i$);"
    to = "VECEL(v_$A$, $i$) *= $scaling_factor$;"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "MATEL($A$, $arg2$, $i$) = $c$ + $scaling_factor$*MATEL($B$, $arg3$, $arg4$);"
    to = "VECEL(v_$A$, $i$) = $c$ + $scaling_factor$*t_VECEL(v_$B$, $arg4$);" 
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "SYRK_LN_MN( $arg1$ + 1, $dim1$,   $dim2$,   $alpha$,   $B$, 0, 0, $A$, 0, 0, $beta$, $C$ , 0, 0, $D$, 0, 0);"
    to = "GEMV_T($dim1$, $dim2$, $alpha$, $A$, 0,0, v_$B$, 0,0, $beta$, v_$C$, 0, v_$D$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GECP($arg1$ + 1, $dim$, $A$, $arg2$, 0, $B$, 0, 0);"
    to = "VECCP($dim$, v_$A$, 0, v_$B$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GETR($arg1$+ 1, $dim$, $A$, $arg2$, 0, $B$, 0, 0);"
    to = "VECCP($dim$, v_$A$, 0, v_$B$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GEAD(1, $dim$, $alpha$, $A$, $arg1$, 0, $B$, $arg2$, 0);"
    to = "AXPY($dim$, $alpha$, v_$A$, 0, v_$B$, 0, v_$B$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GEADTR(1, $dim$, 1.0, $A$, 0, $offsA$, $B$, $arg1$, $offsB$);"
    to = "AXPY($dim$, 1.0, v_$A$, $offsA$, v_$B$, $offsB$, v_$B$, $offsB$);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "LU_FACT_transposed($gamma$, $arg2$ + 1, $arg3$, $rank$ , $G$, $Pl$, $Pr$);"
    ### Ggt_stripe [U1_T \ L1T, L2T; U2T, etaT]
    to = """
    $Pl$ -> PV($rank$, v_$G$);
    // L1^-1 g_stipe[:rho]
    TRSV_UTN($rank$, $G$, 0,0, v_$G$, 0, v_$G$, 0);
    // -L2 L1^-1 g_stripe[:rho] + g_stripe[rho:]
    GEMV_T($arg3$ - $rank$ , $gamma$-$rank$, -1.0, $G$, 0, $rank$, v_$G$, 0, v_$G$, $rank$)
    """
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "GETR( $arg1$  + 1, $dim$, $A$, nu, rank_k, $B$, 0, 0);"
    to = "VECCP($dim$, v_$A$, rank_k, v_$B$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "TRSM_RLNN($arg1$ + 1, $dim$, $alpha$, Ggt_stripe_p, 0, 0, Ggt_stripe_p, rank_k, 0, Ggt_tilde_p + k, 0, 0);"
    print(transl.translate(text))
