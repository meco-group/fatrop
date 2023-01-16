import re
class TranslationPattern:
    def __init__(self, fr, to):
        # remove all spaces from fr
        fr = re.sub(r"\s+", "", fr)
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
        line = re.sub(r"\s+", "", line)
        # find characters between prototype[i] and prototype[i+1]
        res = []
        for i in range(len(self.prototype) - 1):
            pat = re.escape(self.prototype[i]) + r"(.*?)" + re.escape(self.prototype[i+1])
            res.append(re.search(pat, line).groups()[0])
            line = ''.join([self.prototype[i+1]] + re.split(pat, line,1)[2:])
        return res
    def translate_line(self, line):
        line = line.group()
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
        res = re.sub(r"[ \t]+", "", input)
        for pattern in self.patterns:
            # get regex of pattern
            regex = pattern.get_regex()
            res = re.sub(regex, pattern.translate_line, res)
        return res

if __name__ == "__main__":
    transl = Translator()
    text = """GEMM_NT(nu + nx + 1, nxp1,   nxp1,   1.0,    BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
COLSC(nx + 1, scaling_factor, Ggt_ineq_temp_p, 0, i);
MATEL(Ggt_ineq_temp_p, nx, i) = grad_barrier + (scaling_factor)*MATEL(Ggt_ineq_p + K - 1, nu + nx, i);
    """

    fr = "GEMM_NT( $arg1$ + 1, $dim1$,   $dim2$,   $alpha$,   $B$, 0, 0, $A$, 0, 0, $beta$, $C$ , 0, 0, $D$, 0, 0);"
    to = "t_GEMV_N($dim1$, $dim2$, $alpha$, $A$, 0,0, v_$B$, 0,0, $beta$, v_$C$, 0, v_$D$, 0);"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "COLSC($arg1$+ 1, $scaling_factor$, $A$, 0, $i$);"
    to = "t_VECEL(v_$A$, $i$) *= $scaling_factor$;"
    transl.add_pattern(TranslationPattern(fr, to))
    fr = "MATEL($A$, $arg2$, $i$) = $c$ + $scaling_factor$*MATEL($B$, $arg3$, $arg4$);"
    to = "t_VECEL(v_$A$, $i$) = $c$ + $scaling_factor$*t_VECEL(_v$B$, $arg4$);" 
    transl.add_pattern(TranslationPattern(fr, to))

    print(transl.translate(text))
