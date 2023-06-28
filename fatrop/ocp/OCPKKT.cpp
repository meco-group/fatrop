#include "ocp/OCPKKT.hpp"
using namespace fatrop;

OCPKKTMemory::OCPKKTMemory(const OCPDims &dims) : K(dims.K),
                                                  nu(dims.nu),
                                                  nx(dims.nx),
                                                  ng(dims.ng),
                                                  ng_ineq(dims.ng_ineq),
                                                  RSQrqt(dims.nu + dims.nx + 1, dims.nu + dims.nx, dims.K),
                                                  BAbt(dims.nu + dims.nx + 1, rotate(dims.nx, 1), dims.K),
                                                  Ggt(dims.nu + dims.nx + 1, dims.ng, dims.K),
                                                  Ggt_ineq(dims.nu + dims.nx + 1, dims.ng_ineq, dims.K),
                                                  aux(dims){};
OCPKKTMemory::OCPAux::OCPAux(const OCPDims &dims) : ux_offs(offsets(dims.nx + dims.nu)),
                                                    g_offs(offsets(dims.ng)),
                                                    dyn_offs(offsets(rotate(dims.nx,1))),
                                                    dyn_eq_offs(offsets(rotate(dims.nx, 1)) + sum(dims.ng)),
                                                    g_ineq_offs(offsets(dims.ng_ineq) + (sum(dims.nx) - dims.nx.get(0) + sum(dims.ng))),
                                                    ineq_offs(offsets(dims.ng_ineq)),
                                                    max_nu(max(dims.nu)), max_nx(max(dims.nx)),
                                                    max_ng(max(dims.ng)),
                                                    max_ngineq(max(dims.ng_ineq)),
                                                    n_ineqs(sum(dims.ng_ineq)){};