/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
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
                                                    max_nu(maxel(dims.nu)), max_nx(maxel(dims.nx)),
                                                    max_ng(maxel(dims.ng)),
                                                    max_ngineq(maxel(dims.ng_ineq)),
                                                    n_ineqs(sum(dims.ng_ineq)){};