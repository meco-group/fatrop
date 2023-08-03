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
#include <iostream>
#include <fstream>
#include <string>
#include <ocp/StageOCPApplication.hpp>
using namespace fatrop;
using namespace std;
int main(int argc, char **argv)
{
    if (argc == 3)
    {
        //// dynamic memory allocation  
        StageOCPApplication app = StageOCPApplicationFactory::from_rockit_interface(argv[1], argv[2]);
        auto eval_expression = app.get_expression("control_u").at_t0();
        vector<double> u0_result(eval_expression.size());
        ///  no dynamic memory allocation
        app.optimize();
        // ///  retrieve solution
        app.last_solution().evaluate(eval_expression, u0_result);
        ///  initialize next solver run 
        app.set_initial(app.last_solution());
        ///  change solver options
        app.set_option("tol", 1e-6);
        app.set_option("mu_init", 1e-5);
        app.set_option("bound_push", 1e-7); // all *_bound_push variables
        app.set_option("warm_start_mult_bound_push", 1e-7); 
        app.set_option("accept_every_trial_step", false);
        app.set_option("iterative_refinement", false); // fast_step_computation
        app.set_option("warm_start_init_point", true);
        app.set_option("print_level", 0);
        // cout << app -> GetOptions();
        app.optimize();
    }
    else
    {
            cout << "run me as RunFatrop f.so f.json" << endl;
    }
}