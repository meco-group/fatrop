int get_horizon_length();
int get_nxk(const int k);
int get_initial_xk(double *xk, const int k);
int get_initial_uk(double *uk, const int k);
int get_nuk(const int k);
int get_ngk(const int k);
int get_n_stage_params_k(const int k);
int get_default_stage_paramsk(double *stage_params, const int k);
int get_n_global_parmas();
int get_default_global_params(double *global_params);
int get_ng_ineq_k(const int k);
int get_lower_boundsk(double *lower_bounds, const int k);
int get_upper_boundsk(double *upper_bounds, const int k);
int eval_BAbtk(
    const double *states_kp1,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k);
int eval_RSQrqtk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *lam_dyn_k,
    const double *lam_eq_k,
    const double *lam_eq_ineq_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k);
int eval_Ggtk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k);
int eval_Ggt_ineqk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k);
int eval_bk(
    const double *states_kp1,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k);
int eval_gk(
    const double *states_k,
    const double *inputs_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k);
int eval_gineqk(
    const double *states_k,
    const double *inputs_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k);
int eval_rqk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k);
int eval_Lk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k);

// function pointers
// typedef int (*intf_get_horizon_lenght)(void);

// typedef int (*intf_get_nxk)(const int k);

// typedef int (*intf_get_initial_xk)(double *xk, const int k);

// typedef int (*intf_get_initial_uk)(double *uk, const int k);

// typedef int (*intf_get_nuk)(const int k);

// typedef int (*intf_get_ngk)(const int k);

// typedef int (*intf_get_n_stage_params_k)(const int k);

// typedef int (*intf_get_default_stage_paramsk)(double *stage_params, const int k);

// typedef int (*intf_get_n_global_parmas)();

// typedef int (*intf_get_default_global_params)(double *global_params);

// typedef int (*intf_get_ng_ineq_k)(const int k);

// typedef int (*intf_get_lower_boundsk)(double *lower_bounds, const int k);

// typedef int (*intf_get_upper_boundsk)(double *upper_bounds, const int k);

// typedef int (*intf_get_horizon_length)();

// typedef int (*intf_eval_BAbtk)(
//     const double *states_kp1,
//     const double *inputs_k,
//     const double *states_k,
//     const double *stage_params_k,
//     const double *global_params_k,
//     double *res,
//     const int k);

// typedef int (*intf_eval_RSQrqtk)(
//     const double *objective_scale,
//     const double *inputs_k,
//     const double *states_k,
//     const double *lam_dyn_k,
//     const double *lam_eq_k,
//     const double *lam_eq_ineq_k,
//     const double *stage_params_k,
//     const double *global_params_k,
//     double *res,
//     const int k);

// typedef int (*intf_eval_Ggtk)(
//     const double *inputs_k,
//     const double *states_k,
//     const double *stage_params_k,
//     const double *global_params_k,
//     double *res,
//     const int k);

// typedef int (*intf_eval_Ggt_ineqk)(
//     const double *inputs_k,
//     const double *states_k,
//     const double *stage_params_k,
//     const double *global_params_k,
//     double *res,
//     const int k);

// typedef int (*intf_eval_bk)(
//     const double *states_kp1,
//     const double *inputs_k,
//     const double *states_k,
//     const double *stage_params_k,
//     const double *global_params_k,
//     double *res,
//     const int k);

// typedef int (*intf_eval_gk)(
//     const double *states_k,
//     const double *inputs_k,
//     const double *stage_params_k,
//     const double *global_params_k,
//     double *res,
//     const int k);

// typedef int (*intf_eval_gineqk)(
//     const double *states_k,
//     const double *inputs_k,
//     const double *stage_params_k,
//     const double *global_params_k,
//     double *res,
//     const int k);

// typedef int (*intf_eval_rqk)(
//     const double *objective_scale,
//     const double *inputs_k,
//     const double *states_k,
//     const double *stage_params_k,
//     const double *global_params_k,
//     double *res,
//     const int k);

// typedef int (*intf_eval_Lk)(
//     const double *objective_scale,
//     const double *inputs_k,
//     const double *states_k,
//     const double *stage_params_k,
//     const double *global_params_k,
//     double *res,
//     const int k);