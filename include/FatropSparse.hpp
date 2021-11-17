/**
 * @file FatropSparse.hpp
 * @author your name (you@domain.com)
 * @brief this file contains functions for representing a block-sparse KKT matrix, only used for DEBUG and TESTING purposes.
 * @version 0.1
 * @date 2021-11-17
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef FATROPSPARSEINCLUDED
#define FATROPSPARSEINCLUDED
#include "FatropOCPKKT.hpp"
#include "FatropLinearAlgebraEigen.hpp"
#include <vector>
#include <memory>
#include <cstdlib>
#include <string>
#include <iomanip> // std::setprecision
#include <eigen3/Eigen/Sparse>
using namespace std;
namespace fatrop
{
    class variable
    {
    public:
        variable(int size) : size(size){};
        int size;
        int offset;
        void set_grad(vector<double> &grad_)
        {
            grad.assign(size, 0.0);
            for (int i = 0; i < size; i++)
            {
                grad.at(i) = grad_.at(i);
            }
        };
        void add_rhs(vector<double> &rhs)
        {
            for (int i = 0; i < size; i++)
            {
                rhs.at(offset + i) = grad.at(i);
            }
        }
        vector<double> grad;
    };
    typedef std::shared_ptr<variable> var_sp;
    class KKT_matrix_base
    {
    };
    class matrix_vector_base
    {
    public:
        matrix_vector_base(Eig &mat, var_sp var) : fsm(mat), var(var){};
        Eig fsm;
        var_sp var;
    };

    class fatrop_expression
    {
    public:
        virtual bool is_matrix_vector() { return false; };
        virtual void add_to_mv_vec(vector<matrix_vector_base> &mv_vec) = 0;
        virtual int get_size() = 0;
    };
    typedef std::shared_ptr<fatrop_expression> fe_sp;
    class equation
    {
    public:
        equation(int size) : size(size){};
        int size;
        int offset;
        void add_expression(const fe_sp &express, vector<double> &rhs_)
        {
            rhs = rhs_;
            express->add_to_mv_vec(mv_veceq);
        }
        void add_triplets(vector<Eigen::Triplet<double>> &tripl, int offs_H)
        {
            for (unsigned long int i = 0; i < mv_veceq.size(); i++)
            {
                matrix_vector_base mvi = mv_veceq.at(i);
                int offs_var = mvi.var->offset;
                Eig &fsm = mvi.fsm;
                int m = fsm.nrows();
                int n = fsm.ncols();
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        double val = fsm.get_el(i, j);
                        if (val != 0.0)
                        {

                            tripl.push_back(Eigen::Triplet<double>(offs_H + offset + i, offs_var + j, val));
                        }
                    }
                }
            }
        }
        void add_rhs(std::vector<double> &rhs_, int offs_H)
        {
            for (int i = 0; i < size; i++)
            {
                rhs_.at(offs_H + offset + i) = rhs.at(i);
            }
        }
        std::vector<matrix_vector_base> mv_veceq;
        std::vector<double> rhs;
    };
    typename std::shared_ptr<equation> eq_sp;

    class fatrop_sum1 : public fatrop_expression
    {
    public:
        fatrop_sum1(const fe_sp &fe1, const fe_sp &fe2) : child1(fe1), child2(fe2){};
        void add_to_mv_vec(std::vector<matrix_vector_base> &mv_vec)
        {
            child1->add_to_mv_vec(mv_vec);
            child2->add_to_mv_vec(mv_vec);
        };
        int get_size()
        {
            return child1->get_size();
        }
        const fe_sp child1;
        const fe_sp child2;
    };
    fe_sp operator+(const fe_sp &fe1, const fe_sp &fe2)
    {
        fe_sp res = std::make_shared<fatrop_sum1>(fe1, fe2);
        return res;
    }
    class matrix_vector : public fatrop_expression, public matrix_vector_base
    {
    public:
        matrix_vector(Eig &mat, var_sp var) : matrix_vector_base(mat, var){};
        bool is_matrix_vector() { return true; };
        void add_to_mv_vec(std::vector<matrix_vector_base> &mv_vec)
        {
            matrix_vector_base *mvb = static_cast<matrix_vector_base *>(this);
            mv_vec.push_back(*mvb);
        }
        int get_size()
        {
            return fsm.nrows();
        }
    };
    fe_sp operator*(Eig &mat, var_sp var)
    {
        fe_sp res = std::make_shared<matrix_vector>(mat, var);
        return res;
    };

    class hess_block
    {
    public:
        hess_block(Eig &mat, var_sp var1, var_sp var2) : fsm(mat), var1(var1), var2(var2){};
        Eig fsm;
        var_sp var1;
        var_sp var2;
        void add_triplets(std::vector<Eigen::Triplet<double>> &tripl)
        {
            int offs_var1 = var1->offset;
            int offs_var2 = var2->offset;
            int m = fsm.nrows();
            int n = fsm.ncols();
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double val = fsm.get_el(i, j);
                    if (val != 0.0) // only add nonzero's
                    {
                        int ai = i + offs_var1;
                        int aj = i + offs_var2;
                        if (ai >= aj) // only lower triangular matrix
                        {
                            tripl.push_back(Eigen::Triplet<double>(ai, aj, val));
                        }
                        else
                        {
                            if (offs_var1 != offs_var2)
                            {
                                tripl.push_back(Eigen::Triplet<double>(aj, ai, val));
                            }
                        }
                    }
                }
            }
        }
    };

    class KKT_matrix : public KKT_matrix_base
    {
    public:
        var_sp get_variable(int size)
        {
            var_sp var = std::make_shared<variable>(size);
            variable_vec.push_back(var);
            return variable_vec.back();
        };
        equation &get_equation(int size)
        {
            equation eq(size);
            equation_vec.push_back(eq);
            return equation_vec.back();
        };
        void set_equation(fe_sp expr, std::vector<double> rhs)
        {
            equation &eq = get_equation(expr->get_size());
            eq.add_expression(expr, rhs);
        }
        void set_hess_block(Eig &mat, var_sp var1, var_sp var2)
        {
            hess_block hb(mat, var1, var2);
            hess_block_vec.push_back(hb);
        }

        std::vector<var_sp> variable_vec;
        std::vector<equation> equation_vec;
        std::vector<hess_block> hess_block_vec;
        void set_offsets()
        {
            int offs_curr = 0;
            for (long unsigned int i = 0; i < variable_vec.size(); i++)
            {
                variable_vec.at(i)->offset = offs_curr;
                offs_curr += variable_vec.at(i)->size;
            }
            offs_curr = 0;
            for (long unsigned int i = 0; i < equation_vec.size(); i++)
            {
                equation_vec.at(i).offset = offs_curr;
                offs_curr += equation_vec.at(i).size;
            }
        }
        void get_triplets(std::vector<Eigen::Triplet<double>> &tripl_vec)
        {
            this->set_offsets();
            for (long unsigned int i = 0; i < hess_block_vec.size(); i++)
            {
                hess_block_vec.at(i).add_triplets(tripl_vec);
            }
            int offs_vars = variable_vec.back()->offset + variable_vec.back()->size;
            // std::cout << "offs "<< offs_vars << std::endl;
            for (long unsigned int i = 0; i < equation_vec.size(); i++)
            {
                equation_vec.at(i).add_triplets(tripl_vec, offs_vars);
            }
        }
        std::vector<double> get_rhs()
        {
            std::vector<double> rhs;
            this->set_offsets();
            rhs.assign(get_size(), 0.0);
            for (long unsigned int i = 0; i < variable_vec.size(); i++)
            {
                variable_vec.at(i)->add_rhs(rhs);
            }
            int offs_vars = variable_vec.back()->offset + variable_vec.back()->size;
            // std::cout << "offs "<< offs_vars << std::endl;
            for (long unsigned int i = 0; i < equation_vec.size(); i++)
            {
                equation_vec.at(i).add_rhs(rhs, offs_vars);
            }
            return rhs;
        }
        int get_size()
        {
            set_offsets();
            return equation_vec.back().offset + equation_vec.back().size + variable_vec.back()->offset + variable_vec.back()->size;
        }

        void print(char *type_id)
        {

            std::vector<Eigen::Triplet<double>> testvec;
            get_triplets(testvec);
            //printing (matrix)
            if (std::string(type_id) == "matrix")
            {
                std::cout << " " << get_size() << " " << get_size() << std::endl
                          << " ";
                std::cout << std::fixed;

                Eigen::MatrixXd mat(get_size(), get_size() + 1);

                for (long unsigned int i = 0; i < testvec.size(); i++)
                {
                    // std::cout << testvec.at(i).col() + 1 << " " << testvec.at(i).row()  << " " << std::setprecision(20) << testvec.at(i).value() << std::endl
                    mat(testvec.at(i).row(), testvec.at(i).col()) = testvec.at(i).value();
                    mat(testvec.at(i).col(), testvec.at(i).row()) = testvec.at(i).value();
                }
                std::vector<double> rhs = get_rhs();
                for (int i = 0; i < get_size(); i++)
                {
                    mat(i, get_size()) = rhs.at(i);
                }
                std::cout << mat << std::endl;
            }
            //printing (ma57)
            if (std::string(type_id) == "ma57")
            {
                std::cout << " " << get_size() << " " << testvec.size() << std::endl
                          << " ";
                std::cout << std::fixed;

                for (long unsigned int i = 0; i < testvec.size(); i++)
                {
                    std::cout << testvec.at(i).col() + 1 << " " << testvec.at(i).row() + 1 << " " << std::setprecision(20) << testvec.at(i).value() << std::endl
                              << " ";
                }
                std::vector<double> rhs = get_rhs();
                for (int i = 0; i < get_size(); i++)
                {
                    std::cout << rhs.at(i) << " ";
                }
                std::cout << std::endl;
            }
            //printing (MUMPS)
            if (std::string(type_id) == "mumps")
            {
                std::cout << get_size() << "     :N \n"
                          << testvec.size() << "     :NZ ";
                std::cout << std::fixed;

                for (long unsigned int i = 0; i < testvec.size(); i++)
                {
                    std::cout << std::endl
                              << testvec.at(i).col() + 1 << " " << testvec.at(i).row() + 1 << " " << std::setprecision(20) << testvec.at(i).value();
                }
                std::cout << "            :values";
                std::vector<double> rhs = get_rhs();
                for (int i = 0; i < get_size(); i++)
                {
                    std::cout << std::endl
                              << rhs.at(i);
                }
                std::cout << "            :RHS";
                std::cout << std::endl;
            }
        }
    };
} //%namespace fatrop

#endif //FATROPSPARSEINCLUDED