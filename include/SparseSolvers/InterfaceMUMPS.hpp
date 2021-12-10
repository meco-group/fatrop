#ifndef INTERFACEMUMPSINCLUDED
#define INTERFACEMUMPSINCLUDED
#include "Interface.hpp"
extern "C"
{
#include "dmumps_c.h"
#include "mpi.h"
}
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

namespace fatrop
{
    class InterfaceMUMPS : public SparseSolverInterface
    {
    public:
        InterfaceMUMPS(const int nnz, const int dim, const vector<triplet> &tripl) : SparseSolverInterface(nnz, dim, tripl)
        {
            // fortran indexing convention
            ai = ai + 1;
            aj = aj + 1;
        };
        void preprocess()
        {
            MUMPS_INT n = dim_;
            MUMPS_INT8 nnz = nnz_;
            // mumps takes upper part of matrix!
            MUMPS_INT *irn = aj.data();
            MUMPS_INT *jcn = ai.data();

            int argc = 1;
            char *name = "c_example";
            char **argv;
            argv = &name;
            ierr = MPI_Init(&argc, &argv);
            ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

            /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
            id.comm_fortran = USE_COMM_WORLD;
            id.par = 1;
            id.sym = 2;

            id.job = JOB_INIT;
            dmumps_c(&id);
            /* Define the problem on the host */
            if (myid == 0)
            {
                id.n = n;
                id.nnz = nnz;
                id.irn = irn;
                id.jcn = jcn;
            }
#define ICNTL(I) icntl[(I)-1] /* macro sdsdfs.t. indices match documentation */
            /* No outputs */
            id.ICNTL(1) = -1;
            id.ICNTL(2) = -1;
            id.ICNTL(3) = -1;
            id.ICNTL(4) = 0;
            id.ICNTL(7) = 5;
            id.ICNTL(28) = 1;
            id.job = 1;
            dmumps_c(&id);
            if (id.infog[0] < 0)
            {
                printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                       myid, id.infog[0], id.infog[1]);
                error = 1;
            }
        };
        void solve(const vector<triplet> &A, vector<double> &rhsvec)
        {
            MUMPS_INT n = dim_;
            MUMPS_INT8 nnz = nnz_;
            // mumps takes upper part of matrix!
            MUMPS_INT *irn = aj.data();
            MUMPS_INT *jcn = ai.data();
            /* Define the problem on the host */
            if (myid == 0)
            {
                id.n = n;
                id.nnz = nnz;
                id.irn = irn;
                id.jcn = jcn;
            }

            vector<double> a;
            assert(((int)A.size()) == nnz);
            for (int i = 0; i < nnz; i++)
            {
                // check wether right structure
                assert(ai.at(i) == A.at(i).row() + 1);
                assert(aj.at(i) == A.at(i).col() + 1);
                a.push_back(A.at(i).value());
            }
            id.a = a.data();
            id.rhs = rhsvec.data();

            /* Call the MUMPS package (analyse, factorization and solve). */
            id.job = 2;
            dmumps_c(&id);
            if (id.infog[0] < 0)
            {
                printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                       myid, id.infog[0], id.infog[1]);
                error = 1;
            }
            id.job = 3;
            dmumps_c(&id);
            if (id.infog[0] < 0)
            {
                printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                       myid, id.infog[0], id.infog[1]);
                error = 1;
            }

            /* Terminate instance. */
            id.job = JOB_END;
            dmumps_c(&id);
            if (myid == 0)
            {
                if (!error)
                {
                    // printf("Solution is : (%8.2f  %8.2f)\n", rhs[0], rhs[1]);
                }
                else
                {
                    printf("An error has occured, please check error code returned by MUMPS.\n");
                }
            }
        }
        DMUMPS_STRUC_C id;
        int myid, ierr;
        int error = 0;
    };
} // namespace fatrop
#endif //INTERFACEMUMPSINCLUDED