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
        InterfaceMUMPS(const int nnz, const int dim, KKT_matrix &KKT_m) : SparseSolverInterface(nnz, dim, KKT_m){};
        void preprocess()
        {
            MUMPS_INT n = 2;
            MUMPS_INT8 nnz = 2;
            MUMPS_INT irn[] = {1, 2};
            MUMPS_INT jcn[] = {1, 2};
            double a[2];
            double rhs[2];

            int myid, ierr;

            int error = 0;

            int argc = 1;
            char *name = "c_example";
            char **argv;
            argv = &name;
            ierr = MPI_Init(&argc, &argv);
            ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
            rhs[0] = 1.0;
            rhs[1] = 4.0;
            a[0] = 1.0;
            a[1] = 2.0;

            /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
            id.comm_fortran = USE_COMM_WORLD;
            id.par = 1;
            id.sym = 0;
            id.job = JOB_INIT;
            dmumps_c(&id);
            /* Define the problem on the host */
            if (myid == 0)
            {
                id.n = n;
                id.nnz = nnz;
                id.irn = irn;
                id.jcn = jcn;
                id.a = a;
                id.rhs = rhs;
            }
#define ICNTL(I) icntl[(I)-1] /* macro sdsdfs.t. indices match documentation */
            /* No outputs */
            id.ICNTL(1) = -1;
            id.ICNTL(2) = -1;
            id.ICNTL(3) = -1;
            id.ICNTL(4) = 0;

            /* Call the MUMPS package (analyse, factorization and solve). */
            id.job = 6;
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
                    printf("Solution is : (%8.2f  %8.2f)\n", rhs[0], rhs[1]);
                }
                else
                {
                    printf("An error has occured, please check error code returned by MUMPS.\n");
                }
            }
        };
        void solve(const vector<triplet> A, const vector<double> rhs){};
        DMUMPS_STRUC_C id;
    };
} // namespace fatrop
#endif //INTERFACEMUMPSINCLUDED