%module BenchmarkAux
    %{
    #define SWIG_FILE_WITH_INIT
    #include "debug/RandomOCP.hpp"
    #include "aux/SmartPtr.hpp"
    #include <vector>
    using namespace fatrop;
    RefCountPtr<BFOCP> GenerateRandom(int nu, int nx, int ng, int K)
    {
    FatropVector<int> nu_ = vector<int>(K, nu);
    FatropVector<int> nx_ = vector<int>(K, nx);
    FatropVector<int> ng_ = vector<int>(K, ng);
        return new RandomOCP(nu_, nx_, ng_, K);
    }
    double& getEl(int dim2, double* mat, int ai, int aj)
    {
        return mat[dim2*ai+aj];
    }
    %}
    %include "numpy.i"



    %rename(BFOCPp) RefCountPtr<BFOCP>;

    %apply(int DIM1,int DIM2, double* IN_ARRAY2) {(int dim1, int dim2, double* mat)};
    class RefCountPtr<BFOCP> 
    {
        public:
    %extend {
        int GetRSQrqtk(int k, int dim1, int dim2, double* mat)
        {
            int nx = (*self)->get_nxk(k);
            int nu = (*self)->get_nuk(k);
            blasfeo_dmat bfmat;
            blasfeo_allocate_dmat(nu+nx+1, nu+nx, &bfmat);
            double dummy;
            (*self) ->  eval_RSQrqtk(&dummy,
                         &dummy,
                         &dummy,
                         &dummy,
                         &dummy,
                         &bfmat,
                         k) ;
            // blasfeo_print_dmat(nu+nx+1, nu+nx, &bfmat,0,0);
            for(int i = 0; i< nu+nx+1; i++){
            for(int j = 0; j< nu+nx; j++){
            getEl(dim2, mat,i,j) = MATEL(&bfmat,i,j);
            }
            }
            blasfeo_free_dmat(&bfmat);
            return 0;
        }
    }
    };
    RefCountPtr<BFOCP> GenerateRandom(int nu, int nx, int ng, int K);

    %init %{
        import_array();

    %}
