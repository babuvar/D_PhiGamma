// -*- C++ -*-                                                               
//                                                                                                                                            

#include "belle.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cfloat>

#include "basf/module.h"
#include "particle/Particle.h"
#include "tuple/BelleTupleManager.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  class BelleEvent;
  class BasfOutput;
  class BelleHistogram;

  class d0rad: public Module {

  public:

    d0rad( void );
    ~d0rad( void ) {};
    void init( int* );
    void term( void );
    void disp_stat( const char* );
    void hist_def( void );
    void event( BelleEvent*, int* );
    void begin_run( BelleEvent*, int* );
    void end_run( BelleEvent*, int* );
    void other( int*, BelleEvent*, int* );

void FillTuple(Particle& Dst_cand, BelleTuple *nt, int A[100], int B[100], int & n, int IdMatrix[10000][10], int & nIM);

void FillIDHEP(int a, int A[100], int B[100], int & n);

void PrintGenTable(Particle& Dst_cand,int GenTableFlag,int & counts);

void Transfer(int Temp[10],int IdMatrix[10000][10],int & nIM);

void Print(int IdMatrix[10000][10],int & nIM);

void read(int tmpmat[10000][10],int & tmpn);

void Append(int tmpmat[10000][10], int & tmpn, int IdMatrix[10000][10],int & nIM);

    void with_imp_cut(std::vector<Particle> &list);
    HepVector param_at_ip(Particle &p);

  private:
    int dataType;

    //Particle Types                                                                                                                               
    static const Ptype m_ptypeD0;
    static const Ptype m_ptypeD0B;
    static const Ptype m_ptypeDstarP;
    static const Ptype m_ptypeDstarM;

    // Event histograms            
    BelleTuple *nt_check;
    BelleTuple *nt_d0kpkm;

  };



#if defined(BELLE_NAMESPACE)
} // namespace Belle                                                                                                                             
#endif


