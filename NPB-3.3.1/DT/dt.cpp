/*************************************************************************
 *                                                                       *
 *        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3       *
 *                                                                       *
 *                                  D T					 *
 *                                                                       *
 *************************************************************************
 *                                                                       *
 *   This benchmark is part of the NAS Parallel Benchmark 3.3 suite.     *
 *                                                                       *
 *   Permission to use, copy, distribute and modify this software        *
 *   for any purpose with or without fee is hereby granted.  We          *
 *   request, however, that all derived work reference the NAS           *
 *   Parallel Benchmarks 3.3. This software is provided "as is"          *
 *   without express or implied warranty.                                *
 *                                                                       *
 *   Information on NPB 3.3, including the technical report, the         *
 *   original specifications, source code, results and information       *
 *   on how to submit new results, is available at:                      *
 *                                                                       *
 *          http:  www.nas.nasa.gov/Software/NPB                         *
 *                                                                       *
 *   Send comments or suggestions to  npb@nas.nasa.gov                   *
 *   Send bug reports to              npb-bugs@nas.nasa.gov              *
 *                                                                       *
 *         NAS Parallel Benchmarks Group                                 *
 *         NASA Ames Research Center                                     *
 *         Mail Stop: T27A-1                                             *
 *         Moffett Field, CA   94035-1000                                *
 *                                                                       *
 *         E-mail:  npb@nas.nasa.gov                                     *
 *         Fax:     (650) 604-3957                                       *
 *                                                                       *
 *************************************************************************
 *                                                                       *
 *   Author: M. Frumkin							 *						 *
 *                                                                       *
 *************************************************************************/

#include <string>
#include <libdash.h>
#include <iostream>
#include <unistd.h>

using namespace std;

#include "mpi.h"
#include "npbparams.h"
#include "DGraph.h"


/* Added For Dash */
//static int maxInDegree = 0;
#define fielddim 4
#define MAX_FEATURE_LEN ((NUM_SAMPLES) * (fielddim) * (2))

typedef struct DGNodeInfo_s {
  int id = -1;
  int feature_len = 0;
  int inDegree = 0;
  int outDegree = 0;
  int inArc[SMALL_BLOCK_SIZE] = {0};
  int outArc[SMALL_BLOCK_SIZE] = {0};
  char name[SMALL_BLOCK_SIZE];

} DGNodeInfo;

template <typename GraphData>
class Graph
{
public:
  Graph(
    string const name_, size_t numNodes, GraphData &data_):
    m_name(name_),
    m_nodes(numNodes, dash::BLOCKED),
    m_data(data_)
  {}

  dash::Array<DGNodeInfo> & getNodes() {
    return m_nodes;
  }

  GraphData & getData() {
    return m_data;
  }

  string const getName() const {
    return m_name;
  }
private:
  string const m_name;
  dash::Array<DGNodeInfo> m_nodes;
  GraphData &m_data;
};
/* End of this section */

#ifndef CLASS
#define CLASS 'S'
#define NUM_PROCS            1
#endif

int      passed_verification;
extern "C" {
  void    timer_clear( int n );
  void    timer_start( int n );
  void    timer_stop( int n );
  double  timer_read( int n );
  double randlc( double *X, double *A );
  void c_print_results( char const  *name,
                        char   class_name,
                        int    n1,
                        int    n2,
                        int    n3,
                        int    niter,
                        int    nprocs_compiled,
                        int    nprocs_total,
                        double t,
                        double mops,
                        char  const  *optype,
                        int    passed_verification,
                        char const  *npbversion,
                        char const  *compiletime,
                        char  const  *mpicc,
                        char const  *clink,
                        char const  *cmpi_lib,
                        char const  *cmpi_inc,
                        char  const *cflags,
                        char const  *clinkflags );



}
int timer_on=0,timers_tot=64;

int verify(char const* bmname,double rnm2) {
  double verify_value=0.0;
  double epsilon=1.0E-8;
  char cls=CLASS;
  int verified=-1;
  if (cls != 'U') {
    if(cls=='S') {
      if(strstr(bmname,"BH")) {
        verify_value=30892725.0;
      } else if(strstr(bmname,"WH")) {
        verify_value=67349758.0;
      } else if(strstr(bmname,"SH")) {
        verify_value=58875767.0;
      } else {
        fprintf(stderr,"No such benchmark as %s.\n",bmname);
      }
      verified = 0;
    } else if(cls=='W') {
      if(strstr(bmname,"BH")) {
        verify_value = 4102461.0;
      } else if(strstr(bmname,"WH")) {
        verify_value = 204280762.0;
      } else if(strstr(bmname,"SH")) {
        verify_value = 186944764.0;
      } else {
        fprintf(stderr,"No such benchmark as %s.\n",bmname);
      }
      verified = 0;
    } else if(cls=='A') {
      if(strstr(bmname,"BH")) {
        verify_value = 17809491.0;
      } else if(strstr(bmname,"WH")) {
        verify_value = 1289925229.0;
      } else if(strstr(bmname,"SH")) {
        verify_value = 610856482.0;
      } else {
        fprintf(stderr,"No such benchmark as %s.\n",bmname);
      }
      verified = 0;
    } else if(cls=='B') {
      if(strstr(bmname,"BH")) {
        verify_value = 4317114.0;
      } else if(strstr(bmname,"WH")) {
        verify_value = 7877279917.0;
      } else if(strstr(bmname,"SH")) {
        verify_value = 1836863082.0;
      } else {
        fprintf(stderr,"No such benchmark as %s.\n",bmname);
        verified = 0;
      }
    } else if(cls=='C') {
      if(strstr(bmname,"BH")) {
        verify_value = 0.0;
      } else if(strstr(bmname,"WH")) {
        verify_value = 0.0;
      } else if(strstr(bmname,"SH")) {
        verify_value = 0.0;
      } else {
        fprintf(stderr,"No such benchmark as %s.\n",bmname);
        verified = -1;
      }
    } else if(cls=='D') {
      if(strstr(bmname,"BH")) {
        verify_value = 0.0;
      } else if(strstr(bmname,"WH")) {
        verify_value = 0.0;
      } else if(strstr(bmname,"SH")) {
        verify_value = 0.0;
      } else {
        fprintf(stderr,"No such benchmark as %s.\n",bmname);
      }
      verified = -1;
    } else {
      fprintf(stderr,"No such class as %c.\n",cls);
    }
    fprintf(stderr," %s L2 Norm = %f\n",bmname,rnm2);
    if(verified==-1) {
      fprintf(stderr," No verification was performed.\n");
    } else if( rnm2 - verify_value < epsilon &&
               rnm2 - verify_value > -epsilon) {  /* abs here does not work on ALTIX */
      verified = 1;
      fprintf(stderr," Deviation = %f\n",(rnm2 - verify_value));
    } else {
      verified = 0;
      fprintf(stderr," The correct verification value = %f\n",verify_value);
      fprintf(stderr," Got value = %f\n",rnm2);
    }
  } else {
    verified = -1;
  }
  return  verified;
}

int ipowMod(int a,long long int n,int md) {
  int seed=1,q=a,r=1;
  if(n<0) {
    fprintf(stderr,"ipowMod: exponent must be nonnegative exp=%lld\n",n);
    n=-n; /* temp fix */
    /*    return 1; */
  }
  if(md<=0) {
    fprintf(stderr,"ipowMod: module must be positive mod=%d",md);
    return 1;
  }
  if(n==0) return 1;
  while(n>1) {
    int n2 = n/2;
    if (n2*2==n) {
      seed = (q*q)%md;
      q=seed;
      n = n2;
    } else {
      seed = (r*q)%md;
      r=seed;
      n = n-1;
    }
  }
  seed = (r*q)%md;
  return seed;
}

DGraph *buildSH(char cls) {
  /*
    Nodes of the graph must be topologically sorted
    to avoid MPI deadlock.
  */
  DGraph *dg;
  int numSources=NUM_SOURCES; /* must be power of 2 */
  int numOfLayers=0,tmpS=numSources>>1;
  int firstLayerNode=0;
  DGArc *ar=NULL;
  DGNode *nd=NULL;
  int mask=0x0,ndid=0,ndoff=0;
  int i=0,j=0;
  char nm[BLOCK_SIZE];

  sprintf(nm,"DT_SH.%c",cls);
  dg=newDGraph(nm);

  while(tmpS>1) {
    numOfLayers++;
    tmpS>>=1;
  }
  for(i=0; i<numSources; i++) {
    sprintf(nm,"Source.%d",i);
    nd=newNode(nm);
    AttachNode(dg,nd);
  }
  for(j=0; j<numOfLayers; j++) {
    mask=0x00000001<<j;
    for(i=0; i<numSources; i++) {
      sprintf(nm,"Comparator.%d",(i+j*firstLayerNode));
      nd=newNode(nm);
      AttachNode(dg,nd);
      ndoff=i&(~mask);
      ndid=firstLayerNode+ndoff;
      ar=newArc(dg->node[ndid],nd);
      AttachArc(dg,ar);
      ndoff+=mask;
      ndid=firstLayerNode+ndoff;
      ar=newArc(dg->node[ndid],nd);
      AttachArc(dg,ar);
    }
    firstLayerNode+=numSources;
  }
  mask=0x00000001<<numOfLayers;
  for(i=0; i<numSources; i++) {
    sprintf(nm,"Sink.%d",i);
    nd=newNode(nm);
    AttachNode(dg,nd);
    ndoff=i&(~mask);
    ndid=firstLayerNode+ndoff;
    ar=newArc(dg->node[ndid],nd);
    AttachArc(dg,ar);
    ndoff+=mask;
    ndid=firstLayerNode+ndoff;
    ar=newArc(dg->node[ndid],nd);
    AttachArc(dg,ar);
  }
  return dg;

  /*
  for (int idx = 0; idx < dg->numNodes; ++idx) {
    if (dg->node[idx]->inDegree > maxInDegree) {
      maxInDegree = dg->node[idx]->inDegree;
    }
  }
  */
}
DGraph *buildWH(char cls) {
  /*
    Nodes of the graph must be topologically sorted
    to avoid MPI deadlock.
  */
  int i=0,j=0;
  int numSources=NUM_SOURCES,maxInDeg=4;
  int numLayerNodes=numSources,firstLayerNode=0;
  int totComparators=0;
  int numPrevLayerNodes=numLayerNodes;
  int id=0,sid=0;
  DGraph *dg;
  DGNode *nd=NULL,*source=NULL,*tmp=NULL,*snd=NULL;
  DGArc *ar=NULL;
  char nm[BLOCK_SIZE];
  //global maxInDegree for creation of Dash Array
  //maxInDegree = maxInDeg;

  sprintf(nm,"DT_WH.%c",cls);
  dg=newDGraph(nm);

  for(i=0; i<numSources; i++) {
    sprintf(nm,"Sink.%d",i);
    nd=newNode(nm);
    AttachNode(dg,nd);
  }
  totComparators=0;
  numPrevLayerNodes=numLayerNodes;
  while(numLayerNodes>maxInDeg) {
    numLayerNodes=numLayerNodes/maxInDeg;
    if(numLayerNodes*maxInDeg<numPrevLayerNodes)numLayerNodes++;
    for(i=0; i<numLayerNodes; i++) {
      sprintf(nm,"Comparator.%d",totComparators);
      totComparators++;
      nd=newNode(nm);
      id=AttachNode(dg,nd);
      for(j=0; j<maxInDeg; j++) {
        sid=i*maxInDeg+j;
        if(sid>=numPrevLayerNodes) break;
        snd=dg->node[firstLayerNode+sid];
        ar=newArc(dg->node[id],snd);
        AttachArc(dg,ar);
      }
    }
    firstLayerNode+=numPrevLayerNodes;
    numPrevLayerNodes=numLayerNodes;
  }
  source=newNode("Source");
  AttachNode(dg,source);
  for(i=0; i<numPrevLayerNodes; i++) {
    nd=dg->node[firstLayerNode+i];
    ar=newArc(source,nd);
    AttachArc(dg,ar);
  }

  for(i=0; i<dg->numNodes/2; i++) { /* Topological sorting */
    tmp=dg->node[i];
    dg->node[i]=dg->node[dg->numNodes-1-i];
    dg->node[i]->id=i;
    dg->node[dg->numNodes-1-i]=tmp;
    dg->node[dg->numNodes-1-i]->id=dg->numNodes-1-i;
  }
  return dg;
}
DGraph *buildBH(char cls) {
  /*
    Nodes of the graph must be topologically sorted
    to avoid MPI deadlock.
  */
  int i=0,j=0;
  int numSources=NUM_SOURCES,maxInDeg=4;
  int numLayerNodes=numSources,firstLayerNode=0;
  DGraph *dg;
  DGNode *nd=NULL, *snd=NULL, *sink=NULL;
  DGArc *ar=NULL;
  int totComparators=0;
  int numPrevLayerNodes=numLayerNodes;
  int id=0, sid=0;
  char nm[BLOCK_SIZE];

  //global maxInDegree for creation of Dash Array
  //maxInDegree = maxInDeg;

  sprintf(nm,"DT_BH.%c",cls);
  dg=newDGraph(nm);

  for(i=0; i<numSources; i++) {
    sprintf(nm,"Source.%d",i);
    nd=newNode(nm);
    AttachNode(dg,nd);
  }
  while(numLayerNodes>maxInDeg) {
    numLayerNodes=numLayerNodes/maxInDeg;
    if(numLayerNodes*maxInDeg<numPrevLayerNodes)numLayerNodes++;
    for(i=0; i<numLayerNodes; i++) {
      sprintf(nm,"Comparator.%d",totComparators);
      totComparators++;
      nd=newNode(nm);
      id=AttachNode(dg,nd);
      for(j=0; j<maxInDeg; j++) {
        sid=i*maxInDeg+j;
        if(sid>=numPrevLayerNodes) break;
        snd=dg->node[firstLayerNode+sid];
        ar=newArc(snd,dg->node[id]);
        AttachArc(dg,ar);
      }
    }
    firstLayerNode+=numPrevLayerNodes;
    numPrevLayerNodes=numLayerNodes;
  }
  sink=newNode("Sink");
  AttachNode(dg,sink);
  for(i=0; i<numPrevLayerNodes; i++) {
    nd=dg->node[firstLayerNode+i];
    ar=newArc(nd,sink);
    AttachArc(dg,ar);
  }
  return dg;
}

typedef struct Arr_s {
  int len;
  double* val = nullptr;

  Arr_s(int len_):
    len(len_)
  {
    val = new double[len];
  }

  ~Arr_s()
  {
    delete[] val;
  }

} Arr;

/*
Arr *newArr(int len) {
  Arr *arr = new Arr(len);
  return arr;
}

void arrShow(DGNode_Feat* a) {
  if(!a) fprintf(stderr,"-- NULL array\n");
  else {
    fprintf(stderr,"-- length=%d\n",a->len);
  }
}
*/
double CheckVal(Arr const& feat) {
  double csum=0.0;
  int i=0;
  for(i=0; i<feat.len; i++) {
    csum+=feat.val[i]*feat.val[i]/feat.len; /* The truncation does not work since
                                                  result will be 0 for large len  */
  }
  return csum;
}
int GetFNumDPar(int* mean, int* stdev) {
  *mean=NUM_SAMPLES;
  *stdev=STD_DEVIATION;
  return 0;
}
int GetFeatureNum(char const*mbname,int id) {
  double tran=314159265.0;
  double A=2*id+1;
  double denom=randlc(&tran,&A);
  char cval='S';
  int mean=NUM_SAMPLES,stdev=128;
  int rtfs=0,len=0;
  GetFNumDPar(&mean,&stdev);
  rtfs=ipowMod((int)(1/denom)*(int)cval,(long long int) (2*id+1),2*stdev);
  if(rtfs<0) rtfs=-rtfs;
  len=mean-stdev+rtfs;
  return len;
}
void RandomFeatures(char const* bmname, int fdim, double * const feat, DGNodeInfo& nd ) {
  int len=GetFeatureNum(bmname,nd.id)*fdim;
  nd.feature_len=len;
  int nxg=2,nyg=2,nzg=2,nfg=5;
  int nx=421,ny=419,nz=1427,nf=3527;
  long long int expon=(len*(nd.id+1))%3141592;
  int seedx=ipowMod(nxg,expon,nx),
      seedy=ipowMod(nyg,expon,ny),
      seedz=ipowMod(nzg,expon,nz),
      seedf=ipowMod(nfg,expon,nf);
  int i=0;
  if(timer_on) {
    timer_clear(nd.id+1);
    timer_start(nd.id+1);
  }
  for(i=0; i<len; i+=fdim) {
    seedx=(seedx*nxg)%nx;
    seedy=(seedy*nyg)%ny;
    seedz=(seedz*nzg)%nz;
    seedf=(seedf*nfg)%nf;
    feat[i]=seedx;
    feat[i+1]=seedy;
    feat[i+2]=seedz;
    feat[i+3]=seedf;
  }
  if(timer_on) {
    timer_stop(nd.id+1);
    fprintf(stderr,"** RandomFeatures time in node %d = %f\n",nd.id,timer_read(nd.id+1));
  }
}
void Resample(Arr *a,int blen) {
  long long int i=0,j=0,jlo=0,jhi=0;
  double avval=0.0;
  double* nval= new double[blen];
  for(i=0; i<blen; i++) nval[i]=0.0;
  for(i=1; i<a->len-1; i++) {
    jlo=(int)(0.5*(2*i-1)*(blen/a->len));
    jhi=(int)(0.5*(2*i+1)*(blen/a->len));

    avval=a->val[i]/(jhi-jlo+1);
    for(j=jlo; j<=jhi; j++) {
      nval[j]+=avval;
    }
  }
  nval[0]=a->val[0];
  nval[blen-1]=a->val[a->len-1];
  delete[] a->val;
  a->val = nval;
  a->len=blen;
}
void WindowFilter(Arr *a, Arr *b,int w) {
  int i=0,j=0,k=0;
  double rms0=0.0,rms1=0.0,rmsm1=0.0;
  double weight=((double) (w+1))/(w+2);

  w+=1;
  if(timer_on) {
    timer_clear(w);
    timer_start(w);
  }
  if(a->len<b->len) Resample(a,b->len);
  if(a->len>b->len) Resample(b,a->len);
  for(i=fielddim; i<a->len-fielddim; i+=fielddim) {
    rms0=(a->val[i]-b->val[i])*(a->val[i]-b->val[i])
         +(a->val[i+1]-b->val[i+1])*(a->val[i+1]-b->val[i+1])
         +(a->val[i+2]-b->val[i+2])*(a->val[i+2]-b->val[i+2])
         +(a->val[i+3]-b->val[i+3])*(a->val[i+3]-b->val[i+3]);
    j=i+fielddim;
    rms1=(a->val[j]-b->val[j])*(a->val[j]-b->val[j])
         +(a->val[j+1]-b->val[j+1])*(a->val[j+1]-b->val[j+1])
         +(a->val[j+2]-b->val[j+2])*(a->val[j+2]-b->val[j+2])
         +(a->val[j+3]-b->val[j+3])*(a->val[j+3]-b->val[j+3]);
    j=i-fielddim;
    rmsm1=(a->val[j]-b->val[j])*(a->val[j]-b->val[j])
          +(a->val[j+1]-b->val[j+1])*(a->val[j+1]-b->val[j+1])
          +(a->val[j+2]-b->val[j+2])*(a->val[j+2]-b->val[j+2])
          +(a->val[j+3]-b->val[j+3])*(a->val[j+3]-b->val[j+3]);
    k=0;
    if(rms1<rms0) {
      k=1;
      rms0=rms1;
    }
    if(rmsm1<rms0) k=-1;
    if(k==0) {
      j=i+fielddim;
      a->val[i]=weight*b->val[i];
      a->val[i+1]=weight*b->val[i+1];
      a->val[i+2]=weight*b->val[i+2];
      a->val[i+3]=weight*b->val[i+3];
    } else if(k==1) {
      j=i+fielddim;
      a->val[i]=weight*b->val[j];
      a->val[i+1]=weight*b->val[j+1];
      a->val[i+2]=weight*b->val[j+2];
      a->val[i+3]=weight*b->val[j+3];
    } else { /*if(k==-1)*/
      j=i-fielddim;
      a->val[i]=weight*b->val[j];
      a->val[i+1]=weight*b->val[j+1];
      a->val[i+2]=weight*b->val[j+2];
      a->val[i+3]=weight*b->val[j+3];
    }
  }
  if(timer_on) {
    timer_stop(w);
    fprintf(stderr,"** WindowFilter time in node %d = %f\n",(w-1),timer_read(w));
  }
}

int SendResults(DGNodeInfo const& nd) {
  int i;
  for(i=0; i < nd.outDegree; ++i) {
    MPI_Send(&nd.feature_len, 1, MPI_INT, nd.outArc[i], 10, MPI_COMM_WORLD);
  }
  return 1;
}

void CombineStreams(Graph<dash::Array<double>> &graph, DGNodeInfo & nd) {

  if (nd.inDegree == 0) return;

  int len, i, predRank;

  auto &array = graph.getData();
  const std::array<dash::default_index_t, 1> coords { {0} };

  Arr *pred_feat;
  Arr *resfeat = new Arr(NUM_SAMPLES * fielddim);

  for(i = 0; i < nd.inDegree; ++i) {
    predRank = nd.inArc[i];
    if(predRank != nd.id) {
      MPI_Recv(&len, 1 , MPI_INT, predRank , 10 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      pred_feat = new Arr(len);

      //Copy Array
      auto gindex = array.pattern().global(predRank, coords);
      auto copy_start_idx = gindex[0];
      auto copy_end_idx = copy_start_idx + len;

      dash::copy(array.begin() + copy_start_idx, array.begin() + copy_end_idx, &(pred_feat->val[0]));

      WindowFilter(resfeat, pred_feat, nd.id);

      delete pred_feat;
    } else {
      pred_feat = new Arr(nd.feature_len);
      memcpy(pred_feat->val, array.lbegin(), nd.feature_len*sizeof(double));
      WindowFilter(resfeat, pred_feat, nd.id);
      delete pred_feat;
    }
  }

  const int inD = nd.inDegree;

  for(i = 0; i < resfeat->len; ++i) {
    resfeat->val[i]=((int)resfeat->val[i])/inD;
  }

  std::copy(&(resfeat->val[0]), &(resfeat->val[resfeat->len]), array.lbegin());
  nd.feature_len = resfeat->len;
}

double Reduce(Arr const& a,int w) {
  double retv=0.0;
  if(timer_on) {
    timer_clear(w);
    timer_start(w);
  }
  retv=(int)(w*CheckVal(a));/* The casting needed for node
                               and array dependent verifcation */
  if(timer_on) {
    timer_stop(w);
    fprintf(stderr,"** Reduce time in node %d = %f\n",(w-1),timer_read(w));
  }
  return retv;
}

double ReduceStreams(Graph<dash::Array<double>> & graph, DGNodeInfo const& nd) {
  double csum=0.0;
  int i=0;
  int len;
  double retv=0.0;
  auto &array = graph.getData();
  Arr* feat = nullptr;
  const std::array<dash::default_index_t, 1> coords { {0} };
  

  for(i=0; i<nd.inDegree; ++i) {
    auto predRank = nd.inArc[i];
    if(predRank != nd.id) {
      MPI_Recv(&len, 1 , MPI_INT, predRank , 10 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      feat = new Arr(len);

      //Copy Array
      auto gindex = array.pattern().global(predRank, coords);
      auto copy_start_idx = gindex[0];
      auto copy_end_idx = copy_start_idx + len;
      dash::copy(array.begin() + copy_start_idx, array.begin() + copy_end_idx, feat->val);

      //Calculate Checksum
      csum+=Reduce(*feat,(nd.id+1));

      delete feat;
    } else {
      csum+=Reduce(*feat,(nd.id+1));
    }
  }

  if(nd.inDegree > 0) csum=(((long long int)csum)/nd.inDegree);
  retv=(nd.id+1)*csum;

  return retv;
}

int ProcessNodes(Graph<dash::Array<double>> & graph) {
  double chksum=0.0;
  int verified=0,tag;
  double rchksum=0.0;
  MPI_Status status;

  auto &nodes = graph.getNodes();
  auto const numNodes = nodes.size();
  DGNodeInfo &me = nodes.local[0];

  if (strlen(me.name) == 0) return 0;

  if(strstr(me.name,"Source")) {
    double * lfeat = graph.getData().lbegin();
    RandomFeatures(graph.getName().c_str(), fielddim, lfeat, me);
    SendResults(me);
  }
  else if(strstr(me.name, "Sink")) {
    chksum=ReduceStreams(graph, me);
    tag = numNodes + me.id;
    MPI_Send(&chksum,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
  }
  else {
    CombineStreams(graph, me);
    SendResults(me);
  }

  if(me.id == 0) { // Report node
    rchksum=0.0;
    for(size_t idx = 0; idx < numNodes; ++idx) {
      DGNodeInfo nd = nodes[idx];
      if(!(strstr(nd.name, "Sink"))) continue;
      tag=numNodes+nd.id; // make these to avoid clash with arc tags
      MPI_Recv(&rchksum,1,MPI_DOUBLE,nd.id,tag,MPI_COMM_WORLD,&status);
      chksum+=rchksum;
    }
    verified=verify(graph.getName().c_str(),chksum);
  }

  return verified;
}

int main(int argc,char **argv ) {
  int my_rank,comm_size;
  int i;
  DGraph *dg=NULL;
  int verified=0, featnum=0;
  double bytes_sent=2.0,tot_time=0.0;
  dash::init(&argc, &argv);

  my_rank = dash::myid();
  comm_size = dash::size();

  if(argc!=2||
      (  strncmp(argv[1],"BH",2)!=0
         &&strncmp(argv[1],"WH",2)!=0
         &&strncmp(argv[1],"SH",2)!=0
      )
    ) {
    if(my_rank==0) {
      fprintf(stderr,"** Usage: mpirun -np N ../bin/dt.S GraphName\n");
      fprintf(stderr,"** Where \n   - N is integer number of MPI processes\n");
      fprintf(stderr,"   - S is the class S, W, or A \n");
      fprintf(stderr,"   - GraphName is the communication graph name BH, WH, or SH.\n");
      fprintf(stderr,"   - the number of MPI processes N should not be be less than \n");
      fprintf(stderr,"     the number of nodes in the graph\n");
    }
    dash::finalize();
    exit(0);
  }
  if(strncmp(argv[1],"BH",2)==0) {
    dg=buildBH(CLASS);
  } else if(strncmp(argv[1],"WH",2)==0) {
    dg=buildWH(CLASS);
  } else if(strncmp(argv[1],"SH",2)==0) {
    dg=buildSH(CLASS);
  }

  if(dg->numNodes>comm_size) {
    if(my_rank==0) {
      fprintf(stderr,"**  The number of MPI processes should not be less than \n");
      fprintf(stderr,"**  the number of nodes in the graph\n");
      fprintf(stderr,"**  Number of MPI processes = %d\n",comm_size);
      fprintf(stderr,"**  Number nodes in the graph = %d\n",dg->numNodes);
    }
    dash::finalize();
    exit(0);
  }

  if(timer_on&&dg->numNodes+1>timers_tot) {
    timer_on=0;
    if(my_rank==0)
      fprintf(stderr,"Not enough timers. Node timeing is off. \n");
  }

  for(i=0; i<dg->numNodes; i++) {
    dg->node[i]->address=i;
  }

  //Create Dash Array
  //Each process represents 1 node
  //Each node has maxInDegree predecessors
  //We need one feature array for each predecessor
  dash::Array<double> features(dash::size() * MAX_FEATURE_LEN);
  //dash::Array<DGNodeInfo> nodes(dash::size());
  Graph<dash::Array<double>> graph(dg->name, dash::size(), features);


  for(i = 0; i < dg->numNodes; i++) {
    if (my_rank == i) {
      int j;
      DGNode *nd = dg->node[i];
      //copy required info
      DGNodeInfo ndInfo;
      string nodeName(nd->name);
      ndInfo.id = nd->address;
      strcpy(ndInfo.name, nd->name);;
      ndInfo.inDegree = nd->inDegree;
      ndInfo.outDegree = nd->outDegree;

      for (j = 0; j < nd->inDegree; ++j)
      {
        ndInfo.inArc[j] = nd->inArc[j]->tail->address;
      }
      for (j = 0; j < nd->outDegree; ++j)
      {
        ndInfo.outArc[j] = nd->outArc[j]->head->address;
      }

      graph.getNodes().local[0] = ndInfo;
    }
  }



  dash::barrier();

  if( my_rank == 0 ) {
    printf( "\n\n NAS Parallel Benchmarks 3.3 -- DT Benchmark\n\n" );
    graphShow(dg,0);
    timer_clear(0);
    timer_start(0);
  }

  verified=ProcessNodes(graph);

  featnum=NUM_SAMPLES*fielddim;
  bytes_sent=featnum*dg->numArcs;
  bytes_sent/=1048576;
  if(my_rank==0) {
    timer_stop(0);
    tot_time=timer_read(0);

    c_print_results( dg->name,
                     CLASS,
                     featnum,
                     0,
                     0,
                     dg->numNodes,
                     0,
                     comm_size,
                     tot_time,
                     bytes_sent/tot_time,
                     "bytes transmitted",
                     verified,
                     NPBVERSION,
                     COMPILETIME,
                     MPICC,
                     CLINK,
                     CMPI_LIB,
                     CMPI_INC,
                     CFLAGS,
                     CLINKFLAGS );
  }
  /*
  printf("Before finalize, rank: %d\n", my_rank);
  if (my_rank == 9 || my_rank == 8) {
    int wait = 0;
    printf("Comparator processId is: %d\n", getpid());
    while(wait);
  }
  */
  dash::barrier();
  dash::finalize();
  return 0;
}
