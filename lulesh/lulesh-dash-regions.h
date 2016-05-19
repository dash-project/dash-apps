#ifndef LULESH_DASH_REGIONS_H_INCLUDED
#define LULESH_DASH_REGIONS_H_INCLUDED

struct RegionIndexSet
{
private:
  // Region information
  Int_t      m_numReg;
  Int_t      m_cost;          // imbalance cost
  Index_t   *m_regElemSize;   // Size of region sets
  Index_t   *m_regNumList;    // Region number per domain element
  Index_t  **m_regElemlist;   // region indexset

public:
  RegionIndexSet(Int_t nr, Int_t balance, Index_t numElem);
  ~RegionIndexSet();

  Index_t&  numReg()             { return m_numReg ; }
  Int_t&    cost()               { return m_cost ; }

  //
  // Element-centered
  //
  Index_t&  regElemSize(Index_t idx) { return m_regElemSize[idx] ; }
  Index_t&  regNumList(Index_t idx)  { return m_regNumList[idx] ; }
  Index_t*  regNumList()             { return &m_regNumList[0] ; }
  Index_t*  regElemlist(Int_t r)     { return m_regElemlist[r] ; }
  Index_t&  regElemlist(Int_t r, Index_t idx) { return m_regElemlist[r][idx] ; }
};

#endif /* LULESH_DASH_REGIONS_H_INCLUDED */
