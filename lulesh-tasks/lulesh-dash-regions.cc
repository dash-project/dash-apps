
#include <mpi.h>
#include "lulesh.h"
#include "lulesh-dash.h"
#include "lulesh-dash-regions.h"

RegionIndexSet::RegionIndexSet(Int_t numReg, Int_t cost,
			       Int_t balance, Index_t numElem)
{
  Index_t myRank = 0;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
#endif

#if USE_DASH
  myRank = dash::myid();
#endif


  srand(myRank);

  m_numReg      = numReg;
  m_cost        = cost;
  m_regElemSize = new Index_t[numReg];
  m_regElemlist = new Index_t*[numReg];
  m_regNumList  = new Index_t[numElem] ; // material indexset

  Index_t nextIndex = 0;

  //if we only have one region just fill it
  // Fill out the regNumList with material numbers, which are always
  // the region index plus one
  if(numReg == 1) {
    while (nextIndex < numElem) {
      this->regNumList(nextIndex) = 1;
      nextIndex++;
    }
    regElemSize(0) = 0;
  }
  //If we have more than one region distribute the elements.
  else {
    Int_t regionNum;
    Int_t regionVar;
    Int_t lastReg = -1;
    Int_t binSize;
    Index_t elements;
    Index_t runto = 0;
    Int_t costDenominator = 0;
    Int_t* regBinEnd = new Int_t[numReg];
    //Determine the relative weights of all the regions.  This is
    //based off the -b flag.  Balance is the value passed into b.
    for (Index_t i=0 ; i<numReg ; ++i) {
      regElemSize(i) = 0;
      costDenominator += pow((i+1), balance);  //Total sum of all
					       //regions weights
      regBinEnd[i] = costDenominator;  //Chance of hitting a given
				       //region is (regBinEnd[i] -
				       //regBinEdn[i-1])/costDenominator
    }

    //Until all elements are assigned
    while (nextIndex < numElem) {
      //pick the region
      regionVar = rand() % costDenominator;
      Index_t i = 0;
      while(regionVar >= regBinEnd[i])
	i++;

      //rotate the regions based on MPI rank.  Rotation is Rank %
      //NumRegions this makes each domain have a different region with
      //the highest representation
      regionNum = ((i + myRank) % numReg) + 1;
      // make sure we don't pick the same region twice in a row
      while(regionNum == lastReg) {
	regionVar = rand() % costDenominator;
	i = 0;
	while(regionVar >= regBinEnd[i])
	  i++;
	regionNum = ((i + myRank) % numReg) + 1;
      }
      //Pick the bin size of the region and determine the number of
      //elements.
      binSize = rand() % 1000;
      if(binSize < 773) {
	elements = rand() % 15 + 1;
      }
      else if(binSize < 937) {
	elements = rand() % 16 + 16;
      }
      else if(binSize < 970) {
	elements = rand() % 32 + 32;
      }
      else if(binSize < 974) {
	elements = rand() % 64 + 64;
      }
      else if(binSize < 978) {
	elements = rand() % 128 + 128;
      }
      else if(binSize < 981) {
	elements = rand() % 256 + 256;
      }
      else
	elements = rand() % 1537 + 512;
      runto = elements + nextIndex;

      //Store the elements.  If we hit the end before we run out of elements then just stop.
      while (nextIndex < runto && nextIndex < numElem) {
	this->regNumList(nextIndex) = regionNum;
	nextIndex++;
      }
      lastReg = regionNum;
    }
  }

  // Convert regNumList to region index sets
  // First, count size of each region
  for (Index_t i=0 ; i<numElem ; ++i) {
    int r = this->regNumList(i)-1; // region index == regnum-1
    regElemSize(r)++;
  }

  // Second, allocate each region index set
  for (Index_t i=0 ; i<numReg ; ++i) {
    m_regElemlist[i] = new Index_t[regElemSize(i)];
    regElemSize(i) = 0;
  }

  // Third, fill index sets
  for (Index_t i=0 ; i<numElem ; ++i) {
    Index_t r = regNumList(i)-1;       // region index == regnum-1
    Index_t regndx = regElemSize(r)++; // Note increment
    regElemlist(r,regndx) = i;
  }
}

RegionIndexSet::~RegionIndexSet()
{
  for (Index_t i=0 ; i<numReg() ; ++i) {
    delete[](m_regElemlist[i]);
  }
  delete[](m_regElemSize);
  delete[](m_regElemlist);
}
