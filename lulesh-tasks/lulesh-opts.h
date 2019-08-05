#ifndef LULESH_OPTS_H_INCLUDED
#define LULESH_OPTS_H_INCLUDED

#include <iostream>
#include "lulesh.h"

class CmdLineOpts
{
private:
  Int_t m_its;      // -i
  Int_t m_nx;       // -s
  Int_t m_numReg;   // -r
  Int_t m_numFiles; // -f
  Int_t m_showProg; // -p
  Int_t m_quiet;    // -q
  Int_t m_viz;      // -v
  Int_t m_cost;     // -c
  Int_t m_balance;  // -b

  // DASH additions:
  Int_t m_numRanks;
  Int_t m_myRank;

  Int_t m_px; // processes in (x, y, z) dimension
  Int_t m_py;
  Int_t m_pz;

  bool m_valid;

  Int_t m_chunksize;

public:
  CmdLineOpts(Int_t numRanks, Int_t myRank)
  {
    // Set defaults that can be overridden by command line opts
    m_its      = 9999999;
    m_nx       = 3; // 30
    m_numReg   = 11;
    m_numFiles = (int)(numRanks+10)/9;
    m_showProg = 0;
    m_quiet    = 0;
    m_viz      = 0;
    m_balance  = 1;
    m_cost     = 1;

    m_numRanks = numRanks;
    m_myRank   = myRank;
    m_px       = 0;
    m_py       = 0;
    m_pz       = 0;
    m_valid    = true;

    m_chunksize = 6;
  }

  void parseCommandLineOptions(int argc, char *argv[]);

  void printBanner(std::ostream& os);
  void printHelp(std::ostream& os);

  Int_t numRanks() const { return m_numRanks; }
  Int_t quiet()    const { return m_quiet; }
  Int_t nx()       const { return m_nx; }
  Int_t px()       const { return m_px; }
  Int_t py()       const { return m_py; }
  Int_t pz()       const { return m_pz; }
  Int_t showProg() const { return m_showProg; }
  Int_t numReg()   const { return m_numReg; }
  Int_t balance()  const { return m_balance; }
  Int_t cost()     const { return m_cost; }

  Int_t chunksize() const { return m_chunksize; }

  Int_t its()      const { return m_its; }
  bool valid()     const { return m_valid; }
};


static void PrintCommandLineOptions(char *execname, int myRank);

#endif /* LULESH_OPTS_H_INCLUDED */
