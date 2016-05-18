#ifndef LULESH_DASH_PARAMS_H_INCLUDED
#define LULESH_DASH_PARAMS_H_INCLUDED

//
// factored out all scalar parameters and constants into a separate
// struct for better code readability
//
struct Parameters
{
  // Cutoffs (treat as constants)
  const Real_t m_e_cut;             // energy tolerance
  const Real_t m_p_cut;             // pressure tolerance
  const Real_t m_q_cut;             // q tolerance
  const Real_t m_v_cut;             // relative volume tolerance
  const Real_t m_u_cut;             // velocity tolerance

  // Other constants (usually setable, but hardcoded in this proxy app)
  const Real_t m_hgcoef;            // hourglass control
  const Real_t m_ss4o3;
  const Real_t m_qstop;             // excessive q indicator
  const Real_t m_monoq_max_slope;
  const Real_t m_monoq_limiter_mult;
  const Real_t m_qlc_monoq;         // linear term coef for q
  const Real_t m_qqc_monoq;         // quadratic term coef for q
  const Real_t m_qqc;
  const Real_t m_eosvmax;
  const Real_t m_eosvmin;
  const Real_t m_pmin;              // pressure floor
  const Real_t m_emin;              // energy floor
  const Real_t m_dvovmax;           // maximum allowable volume change
  const Real_t m_refdens;           // reference density

  // Variables to keep track of timestep, simulation time, and cycle
  Real_t  m_dtfixed;           // fixed time increment
  Real_t  m_stoptime;          // end time for simulation

  Real_t  m_deltatimemultlb;
  Real_t  m_deltatimemultub;
  Real_t  m_dtcourant;         // courant constraint
  Real_t  m_dthydro;           // volume change constraint
  Real_t  m_dtmax;             // maximum allowable time increment
  Int_t   m_cycle;             // iteration count for simulation

  Real_t  m_time;              // current time
  Real_t  m_deltatime;         // variable time increment

public:
  Parameters() :
    m_e_cut(              Real_t(1.0e-7)          ),
    m_p_cut(              Real_t(1.0e-7)          ),
    m_q_cut(              Real_t(1.0e-7)          ),
    m_v_cut(              Real_t(1.0e-10)         ),
    m_u_cut(              Real_t(1.0e-7)          ),
    m_hgcoef(             Real_t(3.0)             ),
    m_ss4o3(              Real_t(4.0)/Real_t(3.0) ),
    m_qstop(              Real_t(1.0e+12)         ),
    m_monoq_max_slope(    Real_t(1.0)             ),
    m_monoq_limiter_mult( Real_t(2.0)             ),
    m_qlc_monoq(          Real_t(0.5)             ),
    m_qqc_monoq(          Real_t(2.0)/Real_t(3.0) ),
    m_qqc(                Real_t(2.0)             ),
    m_eosvmax(            Real_t(1.0e+9)          ),
    m_eosvmin(            Real_t(1.0e-9)          ),
    m_pmin(               Real_t(0.)              ),
    m_emin(               Real_t(-1.0e+15)        ),
    m_dvovmax(            Real_t(0.1)             ),
    m_refdens(            Real_t(1.0)             )
  {
    // Setup defaults

    // These can be changed (requires recompile) if you want to run
    // with a fixed timestep, or to a different end time, but it's
    // probably easier/better to just run a fixed number of timesteps
    // using the -i flag in 2.x

    m_dtfixed   = Real_t(-1.0e-6);  // Negative means use courant condition
    m_stoptime  = Real_t( 1.0e-2);  // *Real_t(edgeElems*tp/45.0);

    // Initial conditions
    m_deltatimemultlb  = Real_t(1.1);
    m_deltatimemultub  = Real_t(1.2);
    m_dtcourant        = Real_t(1.0e+20);
    m_dthydro          = Real_t(1.0e+20);
    m_dtmax            = Real_t(1.0e-2);
    m_cycle            = Int_t(0);

    m_time             = Real_t(0.0);
    m_deltatime        = Real_t(0.0);
  }
};

#endif /* LULESH_DASH_PARAMS_H_INCLUDED */
