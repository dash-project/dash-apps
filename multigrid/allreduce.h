#ifndef ALLREDUCE_H
#define ALLREDUCE_H

#include <libdash.h>

#ifdef SCOREP_USER_ENABLE
#include <scorep/SCOREP_User.h>
#define SCOREP_USER_FUNC() SCOREP_USER_REGION(__PRETTY_FUNCTION__, \
                                                SCOREP_USER_REGION_TYPE_FUNCTION)
#else
#define SCOREP_USER_FUNC()
#endif

/* The class is used for the async lazy residual computation.

Make this a template class with template arguments for the element type
and for the number of elements per unit.
*/

class Allreduce {

    /* distributed array for all local residual values of every unit.
    It is still going to be allocated in a way, that all elements are
    owned by unit 0.
    It can also be used by a subteam of the team, that created it. */
    dash::Array<double> centralized;

    /* truely distributed version where the global residula is copied to
    so that all units can access quickly. */
    dash::Array<double> distributed;

public:
    Allreduce( dash::Team& team ) :
      centralized( team.size(), dash::BLOCKCYCLIC(team.size()), team),
      distributed(team.size(), dash::BLOCKED, team) {
        reset(team);
    }

    /* can be used with a subteam of the team used in the constructor */
    void reset( dash::Team& team ) {
        SCOREP_USER_FUNC()
        /* really only unit 0 in the given team is doing something.
        std::fill is the really the correct algorithm */
        std::fill( centralized.lbegin(), centralized.lend(), std::numeric_limits<double>::max() );
        std::fill( distributed.lbegin(), distributed.lend(), std::numeric_limits<double>::max() );
    }


    /* can be used with a subteam of the team used in the constructor,
    does a barrier in the given team */
    void collect_and_spread( dash::Team& team ) {
        SCOREP_USER_FUNC()
        centralized.flush();
        team.barrier();

        /* now all new values are there, get the maximum on unit 0 */
        if ( 0 == team.myid() ) {
          /*for(const auto& elem : centralized.local)
            std::cout << elem << ",";
          std::cout << std::endl;*/
          distributed.local[0] =
            *std::max_element(centralized.lbegin(), centralized.lbegin() + team.size());

          for(auto i = 1; i < team.size(); ++i)
            distributed.async[i].set(distributed.local[0]);
        }
    }
#if 0
    /* broadcast the value from unit 0 to all units in the given team
    which might be a subteam of the team from constuction time */
    void asyncbroadcast( dash::Team& team ) {
        SCOREP_USER_FUNC()

        /* all but unit 0 fetch the value because unit 0 cannot know
        everybody else's local element in 'distributed' */
        if ( 0 == team.myid() ) {
          for(auto i = 1; i < team.size(); ++i)
            distributed.async[i].set(distributed.local[0]);
        }
    }

#endif
    void wait( dash::Team& team ) {
        SCOREP_USER_FUNC()

        distributed.async.flush();
        team.barrier();
    }

    /* send local residual to the correct place in unit 0's array,
    need to be followed by collect() eventually */
    void set( double* res, dash::Team& team ) {
        SCOREP_USER_FUNC()
      //std::cout << "TEST - " << team.myid() << " -> " << res << std::endl;
      if(team.myid() != 0)
        centralized.async[team.myid()].set(res);
      else
        centralized.local[0] = *res;
    }

    double get() const {
        SCOREP_USER_FUNC()
        return distributed.local[0];
    }
};

#endif /* ALLREDUCE_H */
