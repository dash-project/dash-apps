#ifndef ALLREDUCE_H
#define ALLREDUCE_H


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
        centralized( team.size(), dash::BLOCKCYCLIC( team.size() ), team ),
        distributed( team.size(), dash::BLOCKED, team ) {
    }

    /* can be used with a subteam of the team used in the constructor */
    void reset( dash::Team& team ) {

        /* really only unit 0 in the given team is doing something.
        std::fill is the really the correct algorithm */
        std::fill( centralized.lbegin(), centralized.lend(), std::numeric_limits<double>::max() );
        std::fill( distributed.lbegin(), distributed.lend(), std::numeric_limits<double>::max() );
    }


    /* can be used with a subteam of the team used in the constructor, 
    does a barrier in the given team */
    void collect( dash::Team& team ) {

        centralized.async.flush();
        team.barrier();

        /* now all new values are there, get the maximum on unit 0 */
        if ( 0 == team.myid() ) {

            double* ptr= &centralized.local[0];
            double max= *ptr;
            size_t n= team.size();
            for ( size_t i= 1; i < n; ++i ) {

                max= ( max >= ptr[i] ) ? max : ptr[i];
            }
            distributed.local[0]= max;
        }
    }

    /* broadcast the value from unit 0 to all units in the given team
    which might be a subteam of the team from constuction time */
    void asyncbroadcast( dash::Team& team ) {

        /* all but unit 0 fetch the value because unit 0 cannot know
        everybody else's local element in 'distributed' */
        if ( 0 != team.myid() ) {

            distributed.local[0]= distributed.async[0];
        }
    }

    void waitbroadcast( dash::Team& team ) {

        distributed.async.flush();
        team.barrier();
    }

    /* send local residual to the correct place in unit 0's array,
    need to be followed by collect() eventually */
    void asyncset( double res, dash::Team& team ) {

        centralized.async[ team.myid() ]= res;
    }

    double get() const {

        return distributed.local[0];
    }
};

#endif /* ALLREDUCE_H */
