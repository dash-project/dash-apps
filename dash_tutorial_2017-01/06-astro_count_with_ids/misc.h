#ifndef MISC_H
#define MISC_H

#include <unistd.h>
#include <iostream>
#include <set>
#include <cstddef>
#include <iomanip>
#include <cstdlib>

#include <libdash.h>

#ifdef WITHSDL
#include <SDL/SDL.h>
#endif


struct RGB {

    uint8_t r;
    uint8_t g;
    uint8_t b;

    RGB( uint8_t rr, uint8_t gg= 0, uint8_t bb= 0 ) : r(rr), g(gg), b(bb) {};
    RGB() = default;

    bool operator==( const RGB& other ) const {
        return ( r == other.r ) && ( g == other.g ) && ( b == other.b );
    };

    bool operator!=( const RGB& other ) const {
        return ( r != other.r ) || ( g != other.g ) || ( b != other.b );
    };

    uint32_t brightness() const { return (uint32_t) r + (uint32_t) g + (uint32_t) b; }
};


struct ImageSize {
    uint32_t height;
    uint32_t width;
};


/** show RGB matrix values graphically with SDL, show a window of w*h at most */
template<class MatrixT>
void show_matrix( MatrixT & matrix, uint32_t w= 400, uint32_t h= 300, uint32_t startx= 0, uint32_t starty= 0 ) {

#ifdef WITHSDL

    /* only first unit may do graphical output */
    if ( 0 != dash::myid() ) return;

    static SDL_Surface* pic= NULL;
    static uint32_t width= 0;
    static uint32_t height= 0;
    SDL_Event event;

    /* init SDL only on first time this is called */
    if ( NULL == pic ) {

        width  = ( w > 2000 ) ? 2000 : w;
        height = ( h > 2000 ) ? 2000 : h;

        SDL_Init(SDL_INIT_EVERYTHING);
        pic = SDL_SetVideoMode( width, height, 24, SDL_HWSURFACE );
        SDL_WM_SetCaption( "DASH RGB Matrix", "matrix" );
    }

    auto mw = matrix.extent(1);
    auto mh = matrix.extent(0);

    w = (mw < width) ? mw : width;
    h = (mh < height) ? mh : height;

    auto range = matrix.cols(startx,w).rows(starty,h);
    RGB* pixels = (RGB*) pic->pixels;

    /* copy only the selected range to the raw pointer of the SDL pic */
    std::copy( range.begin(), range.end(), pixels );

    SDL_Flip( pic );

    cout << endl << "Wait, please press any key ..." << flush;
    /* wait for key pressed before going on */
    do {
        SDL_Delay(50);
        SDL_PollEvent(&event);
    } while( event.type != SDL_QUIT && event.type != SDL_KEYDOWN );
    cout << " done" << endl << endl;

#endif /* WITHSDL */
}


template<class Iter>
void print_histogram( Iter first, Iter last ) {

    constexpr uint64_t HISTOGRAMWIDTH = 60;
    uint64_t max = *(std::max_element(first, last));
    uint64_t sum= 0;
    for ( auto it = first; it != last; ++it ) {

        uint32_t len = *it * HISTOGRAMWIDTH / max;
        cout << setw(len) << std::setfill('#') << "|" <<
            setw(12) << std::setfill(' ') << *it << endl;
        sum += *it;
    }
    cout << "sum is " << sum << endl;
}


constexpr dash::StencilSpec<2,4> stencil_spec({
    dash::Stencil<2>(-1, 0), dash::Stencil<2>( 1, 0),
    dash::Stencil<2>( 0,-1), dash::Stencil<2>( 0, 1)});

constexpr dash::CycleSpec<2> cycle_spec(
    dash::Cycle::NONE,
    dash::Cycle::NONE );


/* counts the objects in the local part of the image, it doesn't use a marker color but a separate array 
'indexes' where a object index per object is set. object index 0 means it is not part of an object. */
uint32_t countobjects( dash::NArray<RGB, 2>& matrix, dash::NArray<uint32_t,2>& indexes, uint32_t startindex, uint32_t limit ) {

    uint32_t found= 0;
    std::set< std::pair<uint32_t,uint32_t> > queue;

    /* step 1: go through all pixels and mark objects in the 'indexes' helper array */

    uint32_t lw= matrix.local.extent(1);
    uint32_t lh= matrix.local.extent(0);
    for ( uint32_t y= 0; y < lh ; y++ ){
        for ( uint32_t x= 0; x < lw ; x++ ){

            RGB& pixel= matrix.local[y][x];
            if ( ( 0 != indexes.local[y][x] ) || ( pixel.brightness() < limit ) ) {

                /* was already marked or not bright enough, do not count*/
                continue;
            }

            found++;

            queue.clear();
            queue.insert( { x, y } );

            uint32_t objid= startindex++;
            indexes.local[y][x]= objid;

            while ( ! queue.empty() ) {

                std::pair<uint32_t,uint32_t> next= *queue.begin();
                uint32_t x= next.first;
                uint32_t y= next.second;

                RGB& pixel= matrix.local[y][x];
                queue.erase( queue.begin() );

                if ( pixel.brightness() < limit ) continue;

                indexes.local[y][x]= objid;
                if ( ( 0 < x    ) && ( 0 == indexes.local[y][x-1] ) ) queue.insert( {x-1,y} );
                if ( ( x+1 < lw ) && ( 0 == indexes.local[y][x+1] ) ) queue.insert( {x+1,y} );
                if ( ( 0 < y    ) && ( 0 == indexes.local[y-1][x] ) ) queue.insert( {x,y-1} );
                if ( ( y+1 < lh ) && ( 0 == indexes.local[y+1][x] ) ) queue.insert( {x,y+1} );
            }
        }
    }

    /* step 2: create halo for helper array 'indexes', then do halo exchange */
    dash::HaloMatrixWrapper< dash::NArray<uint32_t,2>, dash::StencilSpec<2,4> > halo( indexes, stencil_spec, cycle_spec );
    halo.update();

    /* step 3: go along the local borders, if there is a local object (indexes[y][x] != 0) and the
    halo also has an object there ... if so, add the two object ids in a map (actually only keep 
    those with local object id < remote object id, otherwise one would count it twice). */

    std::array<long int,2> upperleft= indexes.pattern().global( {0,0} );
    std::array<long int,2> bottomright= indexes.pattern().global( {lh-1,lw-1} );

    std::set<std::pair<uint32_t,uint32_t>> touchmap;

    if ( upperleft[0] > 0 ) {

        /* there is an upper neighbor */
        for ( uint32_t x= 0; x < lw ; x++ ) {

            uint32_t localid= indexes.local[0][x];
            uint32_t remoteid= *halo.halo_element_at( {upperleft[0]-1,x+upperleft[1]} );

            if ( ( 0 != localid ) && ( 0 != remoteid ) && ( localid < remoteid ) ) {

                touchmap.insert( {localid, remoteid} );
            }
        }
    }

    if ( upperleft[1] > 0 ) {

        /* there is a left neighbor */
        cout << "  unit " << dash::myid() << " has left neigbor -- still to do" << endl;
    }

    if ( bottomright[0] < indexes.extent(0) -1 ) {

        /* there is a bottom neighbor */
        for ( uint32_t x= 0; x < lw ; x++ ) {

            uint32_t localid= indexes.local[lh-1][x];
            uint32_t remoteid= *halo.halo_element_at( {bottomright[0]+1,x+upperleft[1]} );

            if ( ( 0 != localid ) && ( 0 != remoteid ) && ( localid < remoteid ) ) {

                 touchmap.insert( {localid, remoteid} );
            }
        }
    }

    if ( bottomright[1] < indexes.extent(1) -1 ) {

        /* there is a right neighbor */
        cout << "  unit " << dash::myid() << " has right neigbor -- still to do" << endl;

    }

    /* step 4: subtract the number of entries in the map from the number of found objects */

    cout << "  unit " << dash::myid() << " touchmap " << touchmap.size() << endl;

    return found;
}


#endif /* MISC_H */
