#ifndef MISC_H
#define MISC_H

#include <unistd.h>
#include <iostream>
#include <set>
#include <cstddef>
#include <iomanip>
#include <cstdlib>

using std::cout;
using std::endl;
using std::setw;
using std::flush;

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

    if ( NULL == pic ) {

        /* init SDL only on first time this is called */

        if ( w > 2000 ) w= 2000;
        if ( h > 2000 ) h= 2000;

        SDL_Init(SDL_INIT_EVERYTHING);
        pic = SDL_SetVideoMode( w, h, 24, SDL_HWSURFACE );
        SDL_WM_SetCaption( "DASH RGB Matrix", "matrix" );

        /* allow to set width and height only the first time this is called */
        width= w;
        height= h;
    }

    auto mw = matrix.extent(1);
    auto mh = matrix.extent(0);

    w= width;
    h= height;

    if ( mw < w ) w= mw;
    if ( mh < h ) h= mh;

    auto range = matrix.cols(startx,w).rows(starty,h);
    RGB* pixels = (RGB*) pic->pixels;

    /* copy only the selected range to the raw pointer of the SDL pic */
    /* using dash::copy here causes a strange image or a crash! */
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


template<class LocalHistIter>
void print_histogram( LocalHistIter first, LocalHistIter last ) {

    uint64_t max= 1;
    for ( auto it= first; it != last; ++it ) {

        if ( *it > max ) max= *it;
    }

#define HISTOGRAMWIDTH 60

    for ( auto it= first; it != last; ++it ) {

        uint32_t len= ( (uint64_t) HISTOGRAMWIDTH * (uint64_t) *it ) / max;

        cout << setw(len) << std::setfill('#') << "|" <<
            setw(12) << std::setfill(' ') << *it << endl;
    }
    cout << endl;
}


/* The checkobject function checks, if a pixel is brighter than the given limit.
If so, mark it with the marker color and return 1. Also mark all adjacent bright pixels
so that they will not counted again. This is done with a flood-fill algorithm.
Note: Usually, one would implement the flood-fill recursively but this breaks for very
large images! While it would work fine and with much less lines of code for small
images, for very large images the  call stack gets very deep, causing stray segfaults. */
uint32_t checkobject( RGB* ptr,
        uint32_t x, uint32_t y,
        uint32_t w, uint32_t h,
        uint32_t limit, RGB marker ) {

    std::set< std::pair<uint32_t,uint32_t> > queue;
    queue.insert( { x, y } );

    RGB* pixel= ptr + y*w+ x;
    if ( *pixel == marker ) return 0;
    if ( pixel->brightness() < limit ) return 0;

    *pixel= marker;

    while ( ! queue.empty() ) {

        std::pair<uint32_t,uint32_t> next= *queue.begin();
        uint32_t x= next.first;
        uint32_t y= next.second;
        RGB* pixel= ptr + y*w+ x;
        queue.erase( queue.begin() );

        if ( pixel->brightness() < limit ) continue;

        *pixel= marker;
        if ( 0 < x )   if ( *(pixel-1) != marker ) queue.insert( {x-1,y} );
        if ( x+1 < w ) if ( *(pixel+1) != marker ) queue.insert( {x+1,y} );
        if ( 0 < y )   if ( *(pixel-w) != marker ) queue.insert( {x,y-1} );
        if ( y+1 < h ) if ( *(pixel+w) != marker ) queue.insert( {x,y+1} );
    }

    return 1;
}

#endif /* MISC_H */
