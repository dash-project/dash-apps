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



#endif /* MISC_H */
