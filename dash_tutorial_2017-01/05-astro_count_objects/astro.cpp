#include <unistd.h>
#include <iostream>
#include <set>
#include <cstddef>
#include <iomanip>
#include <chrono>
#include <cstdlib>

#include <libdash.h>

//#define WITHSDL
#ifdef WITHSDL
#include <SDL/SDL.h>
#endif

extern "C" {
    #include "tiffio.h"
}

using std::cout;
using std::endl;
using std::setw;
using std::flush;


/* gloal var for simplicity */
dart_unit_t myid;
size_t team_size;


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
    if ( 0 != myid ) return;

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

    auto range = matrix.rows(startx,w).cols(starty,h);
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


int main( int argc, char* argv[] ) {

    dash::init(&argc, &argv);

    /* global vars */
    myid= dash::myid();
    team_size= dash::Team::All().size();

    dash::TeamSpec<2> teamspec(team_size, 1);

    uint32_t w= 0;
    uint32_t h= 0;
    uint32_t rowsperstrip= 0;

    TIFF* tif= NULL;


    /* *** part 1: open TIFF input, find out image dimensions, create distributed DASH matrix *** */


    /* just a fixed size DASH array to distribute the image dimensions.
    Could also use dash::Shared<pair<...>> but we'd need more code for sure. */
    dash::Array<uint32_t> arr(2);

    if ( 0 == myid ) {

        if ( 1 >= argc ) {

            printf( "no file given\n" );
            dash::finalize();
            return EXIT_FAILURE;

        } else {

            tif = TIFFOpen( argv[1], "r" );
            if ( tif ) {

                uint32_t bitsps= 0;
                uint32_t samplespp= 0;

                TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
                TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
                TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsps );
                TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplespp );
                TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip );

                cout << "input file " << argv[1] << endl <<
                    "    image size: " << w << " x " << h << " pixels "
                    "with " << samplespp << " channels x " << bitsps << " bits per pixel, " <<
                    rowsperstrip << " rows per strip" << endl;
            }
        }

        arr[0]= w;
        arr[1]= h;
    }

    dash::barrier();

    if ( 0 != myid ) {

        w= arr[0];
        h= arr[1];
    }

    // cout << "unit " << myid << " thinks the image is " << w << "*" << h << endl;

    dash::barrier();

    dash::Matrix<RGB, 2> matrix( dash::SizeSpec<2>( h, w ),
        dash::DistributionSpec<2>( dash::BLOCKED, dash::NONE ),
        dash::Team::All(), teamspec );


    /* *** part 2: load image strip by strip on unit 0, copy to distributed matrix from there, then show it *** */


    if ( 0 == myid ) {

        std::chrono::time_point<std::chrono::system_clock> start, end;

        uint32_t numstrips= TIFFNumberOfStrips(tif);
        tdata_t buf= _TIFFmalloc( TIFFStripSize(tif) );
        uint32_t line= 0;
        start = std::chrono::system_clock::now();
        for ( uint32_t strip = 0; strip < numstrips; strip++ ) {

            TIFFReadEncodedStrip( tif, strip, buf, (tsize_t) -1);
            RGB* rgb= (RGB*) buf;

            /* it turns out that libtiff and SDL have different oponions about
            the RGB byte order, fix it by swaping red and blue. Do it here where
            the line was freshly read in. This is qicker than an extra sequential or parallel
            swap operation over the matrix when the lines have been out of cache already.
            Then again, we don't care about the colors for the rest of the code,
            and we can leave this out entirely. Histogram and counting objects won't
            produce any difference*/
            std::for_each( rgb, rgb+w*rowsperstrip, [](RGB& rgb){
                std::swap<uint8_t>( rgb.r, rgb.b );
            } );

            for ( uint32_t l= 0; ( l < rowsperstrip ) && ( line < h ) ; l++, line++, rgb += w ) {

                auto range = matrix.cols( line, 1 );
                dash::copy( rgb, rgb+w, range.begin() );
            }

            if ( 0 == ( strip % 100 ) ) {
                cout << "    strip " << strip << "/" << numstrips << "\r" << flush;
            }
        }
        end = std::chrono::system_clock::now();
        cout << "read image in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
    }

    dash::barrier();

    if ( 0 == myid ) {

        show_matrix( matrix, 1600, 1200 );
    }

    dash::barrier();


    /* *** part 3: compute historgramm in parallel *** */


    {
        // really need 64 bits for very large images
        const uint64_t MAXKEY= 255*3;
        const uint64_t BINS= 17;

        dash::Array<uint32_t> histogram( BINS * team_size, dash::BLOCKED );
        dash::fill( histogram.begin(), histogram.end(), (uint32_t) 0 );

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        for ( auto it= matrix.lbegin(); it != matrix.lend(); ++it ) {

                histogram.local[ it->brightness() * BINS / MAXKEY ]++;
        }

        /* is this barrier needed? The example 'bench.04.histo-tf' doesn't have it */
        dash::barrier();

        if ( 0 != myid ) {

            // is this using atomic accesses?
            dash::transform<uint32_t>(
                histogram.lbegin(), histogram.lend(), // first source
                histogram.begin(), // second source
                histogram.begin(), // destination
                dash::plus<uint32_t>() );
        }

        dash::barrier();

        end = std::chrono::system_clock::now();

        if ( 0 == myid ) {

            cout << "computed parallel historgram in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl << endl;
            print_histogram<uint32_t*>( histogram.lbegin(), histogram.lend() );
        }
    }

    /* from the brightness histogram we learned, that we should define all but the first two histogram bins
    as bright pixels */
    const uint32_t limit= 255*3*2/17;


    dash::barrier();


    /* *** part 4: define a marker color and check that it is not appearing in the image yet *** */


    RGB marker(0,255,0);
    {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        uint64_t count= 0;
        for ( auto it= matrix.lbegin(); it != matrix.lend(); ++it ) {

            if ( marker ==  *it ) count++;
        }

        end = std::chrono::system_clock::now();
        if ( 0 == myid ) cout << "computed parallel count of pixel value " <<
            (uint32_t) marker.r << ":" <<
            (uint32_t) marker.g << ":" <<
            (uint32_t) marker.b << ":" <<
            "in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
        dash::barrier();
        cout << "    unit " << myid << " found it " << count << " times" << endl;
    }

    dash::barrier();


    /* *** part 5: count bright objects in the image by finding a bright pixel, then flood-filling all
    bright neighbor pixels with marker color *** */


    {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        uint32_t lw= matrix.local.extent(1);
        uint32_t lh= matrix.local.extent(0);

        /* warning */
        if ( lh * team_size != h ) {

            cout << "            unit " << myid << " of " << team_size <<
            ": local height " << lh << " * teamsize " << team_size <<
            " != total height " << h << endl <<
            "           ??? shouldn't unit's local heights differ, "
            "if total height is not a multiple of the number of units?" << endl;
        }

        uint32_t foundobjects= 0;

        for ( uint32_t y= 0; y < lh ; y++ ){
        for ( uint32_t x= 0; x < lw ; x++ ){

            foundobjects += checkobject( matrix.lbegin(),
                x, y, lw, lh,
                limit, marker );
        }
        }


        /* now we should correct for objects countet twice because they are at the
        border between units. This is left as an exercise ... */


        dash::barrier();
        end = std::chrono::system_clock::now();
        if ( 0 == myid ) {
            cout << "marked pixels in parallel in " << std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
        }

        /* combine the local number of objects found into the global result */

        dash::Array<uint32_t> sums( team_size, dash::BLOCKED );
        sums.local[0]= foundobjects;

        /* is this barrier needed? The example 'bench.04.histo-tf' doesn't have it */
        dash::barrier();

        if ( 0 != myid ) {

            // is this using atomic accesses?
            dash::transform<uint32_t>(
                sums.lbegin(), sums.lend(), // first source
                sums.begin(), // second source
                sums.begin(), // destination
                dash::plus<uint32_t>() );
        }

        dash::barrier();

        if ( 0 == myid ) {

            cout << "found " << sums.local[0] << " objects in total" << endl;

        }
    }

    dash::barrier();


    /* show again an part of the large image */

    if ( 0 == myid ) {

        show_matrix( matrix, 1600, 1200 );
    }

    dash::finalize();
}
