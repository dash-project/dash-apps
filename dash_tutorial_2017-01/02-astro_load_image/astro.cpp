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


int main( int argc, char* argv[] ) {

    dash::init(&argc, &argv);

    /* global vars */
    myid= dash::myid();
    team_size= dash::Team::All().size();

    dash::TeamSpec<2> teamspec(team_size, 1);
    teamspec.balance_extents();
    
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
        dash::DistributionSpec<2>( dash::BLOCKED, dash::BLOCKED ),
        dash::Team::All(), teamspec );


    /* *** part 2: load image strip by strip on unit 0, copy to distributed matrix from there, then show it *** */


    if ( 0 == myid ) {

        std::chrono::time_point<std::chrono::system_clock> start, end;

        uint32_t numstrips= TIFFNumberOfStrips(tif);
        tdata_t buf= _TIFFmalloc( TIFFStripSize(tif) );
        uint32_t line= 0;
        auto iter= matrix.begin();
        start = std::chrono::system_clock::now();
        for ( uint32_t strip = 0; strip < numstrips; strip++, line += rowsperstrip ) {

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

            // in the last iteration we can overwrite 'rowsperstrip'
            if ( line + rowsperstrip > h ) rowsperstrip= h - line;

            iter= dash::copy( rgb, rgb+w*rowsperstrip, iter );

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

    dash::finalize();
}
