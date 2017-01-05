#include <unistd.h>
#include <iostream>
#include <cstddef>
#include <iomanip>
#include <chrono>
#include <cstdlib>

using std::cout;
using std::endl;
using std::setw;
using std::flush;

#include <libdash.h>
extern "C" {
    #include <tiffio.h>
}
#include "misc.h"


int main( int argc, char* argv[] ) {

    dash::init( &argc, &argv );

    dart_unit_t myid= dash::myid();
    size_t numunits= dash::Team::All().size();
    dash::TeamSpec<2> teamspec( numunits, 1 );
    teamspec.balance_extents();

    uint32_t w= 0;
    uint32_t h= 0;
    uint32_t rowsperstrip= 0;
    TIFF* tif= NULL;
    std::chrono::time_point<std::chrono::system_clock> start, end;

    /* *** Part 1: open TIFF input, find out image dimensions, create distributed DASH matrix *** */

    /* Just a fixed size DASH array to distribute the image dimensions.
    Could also use dash::Shared<pair<...>> but we'd need more code for sure. */
    dash::Array<uint32_t> array(2);

    if ( 0 == myid ) {

        if ( 1 >= argc ) {

            cout << "no file given (x_x)" << endl;
            dash::finalize();
            return EXIT_FAILURE;
        }

        tif = TIFFOpen( argv[1], "r" );
        if ( NULL == tif ) {

            cout << "not a TIFF file (x_x)" << endl;
            dash::finalize();
            return EXIT_FAILURE;
        }

        uint32_t bitsps= 0;
        uint32_t samplespp= 0;

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w );
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h );
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsps );
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplespp );
        TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip );

        cout << "input file " << argv[1] << endl << "    image size: " << w << " x " << h << " pixels "
            "with " << samplespp << " channels x " << bitsps << " bits per pixel, " <<
            rowsperstrip << " rows per strip" << endl;

        array[0]= w;
        array[1]= h;
    }

    array.barrier();

    if ( 0 != myid ) {

        w= array[0];
        h= array[1];
    }

    dash::Matrix<RGB, 2> matrix( dash::SizeSpec<2>( h, w ),
        dash::DistributionSpec<2>( dash::BLOCKED, dash::BLOCKED ),
        dash::Team::All(), teamspec );

    /* *** part 2: load image strip by strip on unit 0, copy to distributed matrix from there, then show it *** */

    if ( 0 == myid ) {

        uint32_t numstrips= TIFFNumberOfStrips( tif );
        tdata_t buf= _TIFFmalloc( TIFFStripSize( tif ) );
        auto iter= matrix.begin();
        uint32_t line= 0;
        start= std::chrono::system_clock::now();
        for ( uint32_t strip = 0; strip < numstrips; strip++, line += rowsperstrip ) {

            TIFFReadEncodedStrip( tif, strip, buf, (tsize_t) -1 );
            RGB* rgb= (RGB*) buf;

            /* it turns out that libtiff and SDL have different oponions about
            the RGB byte order, fix it by swaping red and blue. Do it here where
            the line was freshly read in. This is qicker than an extra sequential or parallel
            swap operation over the matrix when the lines have been out of cache already.
            Then again, we don't care about the colors for the rest of the code,
            and we can leave this out entirely. Histogram and counting objects won't
            produce any difference*/
            std::for_each( rgb, rgb+w*rowsperstrip, [](RGB& rgb) {
                std::swap<uint8_t>( rgb.r, rgb.b );
            } );

            // in the last iteration we can overwrite 'rowsperstrip'
            if ( line + rowsperstrip > h ) rowsperstrip= h - line;

            iter= dash::copy( rgb, rgb+w*rowsperstrip, iter );

            if ( 0 == ( strip % 100 ) ) {
                cout << "    strip " << strip << "/" << numstrips << "\r" << flush;
            }
        }
        end= std::chrono::system_clock::now();
        cout << "read image in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
    }

    matrix.barrier();

    if ( 0 == myid ) {
        show_matrix( matrix, 1600, 1200 );
    }

    matrix.barrier();

    /* *** part 3: compute historgramm in parallel *** */

    {
        // really need 64 bits for very large images
        const uint64_t MAXKEY= 255*3;
        const uint64_t BINS= 17;

        /* Assignment: Create a distributed DASH array for the histogram
        of size ( 'BINS' x team size ) which is block distributed
        over all units. Fill it with '0' using a DASH algorithm. */

        start= std::chrono::system_clock::now();
        /* Assignment: Iterate over all local pixels in the
        distributed matrix using the local iterator. For every pixel
        determine the brightness bin in the histogram of 'BINS' bins.
        Then increment the corresponding histogram bin in the local
        part of the histogram by 1. */

        /* Assignment: Add a barrier for the histogram data structure */

        if ( 0 != myid ) {

            /* Assignment: For all units except unit 0, add their histogram
            bins to the values of the histogram of unit 0. Try using a DASH
            algorithm for that. */
        }

        /* Assignment: Add a barrier for the histogram data structure */

        end= std::chrono::system_clock::now();

        if ( 0 == myid ) {
            cout << "computed parallel historgram in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl << endl;

            /* Assignment: Use the given function 'print_histogram' to print out the global histogram. */

        }
    }

    /* from the brightness histogram we learned, that we should define all but the first two histogram bins
    as bright pixels */

    /* Assignment: put here the limit where the pixel color is considered
    part of a bright object instead of the dark background. Look at the
    histogram you produced to find out a good limit */

    const uint32_t limit= 0;

    dash::finalize();
    return 0;
}
