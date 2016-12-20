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

        dash::Array<uint32_t> histogram( BINS * numunits, dash::BLOCKED );
        dash::fill( histogram.begin(), histogram.end(), (uint32_t) 0 );

        start= std::chrono::system_clock::now();
        for ( auto it= matrix.lbegin(); it != matrix.lend(); ++it ) {

                histogram.local[ it->brightness() * BINS / MAXKEY ]++;
        }

        if ( 0 != myid ) {

            dash::transform<uint32_t>(
                histogram.lbegin(), histogram.lend(), // first source
                histogram.begin(), // second source
                histogram.begin(), // destination
                dash::plus<uint32_t>() );
        }

        histogram.barrier();
        end= std::chrono::system_clock::now();

        if ( 0 == myid ) {
            cout << "computed parallel historgram in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl << endl;
            print_histogram<uint32_t*>( histogram.lbegin(), histogram.lend() );
        }
    }

    /* from the brightness histogram we learned, that we should define all but the first two histogram bins
    as bright pixels */
    const uint32_t limit= 255*3*2/17;

    matrix.barrier();

    /* *** part 4: define a marker color and check that it is not appearing in the image yet *** */

    /* Assignment: define a marker color, which is not appearing
    in the image yet. Test with the following code if it is really not used so far. */
    RGB marker(,,);

    {
        start= std::chrono::system_clock::now();

        uint64_t count= 0;
        /* Assignment: Go through all pixels in the local part of the distributed matrix
        and count how often the marker value appears.  */

        end= std::chrono::system_clock::now();
        cout << "    unit " << myid << " found marker color " << count << " times " <<
            "in " << std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
    }

    dash::finalize();
    return 0;
}
