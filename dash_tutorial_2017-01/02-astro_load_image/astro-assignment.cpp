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
    dash::TeamSpec<2> teamspec{};

    uint32_t rowsperstrip= 0;
    TIFF* tif= NULL;
    std::chrono::time_point<std::chrono::system_clock> start, end;

    /* *** Part 1: open TIFF input, find out image dimensions, create distributed DASH matrix *** */

    dash::Shared<ImageSize> imagesize_shared {};

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

        uint32_t width  = 0;
        uint32_t height = 0;

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width );
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height );
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsps );
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplespp );
        TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip );

        cout << "input file " << argv[1] << endl << "    image size: " << width << " x "
             << height << " pixels with " << samplespp << " channels x " << bitsps
             << " bits per pixel, " << rowsperstrip << " rows per strip" << endl;

        imagesize_shared.set( {height, width} );
    }

    imagesize_shared.barrier();

    ImageSize imagesize = imagesize_shared.get();
    cout << "unit " << myid << " thinks image is " << imagesize.width << " x " << imagesize.height << endl;

    auto distspec= dash::DistributionSpec<2>( dash::BLOCKED, dash::NONE );
    dash::NArray<RGB, 2> matrix( dash::SizeSpec<2>( imagesize.height, imagesize.width),
        distspec, dash::Team::All(), teamspec );

    /* Assignment: declare and allocate a distributed 2d DASH matrix
    container with width 'w' and height 'h' with the RGB basic data type.
    Select BLOCKED or NONE distribution in each dimension. In this case
    the TILED distribution is not recommended, though. */

    /* *** part 2: load image strip by strip on unit 0, copy to distributed matrix from there, then show it *** */

    matrix.barrier();

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
            std::for_each( rgb, rgb + imagesize.width * rowsperstrip, [](RGB& rgb) {
                std::swap<uint8_t>( rgb.r, rgb.b );
            } );


            /* Assignment: copy the next 'rowsperstrip' lines
            from pointer 'rgb' with width 'w' each to the target
            position in the distributed matrix. */

            if ( 0 == ( strip % 100 ) ) {
                cout << "    strip " << strip << "/" << numstrips << "\r" << flush;
            }
        }

        end= std::chrono::system_clock::now();
        cout << "read image in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
    }

    matrix.barrier();

    if ( 0 == myid ) {
        show_matrix( matrix, 1200, 1000 );
    }

    dash::finalize();
    return 0;
}
