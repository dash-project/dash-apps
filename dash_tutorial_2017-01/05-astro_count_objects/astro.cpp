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
    teamspec.balance_extents();

    uint32_t rowsperstrip= 0;
    TIFF* tif= NULL;
    std::chrono::time_point<std::chrono::system_clock> start, end;

    /* *** Part 1: open TIFF input, find out image dimensions, create distributed DASH matrix *** */

    /* Just a fixed size DASH array to distribute the image dimensions.
    Could also use dash::Shared<pair<...>> but we'd need more code for sure. */
    dash::Shared<ImageSize> image_size{};

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

        image_size.set( {height, width} );
    }

    image_size.barrier();

    auto distspec= dash::DistributionSpec<2>( dash::BLOCKED, dash::BLOCKED );
    ImageSize image_size_cached = image_size.get();
    dash::NArray<RGB, 2> matrix( dash::SizeSpec<2>( image_size_cached.height, image_size_cached.width),
        distspec, dash::Team::All(), teamspec );

    /* this is a workaround for incorrect local iterators. Local extents and the global iterator are fine, though. */
    RGB black(0,0,0);
    dash::fill( matrix.begin(), matrix.end(), black );

    /* *** part 2: load image strip by strip on unit 0, copy to distributed matrix from there, then show it *** */

    matrix.barrier();

    if ( 0 == myid ) {

        /*uint32_t numstrips= TIFFNumberOfStrips( tif );
        tdata_t buf= _TIFFmalloc( TIFFStripSize( tif ) );
        */
        auto iter= matrix.begin();

        uint32_t line= 0;
        start= std::chrono::system_clock::now();
        uint32_t imagelength;
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
        tdata_t buf = _TIFFmalloc(TIFFScanlineSize(tif));
        for( auto row = 0; row < imagelength; ++row)
        {
          TIFFReadScanline(tif, buf, row);
          RGB* rgb= (RGB*) buf;

          iter = std::copy(rgb, rgb + image_size_cached.width, iter);
          //dash::copy( rgb, rgb+w, iter );
        }
#if 0
        for ( uint32_t strip = 0; strip < numstrips; strip++, line += rowsperstrip) {

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

            iter = dash::copy( rgb, rgb+w*rowsperstrip, iter );

            if ( 0 == ( strip % 100 ) ) {
                cout << "    strip " << strip << "/" << numstrips << "\r" << flush;
            }
        }
#endif
        end= std::chrono::system_clock::now();
        cout << "read image in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
    }

    matrix.barrier();

    if ( 0 == myid ) {
        show_matrix( matrix, 1200, 700 );
    }

    matrix.barrier();

    /* *** part 3: compute historgramm in parallel *** */

    {
        // really need 64 bits for very large images
        constexpr uint64_t MAXKEY = 256*3;
        constexpr uint64_t BINS = 17;
        constexpr double   HIST_FACTOR = (double) BINS / MAXKEY;

        dash::Array<uint32_t> histogram( BINS * numunits, dash::BLOCKED );
        dash::fill( histogram.begin(), histogram.end(), 0 );

        start= std::chrono::system_clock::now();
        for ( const RGB& pixel : matrix.local )
                histogram.local[ pixel.brightness() * HIST_FACTOR ]++;

        histogram.barrier();

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
            print_histogram( histogram.lbegin(), histogram.lend() );
        }
    }

    /* from the brightness histogram we learned, that we should define all but the first two histogram bins
    as bright pixels */

    matrix.barrier();

    /* *** part 4: define a marker color and check that it is not appearing in the image yet *** */

    RGB marker(0,255,0);

    {
        start= std::chrono::system_clock::now();

        auto count = 0;
        for ( const RGB& pixel : matrix.local)
            if ( marker == pixel )
                count++;

        end = std::chrono::system_clock::now();
        cout << "    unit " << myid << " found marker color " << count << " times "
             << "in " << std::chrono::duration_cast<std::chrono::seconds> (end-start).count()
             << " seconds" << endl;
    }

    matrix.barrier();

    /* *** part 5: count bright objects in the image by finding a bright pixel, then flood-filling all
    bright neighbor pixels with marker color *** */

    {
        dash::Array<uint32_t> sums( numunits, dash::BLOCKED );
        start= std::chrono::system_clock::now();
        /*
        constexpr uint32_t limit= 256*3*2/17;
        uint32_t lw= matrix.local.extent(1);
        uint32_t lh= matrix.local.extent(0);
        uint32_t foundobjects= 0;
        for ( uint32_t y= 0; y < lh ; y++ ){
        for ( uint32_t x= 0; x < lw ; x++ ){

            foundobjects += checkobject( matrix.lbegin(), x, y, lw, lh, limit, marker );
        }
        }

        matrix.barrier();
        sums.local[0] = foundobjects;
        */

        sums.local[0]= check_objects_it(matrix, marker);
        end= std::chrono::system_clock::now();
        if ( 0 == myid ) {
            cout << "marked pixels in parallel in " << std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
        }

        /* combine the local number of objects found into the global result */

        sums.barrier();

        if ( 0 != myid ) {

            dash::transform<uint32_t>(
                sums.lbegin(), sums.lend(), // first source
                sums.begin(), // second source
                sums.begin(), // destination
                dash::plus<uint32_t>() );
        }

        sums.barrier();

        if ( 0 == myid ) {
            cout << "found " << sums.local[0] << " objects in total" << endl;
        }
    }

    matrix.barrier();

    if ( 0 == myid ) {
        show_matrix( matrix, 1200, 700 );
    }

    dash::finalize();
    return 0;
}
