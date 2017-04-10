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

    char hostname[100];
    double time_creatematrix;
    double time_readtiff;
    double time_histogram;
    double time_searchmarkercolor;
    double time_checkobjects;


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

    start= std::chrono::system_clock::now();
    auto distspec= dash::DistributionSpec<2>( dash::BLOCKED, dash::NONE );
    dash::NArray<RGB, 2> matrix( dash::SizeSpec<2>( imagesize.height, imagesize.width),
        distspec, dash::Team::All(), teamspec );

    end= std::chrono::system_clock::now();
    time_creatematrix= 0.001*std::chrono::duration_cast<std::chrono::milliseconds> (end-start).count() ;

    matrix.barrier();

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
            std::for_each( rgb, rgb + imagesize.width * rowsperstrip, [](RGB& rgb) {
                std::swap<uint8_t>( rgb.r, rgb.b );
            } );

            for ( uint32_t l= 0; ( l < rowsperstrip ) && (line+l < imagesize.height); l++ ) {

                iter = dash::copy( rgb, rgb+imagesize.width, iter );
                rgb += imagesize.width;
            }

            if ( 0 == ( strip % 100 ) ) {
                cout << "    strip " << strip << "/" << numstrips << "\r" << flush;
            }
        }

        end= std::chrono::system_clock::now();
        cout << "read image in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
        time_readtiff= 0.001*std::chrono::duration_cast<std::chrono::milliseconds> (end-start).count() ;
    }

    matrix.barrier();

    /* *** part 3: compute historgramm in parallel *** */

    {
        // really need 64 bits for very large images
        constexpr uint64_t MAXKEY = 256*3;
        constexpr uint64_t BINS = 17;

        dash::Array<uint32_t> histogram( BINS * numunits, dash::BLOCKED );
        dash::fill( histogram.begin(), histogram.end(), 0 );

        start= std::chrono::system_clock::now();
        for ( auto it= matrix.lbegin(); it != matrix.lend(); ++it ) {
                histogram.local[ it->brightness() * BINS / MAXKEY ]++;
        }

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
            cout << "computed parallel histogram in "<< std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl << endl;
            print_histogram( histogram.lbegin(), histogram.lend() );
        }
        time_histogram= 0.001 * std::chrono::duration_cast<std::chrono::milliseconds> (end-start).count();
    }

    /* from the brightness histogram we learned, that we should define all but the first two histogram bins
    as bright pixels */
    const uint32_t limit= 256*3*2/17;

    matrix.barrier();

    /* *** part 4: define a marker color and check that it is not appearing in the image yet *** */

    RGB marker(0,255,0);

    {
        start= std::chrono::system_clock::now();

        uint64_t count= 0;
        uint64_t visitedpixels= 0;
        uint64_t brightnesssum= 0;
        for ( auto it= matrix.lbegin(); it != matrix.lend(); ++it ) {

            if ( marker == *it ) count++;
	    visitedpixels++;
	    brightnesssum += it->brightness();
        }

        end= std::chrono::system_clock::now();
        cout << "    unit " << myid << " found marker color " << count << " times "
             << "in " << std::chrono::duration_cast<std::chrono::seconds> (end-start).count()
             << " seconds" << ", visited " << visitedpixels << ", sum of brightness "
             << brightnesssum << endl;

        time_searchmarkercolor= 0.001 * std::chrono::duration_cast<std::chrono::milliseconds> (end-start).count();
    }

    matrix.barrier();

    /* *** part 5: count bright objects in the image by finding a bright pixel, then flood-filling all
    bright neighbor pixels with marker color *** */
    uint32_t sum_objects = 0;
    {
        dash::Array<uint32_t> sums( numunits, dash::BLOCKED );
        dash::fill( sums.begin(), sums.end(), 0 );
        start= std::chrono::system_clock::now();

        uint32_t lw= matrix.local.extent(1);
        uint32_t lh= matrix.local.extent(0);

        auto foundobjects = sums.lbegin();
        uint64_t visitedpixels = 0;
        uint64_t brightnesssum = 0;

        for ( uint32_t y= 0; y < lh ; y++ ){
            for ( uint32_t x= 0; x < lw ; x++ ){
                *foundobjects += checkobject( matrix.lbegin(), x, y, lw, lh, limit, marker );
                visitedpixels++;
                RGB* pixel= matrix.lbegin() + y*lw+ x;
                brightnesssum += pixel->brightness();
            }
        }

        cout << "unit " << myid << " has " << lw << " x " << lh << " == " << lw*lh
             << " and visited " << visitedpixels << " pixels, found " << *foundobjects
             << ", sum of brightness " << brightnesssum << endl;

        sums.barrier();
        end= std::chrono::system_clock::now();
        if ( 0 == myid ) {
            cout << "marked pixels in parallel in " << std::chrono::duration_cast<std::chrono::seconds> (end-start).count() << " seconds" << endl;
            time_checkobjects= 0.001 * std::chrono::duration_cast<std::chrono::milliseconds> (end-start).count();
        }

        /* combine the local number of objects found into the global result */
        cout << "unit " << myid << " found " << *foundobjects << " objects locally" << endl;

        if ( 0 == myid ) {
            sum_objects = std::accumulate(sums.begin(), sums.end(), 0);
            cout << "found " << sum_objects << " objects in total" << endl;
        }
    }

    matrix.barrier();

    if ( 0 == myid ) {


        /* benchmark output in one line for convenient */
        gethostname( hostname , 100 );

        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        cout << "astro-benchmark timings with " << numunits << " units " <<
            " with distribution-spec " << distspec <<
            " on " << hostname <<
            " at " << std::put_time( std::localtime(&in_time_t), "%Y-%m-%d %X" ) <<
            " creatematrix " << time_creatematrix << " s," <<
            " readtiff " << time_readtiff << " s," <<
            " histogram " << time_histogram << " s," <<
            " searchcolor " << time_searchmarkercolor << " s," <<
            " checkobjects " << time_checkobjects << " s" <<
            " imagefile " << argv[1] <<
            " foundobjects " << sum_objects << flush << endl;
    }

    dash::finalize();
    return 0;
}
