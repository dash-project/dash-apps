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


uint32_t countobjects( dash::NArray<RGB, 2>& matrix, dash::NArray<uint32_t,2>& indexes, uint32_t startindex, uint32_t limit );


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

    auto distspec= dash::DistributionSpec<2>( dash::BLOCKED, dash::NONE );
    //teamspec.balance_extents();
    dash::NArray<RGB, 2> matrix( dash::SizeSpec<2>( imagesize.height, imagesize.width),
        distspec, dash::Team::All(), teamspec );

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
    }

    matrix.barrier();

    if ( 0 == myid ) {
        show_matrix( matrix, 1200, 1000 );
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
    }

    /* from the brightness histogram we learned, that we should define all but the first two histogram bins
    as bright pixels */
    const uint32_t limit= 256*3*2/17;

    matrix.barrier();


    /* *** part 6: create a matching index array and do object coundting there *** */

    {
        dash::NArray<uint32_t,2> indexes( 
            dash::SizeSpec<2>( imagesize.height, imagesize.width ),
            distspec, dash::Team::All(), teamspec );
        dash::fill( indexes.begin(), indexes.end(), 0 );

        dash::Array<uint32_t> sums( numunits, dash::BLOCKED );

        start= std::chrono::system_clock::now();

        sums.local[0]= countobjects( matrix, indexes, 1 + (1<<30) / numunits, limit );
        sums.barrier();
        end= std::chrono::system_clock::now();

        if ( 0 == myid ) {
            uint32_t sum_objects= std::accumulate(sums.begin(), sums.end(), 0);
            cout << "found " << sum_objects << " objects in total in " << 
            std::chrono::duration_cast<std::chrono::seconds> (end-start).count()<< " seconds" << endl;
        }
    }

    dash::finalize();
    return 0;
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

    return found - touchmap.size();
}
