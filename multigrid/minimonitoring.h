#ifndef MINIMONITORING_H
#define MINIMONITORING_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <chrono>
#include <ostream>
#include <cstring>

#include <libdash.h>

using namespace std;
using time_point_t = std::chrono::time_point<std::chrono::high_resolution_clock>;

const size_t NAMELEN= 12;

struct MiniMonT {

   uint16_t step; /* 0 == start, 1 == stop, other values as you like */
   std::chrono::time_point<std::chrono::high_resolution_clock> time;
   char name[NAMELEN];
   uint32_t param; /* one extra parameter to use as desired */

   static std::vector<MiniMonT>* tape;
   static std::chrono::time_point<std::chrono::high_resolution_clock> inittime;

   MiniMonT( uint8_t s, const char n[10], uint32_t p= 0 ) {

      step= s;
      time= std::chrono::high_resolution_clock::now();
      strncpy( name, n, NAMELEN ); name[NAMELEN-1]= '\0';
      param= p;
   };

   static void MiniMonInit() {

      tape= new vector<MiniMonT>();
      tape->reserve( 1000 );
      inittime= std::chrono::high_resolution_clock::now();
   }

   static void MiniMonRecord( uint8_t s, const char n[10], uint32_t p= 0 ) {

      tape->push_back( MiniMonT( s, n, p ) );
   }

   static void MiniMonPrint( uint32_t id, uint32_t num );

};

struct cmpByNameOnly {
   bool operator()(const MiniMonT& a, const MiniMonT& b) const {
        return ( 0 > strncmp( a.name, b.name, NAMELEN ) );
   }
};


void MiniMonT::MiniMonPrint( uint32_t id, uint32_t num ) {

   /* print out log to individual files */
   char filename[31];
   std::snprintf( filename, 30, "trace%05u.csv", id );

   ofstream file;
   file.open( filename );

   file << "# Unit;Function_name;Step;Param;Timestamp;"
      "Duration;" << endl;

   /* keep a set of past MiniMonT entries, compare only by name.
   If step increases for the same name then also print the duration.
   This does not work for recursive functions but is otherwise good enough. */

   std::set<MiniMonT,cmpByNameOnly> oldtimes;

   for ( auto it= tape->begin(); it != tape->end(); ++it ) {

      auto old= oldtimes.find( *it );
      if ( old != oldtimes.end() ) {

         if ( old->step +1 == it->step ) {

            /* step values match, print with duration */
            file << id << ";" << it->name << ";" << it->step << ";" << it->param << ";" <<
               std::setprecision(12) <<
               1e-9 * chrono::duration_cast<chrono::nanoseconds>( it->time - inittime ).count() << ";" <<
               std::setprecision(12) <<
               1e-9 * chrono::duration_cast<chrono::nanoseconds>( it->time - old->time ).count()<< ";" << endl;

         } else {

            /* step values don't match, print without duration */
            file << id << ";" << it->name << ";" << it->step << ";" << it->param << ";" <<
               std::setprecision(12) <<
               1e-9 * chrono::duration_cast<chrono::nanoseconds>( it->time - inittime ).count() <<
               ";;" << endl;
         }

         /* replace old entry with current one */
         oldtimes.erase( old );
         oldtimes.insert( *it );

      } else {

         /* no old time value found, print none, print wiothout duration */
         file << id << ";" << it->name << ";" << it->step << ";" << it->param << ";" <<
            std::setprecision(12) <<
            1e-9 * chrono::duration_cast<chrono::nanoseconds>( it->time - inittime ).count() <<
            ";;" << endl;

         oldtimes.insert( *it );
      }
   }

   file.close();
}


std::ostream& operator<< ( std::ostream& stream, const MiniMonT& record );

std::ostream& operator<<( std::ostream& stream, const MiniMonT& record ) {

   stream << record.step << "; " << record.time.time_since_epoch().count() << "; " << record.name;

   return stream;
}




#endif /* MINIMONITORING_H */
