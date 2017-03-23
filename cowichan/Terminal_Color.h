#ifndef TERMINAL_COLOR_HEADER
#define TERMINAL_COLOR_HEADER

#define TERMINAL_COLORS

#define BEGIN_COLOR "\x1B["
#define END_COLOR   "\x1B[0m"

// Color codes for colorized output on terminal
enum Code {
  FBLK   = 30,
  FRED       ,
  FGREEN     ,
  FYEL       ,
  FBLUE      ,
  FMAG       ,
  FCYN       ,
  FWHT   = 37,
  BRED   = 41,
  BGREEN     ,
  BYEL       ,
  BBLUE      ,
  BMAG       ,
  BCYN       ,
  BWHT   = 47
};


// Helper struct for fmt function
template<typename T>
struct sPar{
  public:
    const T& in;
    const Code code;
    const int width;
    sPar(const T& in_, const Code code_, const int width_ = 0) 
    :in(in_), code(code_), width(width_){};
};

// Function for formated output and more beautiful code
template<typename T>
inline const sPar<T> fmt( const T& in, const Code code, const int width = 0 ) {
  return sPar<T>( in, code, width );
}

std::ostream& operator<<(std::ostream& os, const sPar<unsigned char>& par) {
  #ifdef TERMINAL_COLORS
    return os << BEGIN_COLOR << par.code << "m" << std::setw(par.width) << static_cast<uint>(par.in) << END_COLOR;
  #else
    return os << std::setw(par.width) << static_cast<uint>(par.in);
  #endif
}

// overloading for transfer of several parameters
template<typename T>
std::ostream& operator<<(std::ostream& os, const sPar<T>& par) {
  #ifdef TERMINAL_COLORS
    return os << BEGIN_COLOR << par.code << "m" << std::setw(par.width) << par.in << END_COLOR;
  #else
    return os << std::setw(par.width) << par.in;
  #endif
}



#endif