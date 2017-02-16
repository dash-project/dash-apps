#include "JUtils.h"
#include "AStar.h"


int main(int argc, char *argv[])
{
  AStar a;
  JDurationManager dm;
	
  dm.start();
	a.run();
  dm.stop();
  
  a.print();
  dm.print();
  
  return 0;
}

