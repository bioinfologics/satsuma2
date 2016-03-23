#include "Eps.h"

void PrintEpsHeader( ostream& out, 
                     const float horizSize, 
                     const float vertSize,
                     const float border )
{
  const float leftBoundbox = 0;
  const float rightBoundbox = border + horizSize + border;
  const float bottomBoundbox = 0;
  const float topBoundbox = border + vertSize + border;

  out.setf( ios::fixed );
  
  out << "%!PS-Adobe-3.0 EPSF-3.0\n"
      << "%%BoundingBox: "
      << leftBoundbox << " " << bottomBoundbox << " "
      << rightBoundbox << " " << topBoundbox
      << "\n";

  // Translate boundbox to origin.
  out << border << " " <<  border << " translate\n";
}
