/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file testEllipsoid.cpp
 * @ingroup Tests
 * @author Anis Benyoub (\c anis.benyoub@insa-lyon.fr )
 *
 * @date 2012/°6/05
 *
 * Functions for testing class Ball3D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <QtGui/QApplication>
#include "DGtal/shapes/parametric/Flower3D.h"
#include "DGtal/helpers/StdDefs.h"
#include <iostream>
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/io/Color.h"
#include "DGtal/kernel/sets/SetPredicate.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/helpers/BoundaryPredicate.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/topology/SCellsFunctors.h"
///////////////////////////////////////////////////////////////////////////////

 using namespace std;
 using namespace DGtal;
 using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :


 int main(int argc, char** argv)
 {

     
   
   // -------------------------------------------------------------------------- Type declaring
   typedef Space::RealPoint RealPoint;
   typedef Flower3D<Space> EuclideanShape;
   typedef GaussDigitizer<Space,EuclideanShape> DigitalShape;
   typedef ImageContainerBySTLVector<Domain,DGtal::uint8_t> Image;

    
   
   // -------------------------------------------------------------------------- Creating the shape
    RealPoint c1(0, 0, 0 );
    EuclideanShape ball1( c1, 20, 4, 2, 0.0); 
	       
   // -------------------------------------------------------------------------- GaussDigitizing
    DigitalShape dshape;
    dshape.attach( ball1 );
    RealPoint p1 =RealPoint( -30.0, -30.0, -30.0 );
    RealPoint p2 =RealPoint( 30.0, 30.0, 30.0 );
    dshape.init( RealPoint( p1 ), RealPoint( p2 ), 1.0);
    Domain domain = dshape.getDomain();

   
   // -------------------------------------------------------------------------- Khalimskhy
    KSpace K;
    bool space_ok = K.init( domain.lowerBound(), domain.upperBound(), true );
    if (!space_ok)
    {
      return 2;
    }

   
    // -------------------------------------------------------------------------- Other types
    typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
    typedef KSpace::Surfel Surfel;
    typedef KSpace::SurfelSet SurfelSet;
    typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
    typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;


    // -------------------------------------------------------------------------- Tracking the boudnadry
    MySurfelAdjacency surfAdj( true ); // interior in all directions.
    MySetOfSurfels theSetOfSurfels( K, surfAdj );
    Surfel bel = Surfaces<KSpace>::findABel( K, dshape, 1000 );
    Surfaces<KSpace>::trackBoundary( theSetOfSurfels.surfelSet(),  K, surfAdj, dshape, bel );

    

    QApplication application(argc,argv);
    Viewer3D viewer;
    viewer.show();
    viewer << SetMode3D( domain.className(), "BoundingBox" ) << domain;





//-----------------------------------------------------------------------
//Specifing a color map

  GradientColorMap<double> cmap_grad( 0, 10 );
  cmap_grad.addColor( Color( 50, 50, 255 ) );
  cmap_grad.addColor( Color( 255, 0, 0 ) );
  cmap_grad.addColor( Color( 255, 255, 10 ) );


//-----------------------------------------------------------------------
//Drawing the ball  && giving a color to the surfels( depending on the
//curvature)

  unsigned int nbSurfels = 0;
  SCellToMidPoint<KSpace> midpoint(K);

  for ( std::set<SCell>::iterator it = theSetOfSurfels.begin(), it_end = theSetOfSurfels.end(); 
	it != it_end; ++it, ++nbSurfels )
  {
    RealPoint A = midpoint( *it );

    DGtal::StarShaped3D<Space>::AngularCoordinates  Angles= ball1.parameter(A);
    double curvature =ball1.meanCurvature(Angles);

//    cout<<"Gaussian"<<curvature<<endl;

//    cout<<"Mean"<<ball1.meanCurvature(Angles)<<endl;
    viewer <<   CustomColors3D( Color::Black, cmap_grad( curvature));
    viewer << *it;
  }



  viewer << Viewer3D::updateDisplay;

  return application.exec();
}
//                                                                           //
/////////////////////////////////////////////////////////////////////////////// 


 


