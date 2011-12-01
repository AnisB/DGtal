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
 * @file imageGridCurveEstimator.cpp
 * @ingroup tutorial-examples
 * @author Tristan Roussillon (tristan.roussillon@liris.cnrs.fr)
 *
 *
 * @date 2010/10/17
 * 
 * @brief An example of extracting a grid curve from an image iso-contour
 * and estimating its length. 
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <algorithm>
///////////////////////////////////////////////////////////////////////////////

//! [imageGridCurveEstimator-basicIncludes]
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "ConfigExamples.h"
//! [imageGridCurveEstimator-basicIncludes]

//! [imageGridCurveEstimator-imageIncludes]
#include "DGtal/io/readers/PNMReader.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/imagesSetsUtils/IntervalForegroundPredicate.h"
//! [imageGridCurveEstimator-imageIncludes]

//! [imageGridCurveEstimator-trackingIncludes]
#include "DGtal/topology/helpers/Surfaces.h"
//! [imageGridCurveEstimator-trackingIncludes]

//! [imageGridCurveEstimator-estimatorIncludes]
#include "DGtal/geometry/2d/estimators/DSSLengthEstimator.h"
//! [imageGridCurveEstimator-estimatorIncludes]

//display
#include "DGtal/io/boards/Board2D.h"

//segmentation
#include "DGtal/geometry/2d/GreedySegmentation.h"

///////////////////////////////////////////////////////////////////////////////

int main()
{
  //image import
  typedef DGtal::ImageContainerBySTLVector< Z2i::Domain, int> Image;
  std::string filename =  examplesPath + "samples/contourS.pgm";
  Image image = DGtal::PNMReader<Image>::importPGMImage(filename); 

  //! [imageGridCurveEstimator-predicate] 
  //predicate from the image
  IntervalForegroundPredicate<Image> predicate(image,0,135); 
  //! [imageGridCurveEstimator-predicate]

  //! [imageGridCurveEstimator-prepareTracking]
  Z2i::KSpace ks;                                            //Khalimsky space 
  ks.init( image.lowerBound(), image.upperBound(), true );
  SurfelAdjacency<2> sAdj( true );                           //adjacency
  //! [imageGridCurveEstimator-prepareTracking]

  //! [imageGridCurveEstimator-tracking]
  //extraction of all the contours
  std::vector< std::vector< Z2i::SCell > > contours;
  Surfaces<Z2i::KSpace>
    ::extractAll2DSCellContours( contours, ks, sAdj, predicate );
  //! [imageGridCurveEstimator-tracking]

  if (contours.size() > 0)
  {
    
    //! [imageGridCurveEstimator-instantiation]
    //init grid curve from the first retrieved contour
    Z2i::Curve c;
    c.initFromSCellsVector( contours.at(1) );  
    //! [imageGridCurveEstimator-instantiation]

    //! [imageGridCurveEstimator-getRange]
    //range of points
    typedef Z2i::Curve::PointsRange Range; 
    Range r = c.getPointsRange(); 
    //! [imageGridCurveEstimator-getRange]

    //! [imageGridCurveEstimator-lengthEstimation]
    //length estimation based on a DSS segmentation
    DSSLengthEstimator< Range::ConstIterator > DSSlength;
    DSSlength.init(1, r.begin(), r.end(), c.isClosed());
    double length = DSSlength.eval();
    trace.info() << "Length: " << length << endl; 
    //! [imageGridCurveEstimator-lengthEstimation]
    
    //DSS segmentation display
    typedef Z2i::Curve::PointsRange::ConstCirculator ConstCirculator; 
    typedef ArithmeticalDSS<ConstCirculator,int,4> SegmentComputer;
    typedef GreedySegmentation<SegmentComputer> Segmentation;

    Segmentation theSegmentation( r.c(), r.c(), SegmentComputer() );
    Segmentation::SegmentComputerIterator i = theSegmentation.begin();
    Segmentation::SegmentComputerIterator end = theSegmentation.end();
    
    DGtal::Board2D aBoard;
    aBoard << SetMode("PointVector", "Grid");
    for ( ; i != end; ++i) {
      aBoard << SetMode(i->className(), "Points") << *i; 
      aBoard << SetMode(i->className(), "BoundingBox") << *i; 
    } 
    aBoard.saveEPS("DisplayDSSSegmentationTuto3.eps");
  
  } else trace.info() << "no contour" << endl; 
  
  return 0;

}

///////////////////////////////////////////////////////////////////////////////
