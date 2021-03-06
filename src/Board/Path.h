/* -*- mode: c++ -*- */
/**
 * @file   Path.h
 * @author Sebastien Fourey <http://www.greyc.ensicaen.fr/~seb>
 * @date   Aug 2009
 * 
 * @brief
 */
/*
 * \@copyright This File is part of the Board library which is
 * licensed under the terms of the GNU Lesser General Public Licence.
 * See the LICENCE file for further details.
 */
#ifndef _BOARD_PATH_H_
#define _BOARD_PATH_H_

#include "Board/Point.h"
#include "Board/Rect.h"
#include "Board/Transforms.h"
#include <vector>
#include <iostream>

#ifdef WITH_CAIRO
// cairo
#include <cairo.h>
#endif

namespace LibBoard {



/**
 * The path structure.
 * @brief A path, according to Postscript and SVG definition.
 */
struct Path { 
  
  Path() : _closed( false ) { }

  Path( const std::vector<Point> & points, bool closedPath )
    : _points( points ), _closed( closedPath ) { }

  Path( bool closedPath ) : _closed( closedPath ) { }    
  
  inline void clear();

  inline bool closed() const;

  inline bool empty() const;

  inline unsigned int size() const;

  inline void setClosed( bool closed  );
  
  /** 
   * Barycenter of the path
   * @return 
   */
  Point center() const;

  /** 
   * Add a point at the end of the path.
   * 
   * @param p 
   * 
   * @return 
   */
  Path & operator<<( const Point & p );

  /** 
   * 
   * 
   * 
   * @return 
   */
  Path & pop_back();

  /** 
   * Returns the n-th point of the polyline.
   * 
   * @param i 
   * 
   * @return 
   */
  Point & operator[]( const unsigned int n ) {
    return _points[ n ];
  }

  /** 
   * Returns the n-th point of the polyline.
   * 
   * @param i 
   * 
   * @return 
   */
  const Point & operator[]( const unsigned int n ) const {
    return _points[ n ];
  }

  /** 
   * 
   * 
   * @param angle 
   * @param center 
   * 
   * @return 
   */
  Path & rotate( double angle, const Point & center );

  /** 
   * 
   * 
   * @param angle 
   * @param center 
   * 
   * @return 
   */
  Path rotated( double angle, const Point & center ) const;

  /** 
   * 
   * 
   * @param angle 
   * 
   * @return 
   */
  Path & rotate( double angle );
  
  /** 
   * 
   * 
   * @param angle 
   * 
   * @return 
   */
  Path rotated( double angle ) const;

  /** 
   * 
   * 
   * @param dx 
   * @param dy 
   * 
   * @return 
   */
  Path & translate( double dx, double dy );
 
  /** 
   * 
   * 
   * @param dx 
   * @param dy 
   * 
   * @return 
   */
  Path translated( double dx, double dy ) const;

  /** 
   * 
   * 
   * @param sx 
   * @param sy 
   * 
   * @return 
   */
  Path & scale( double sx, double sy );

  /** 
   * 
   * 
   * @param s 
   * 
   * @return 
   */
  Path & scale( double s );
  
  /** 
   * 
   * 
   * @param sx 
   * @param sy 
   * 
   * @return 
   */
  Path scaled( double sx, double sy )  const;

  Path scaled( double s )  const;

  /** 
   * Scales all the points.
   * 
   * @param s The scaling factor.
   */
  void scaleAll( double s );

  void flushPostscript( std::ostream & stream,
      const TransformEPS & transform ) const;
  
  void flushFIG( std::ostream & stream,
     const TransformFIG & transform ) const;
  
  void flushSVGPoints( std::ostream & stream,
           const TransformSVG & transform ) const;

  void flushSVGCommands( std::ostream & stream,
       const TransformSVG & transform ) const;

#ifdef WITH_CAIRO
  void flushCairoPoints( cairo_t *cr,
     const TransformCairo & transform ) const;
#endif

  void flushTikZPoints( std::ostream & stream,
     const TransformTikZ & transform ) const;

  Rect boundingBox() const;

protected:
  std::vector<Point> _points;
  bool _closed;
};
  
void
Path::clear()
{
  _points.clear();
}

bool
Path::closed() const
{
  return _closed;
}

bool
Path::empty() const
{
  return _points.empty();
}

unsigned int
Path::size() const
{
  return (unsigned int)_points.size();
}

void
Path::setClosed( bool closedPath )
{
  _closed = closedPath;
}

} // namespace LibBoard  

#endif /* _PATH_H_ */

