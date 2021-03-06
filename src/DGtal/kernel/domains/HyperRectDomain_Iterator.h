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

#pragma once

/**
 * @file HyperRectDomain_Iterator.h
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * @author Guillaume Damiand
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2010/05/31
 *
 * Header file for module HyperRectDomain_Iterator.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(HyperRectDomain_Iterator_RECURSES)
#error Recursive header files inclusion detected in HyperRectDomain_Iterator.h
#else // defined(HyperRectDomain_Iterator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define HyperRectDomain_Iterator_RECURSES

#if !defined HyperRectDomain_Iterator_h
/** Prevents repeated inclusion of headers. */
#define HyperRectDomain_Iterator_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//#include <iterator> // Bug for operator * => dangling reference !!!
// Class allowing to build a reverse iterator of a given iterator.
template<typename _Iterator>
class myreverse_iterator
  : public iterator<typename iterator_traits<_Iterator>::iterator_category,
  typename iterator_traits<_Iterator>::value_type,
  typename iterator_traits<_Iterator>::difference_type,
  typename iterator_traits<_Iterator>::pointer,
  typename iterator_traits<_Iterator>::reference>
{
protected:
  _Iterator current;
  _Iterator prev;

public:
  typedef _Iterator                 iterator_type;
  typedef typename iterator_traits<_Iterator>::difference_type
      difference_type;
  typedef typename iterator_traits<_Iterator>::reference   reference;
  typedef typename iterator_traits<_Iterator>::pointer     pointer;

public:
  explicit
      myreverse_iterator(iterator_type __x) : current(__x),
      prev(current)
  { --prev; }

  myreverse_iterator(const myreverse_iterator& __x)
    : current(__x.current), prev(__x.prev) { }

  iterator_type base() const
  { return current; }

  /*const*/ reference operator*() const
  { return *prev; }

  reference operator*()
  { return *prev; }

  pointer operator->() const
  { return &(operator*()); }

  myreverse_iterator& operator++()
  { --current; --prev;
    return *this;
  }

  myreverse_iterator operator++(int)
  {
    myreverse_iterator __tmp = *this;
    operator++();
    return __tmp;
  }

  myreverse_iterator& operator--()
  {
    ++current; ++prev;
    return *this;
  }

  myreverse_iterator operator--(int)
  {
    myreverse_iterator __tmp = *this;
    operator--();
    return __tmp;
  }

  myreverse_iterator operator+(difference_type __n) const
  { return myreverse_iterator(current - __n); }

  myreverse_iterator& operator+=(difference_type __n)
                                {
    current -= __n; prev = current; --prev;
    return *this;
  }

  myreverse_iterator operator-(difference_type __n) const
  { return myreverse_iterator(current + __n); }

  myreverse_iterator& operator-=(difference_type __n)
                                {
    current += __n; prev = current; --prev;
    return *this;
  }

  reference operator[](difference_type __n) const
  { return *(*this + __n); }
};
template<typename _Iterator>
inline bool
    operator==(const myreverse_iterator<_Iterator>& __x,
               const myreverse_iterator<_Iterator>& __y)
{ return __x.base() == __y.base(); }
template<typename _Iterator>
inline bool
    operator!=(const myreverse_iterator<_Iterator>& __x,
               const myreverse_iterator<_Iterator>& __y)
{ return !(__x == __y); }

//******************************************************************************
namespace DGtal
{
  /////////////////////////////////////////////////////////////////////////////
  // class HyperRectDomain_Iterator
  /**
   * Description of class 'HyperRectDomain_Iterator' <p>
   * Aim:
   */
  template<typename TPoint>
  class HyperRectDomain_Iterator
  {
  public:
    typedef std::bidirectional_iterator_tag iterator_category; ///\todo construct a RANDOM-ACCESS iterator
    typedef TPoint value_type;
    typedef ptrdiff_t difference_type;
    typedef TPoint* pointer;
    typedef TPoint& reference;
    typedef typename TPoint::Dimension Dimension;


    HyperRectDomain_Iterator( const TPoint & p, const TPoint& lower, const TPoint &upper )
      : myPoint( p ), mylower( lower ), myupper( upper ),  myCurrentPos( 0 )
    {
      ASSERT( lower <= upper );
      ASSERT( lower <= p && p <= upper );
    }

    const TPoint & operator*() const
    {
      ASSERT(mylower<=myPoint && myPoint<=myupper); // we must be between [begin,end]
      return myPoint;
    }
    TPoint & operator*()
    {
      ASSERT(mylower<=myPoint && myPoint<=myupper); // we must be between [begin,end]
      return myPoint;
      }

    /**
     * Operator ==
     *
     */
    bool operator== ( const HyperRectDomain_Iterator<TPoint> &it ) const
    {
      return ( myPoint==it.myPoint );
    }

    /**
     * Operator !=
     *
     */
    bool operator!= ( const HyperRectDomain_Iterator<TPoint> &aIt ) const
    {
      return ( myPoint!=aIt.myPoint );
    }

    /**
     * Implements the next() method to scan the domain points dimension by dimension
     * (lexicographic order).
     **/
    void nextLexicographicOrder()
    {
      ++myPoint[myCurrentPos];
      if (( myCurrentPos < TPoint::dimension - 1 ) &&
    ( myPoint[myCurrentPos] > myupper[myCurrentPos] ) )
        {
          do
      {
        myPoint[myCurrentPos] = mylower[myCurrentPos];
        myCurrentPos++;
        if ( myCurrentPos < TPoint::dimension )
    ++myPoint[myCurrentPos];
      }
          while (( myCurrentPos < TPoint::dimension - 1 ) &&
     ( myPoint[myCurrentPos]  >  myupper[ myCurrentPos ] ) );
          myCurrentPos = 0;
        }
    }

    /**
     * Operator ++ (++it)
     */
    HyperRectDomain_Iterator<TPoint> &operator++()
    {
      nextLexicographicOrder();
      return *this;
    }

    /**
     * Operator ++ (it++)
     *
     */
    HyperRectDomain_Iterator<TPoint> operator++ ( int )
    {
      HyperRectDomain_Iterator<TPoint> tmp = *this;
      nextLexicographicOrder();
      return tmp;
    }

    /**
     * Implements the prev() method to scan the domain points dimension by dimension
     * (lexicographic order).
     **/
    void prevLexicographicOrder()
    {
      --myPoint[ myCurrentPos ];
      if (( myCurrentPos < TPoint::dimension - 1 ) &&
    ( myPoint[ myCurrentPos ]  <  mylower[ myCurrentPos ] ) )
        {
          do
      {
        myPoint[ myCurrentPos ] = myupper[ myCurrentPos ];
        ++myCurrentPos;
        if ( myCurrentPos < TPoint::dimension )
    --myPoint[ myCurrentPos ];
      }
          while (( myCurrentPos < TPoint::dimension - 1 ) &&
     ( myPoint[ myCurrentPos ]  <  mylower[ myCurrentPos ] ) );
          myCurrentPos = 0;
        }
    }

    /**
     * Operator -- (--it)
     *
     */
    HyperRectDomain_Iterator<TPoint> &operator--()
    {
      prevLexicographicOrder();
      return *this;
    }

    /**
     * Operator -- (it--)
     */
    HyperRectDomain_Iterator<TPoint> operator-- ( int )
    {
      HyperRectDomain_Iterator<TPoint> tmp = *this;
      prevLexicographicOrder();
      return tmp;
    }

  private:
    ///Current Point in the domain
    TPoint myPoint;
    ///Copies of the Domain limits
    TPoint mylower, myupper;
    ///Second index of the iterator position
    Dimension myCurrentPos;
  };
  /////////////////////////////////////////////////////////////////////////////
  // class HyperRectDomain_Iterator
  /**
   * Description of class 'HyperRectDomain_Iterator' <p>
   * Aim:
   */
  template<typename TPoint>
  class HyperRectDomain_subIterator
  {
  public:
    typedef std::bidirectional_iterator_tag iterator_category; ///\todo construct a RANDOM-ACCESS iterator
    typedef TPoint value_type;
    typedef ptrdiff_t difference_type;
    typedef TPoint* pointer;
    typedef TPoint& reference;
    typedef typename TPoint::Dimension Dimension;

#ifdef CPP11_INITIALIZER_LIST
    HyperRectDomain_subIterator(const TPoint & p, const TPoint& lower,
        const TPoint &upper,
        std::initializer_list<Dimension> subDomain)
      : myPoint( p ), mylower( lower ), myupper( upper ),  myCurrentPos( 0 )
    {
      ASSERT( lower <= upper );
      ASSERT( lower <= p && p <= upper );
      ASSERT( subDomain.size() <= TPoint::dimension );
      mySubDomain.reserve( subDomain.size() );
      for ( const unsigned int *c = subDomain.begin();
            c != subDomain.end(); ++c )
        {
          ASSERT( *c <= TPoint::dimension );
          mySubDomain.push_back( *c );
        }

      // TODO: check the validity of the subDomain ?
    }
#endif
    HyperRectDomain_subIterator(const TPoint & p, const TPoint& lower,
           const TPoint &upper,
           const std::vector<Dimension> &subDomain)
      : myPoint( p ), mylower( lower ), myupper( upper ),  myCurrentPos( 0 )
    {
      ASSERT( lower <= upper );
      ASSERT( lower <= p && p <= upper );
      ASSERT( subDomain.size() <= TPoint::dimension );
      mySubDomain.reserve( subDomain.size() );
      for ( typename std::vector<Dimension>::const_iterator it = subDomain.begin();
            it != subDomain.end(); ++it )
        {
          ASSERT( *it <= TPoint::dimension );
          mySubDomain.push_back( *it );
        }

      // TODO: check the validity of the subDomain ?
    }


    const TPoint & operator*() const
    {
      ASSERT(mylower<=myPoint && myPoint<=myupper); // we must be between [begin,end]
      return myPoint;
    }
    TPoint & operator*()
    {
      ASSERT(mylower<=myPoint && myPoint<=myupper); // we must be between [begin,end]
      return myPoint;
    }

    /**
     * Operator ==
     *
     */
    bool operator== ( const HyperRectDomain_subIterator<TPoint> &it ) const
    {
			for (unsigned int i=0; i<mySubDomain.size(); ++i)
				if ( myPoint[mySubDomain[i]]!=it.myPoint[mySubDomain[i]]) return false;
			return true;
			//  return ( myPoint==it.myPoint );
    }

    /**
     * Operator !=
     *
     */
    bool operator!= ( const HyperRectDomain_subIterator<TPoint> &aIt ) const
    {
      return !operator==(aIt);
			// ( myPoint!=aIt.myPoint );
    }

    /**
     * Implements the next() method to scan the domain points dimension by dimension
     * (by using the subDomain order given by the user).
     **/
    void nextSubDomainOrder()
    {
      ASSERT( myCurrentPos < mySubDomain.size() );
      ++myPoint[ mySubDomain[myCurrentPos] ];

      if ( myCurrentPos < mySubDomain.size() - 1 &&
     myPoint[ mySubDomain[myCurrentPos] ] >
     myupper[ mySubDomain[myCurrentPos] ] )
        {
          do
      {
        myPoint[ mySubDomain[myCurrentPos] ] =
    mylower[ mySubDomain[myCurrentPos] ];
        ++myCurrentPos;
        if ( myCurrentPos < mySubDomain.size() )
    ++myPoint[ mySubDomain[myCurrentPos] ];
      }
          while (( myCurrentPos < mySubDomain.size() - 1  ) &&
     ( myPoint[ mySubDomain[myCurrentPos] ]  >
       myupper[ mySubDomain[myCurrentPos] ] ) );
          myCurrentPos = 0;
        }
    }

    /**
     * Operator ++ (++it)
     */
    HyperRectDomain_subIterator<TPoint> &operator++()
    {
      nextSubDomainOrder();
      return *this;
    }

    /**
     * Operator ++ (it++)
     *
     */
    HyperRectDomain_subIterator<TPoint> operator++ ( int )
    {
      HyperRectDomain_subIterator<TPoint> tmp = *this;
      nextSubDomainOrder();
      return tmp;
    }

    /**
     * Implements the prev() method to scan the domain points dimension by dimension
     * (subDomain order).
     **/
    void prevSubDomainOrder()
    {
      ASSERT( myCurrentPos < mySubDomain.size() );
      --myPoint[ mySubDomain[myCurrentPos] ];

      if (  myCurrentPos < mySubDomain.size() - 1 &&
            myPoint[ mySubDomain[myCurrentPos] ]  <
            mylower[ mySubDomain[myCurrentPos] ] )
        {
          do
      {
        myPoint[ mySubDomain[myCurrentPos] ] =
    myupper[ mySubDomain[myCurrentPos] ];
        ++myCurrentPos;
        if ( myCurrentPos < mySubDomain.size() )
    --myPoint[ mySubDomain[myCurrentPos] ];
      }
          while (( myCurrentPos < mySubDomain.size() - 1 ) &&
     ( myPoint[ mySubDomain[myCurrentPos] ]  <
       mylower[ mySubDomain[myCurrentPos] ] ) );
          myCurrentPos = 0;
        }
    }

    /**
     * Operator -- (--it)
     *
     */
    HyperRectDomain_subIterator<TPoint> &operator--()
    {
      prevSubDomainOrder();
      return *this;
    }

    /**
     * Operator -- (it--)
     */
    HyperRectDomain_subIterator<TPoint> operator-- ( int )
    {
      HyperRectDomain_subIterator<TPoint> tmp = *this;
      prevSubDomainOrder();
      return tmp;
    }

  private:
    ///Current Point in the domain
    TPoint myPoint;
    ///Copies of the Domain limits
    TPoint mylower, myupper;
    ///Second index of the iterator position
    Dimension myCurrentPos;
    ///Vector of subDomain on dimension, to fix the order in which dimensions
    /// are considered.
    std::vector<Dimension> mySubDomain;
  };

} //namespace
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined HyperRectDomain_Iterator_h

#undef HyperRectDomain_Iterator_RECURSES
#endif // else defined(HyperRectDomain_Iterator_RECURSES)
