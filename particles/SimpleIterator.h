/*
 * SimpleIterator.h
 *
 *  Created on: Feb 24, 2015
 *      Author: henrym
 */

#ifndef SIMPLEITERATOR_H_
#define SIMPLEITERATOR_H_

/* A Java-style iterator.  IMO, it is easier to implement and use than STL
 * iterators when iterating over non-traditional (perhaps procedurally-defined
 * collections.
 *
 * The standard iteration idiom is:
 *   SimpleIterator<Foo> i = ...;
 *   while (i.hasNext()) {
 *     Foo& f = i.next();
 *     ...
 *   }
 *   delete i; // (if necessary)
 */
template <typename T>
class SimpleIterator {
public:
  virtual ~SimpleIterator() = 0;
  /* True if there is a value to get. */
  virtual bool hasNext() const = 0;
  /* Returns the current value and advances the iterator.  It is an error to
   * call next() if hasNext() is false. */
  virtual T next() = 0;
};

template <typename T>
inline SimpleIterator<T>::~SimpleIterator() { }

#endif /* SIMPLEITERATOR_H_ */
