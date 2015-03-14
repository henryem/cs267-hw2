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
 *   unique_ptr<SimpleIterator<Foo> > i = ...;
 *   while (i->hasNext()) {
 *     Foo f = i->next();
 *     ...
 *   }
 *   delete i; // (if necessary)
 */
template <typename T>
class SimpleIterator {
public:
  virtual ~SimpleIterator() = 0;
  /* True if there is a value to get.  Logically const. */
  virtual bool hasNext() = 0;
  /* Returns the current value and advances the iterator.  It is an error to
   * call next() if hasNext() is false. */
  virtual T next() = 0;
  /* The number of elements left in this iterator.  Calling this exhausts
   * the iterator, so usually it should be called on a deep copy. */
  int size() {
    int num_elements = 0;
    while (hasNext()) {
      next();
      num_elements++;
    }
    return num_elements;
  }
};

template <typename T>
inline SimpleIterator<T>::~SimpleIterator() { }

#endif /* SIMPLEITERATOR_H_ */
