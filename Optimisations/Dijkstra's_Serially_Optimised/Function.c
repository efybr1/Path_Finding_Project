
// Heap-manipulation functions: push_heap, pop_heap, make_heap, sort_heap,
  // + is_heap and is_heap_until in C++0x.

  template<typename _RandomAccessIterator, typename _Distance, typename _Tp,
	   typename _Compare>
    _GLIBCXX20_CONSTEXPR
    void
    __push_heap(_RandomAccessIterator __first,
		_Distance __holeIndex, _Distance __topIndex, _Tp __value,
		_Compare& __comp)
    {
      _Distance __parent = (__holeIndex - 1) / 2;
      while (__holeIndex > __topIndex && __comp(__first + __parent, __value))
	{
	  *(__first + __holeIndex) = _GLIBCXX_MOVE(*(__first + __parent));
	  __holeIndex = __parent;
	  __parent = (__holeIndex - 1) / 2;
	}
      *(__first + __holeIndex) = _GLIBCXX_MOVE(__value);
    }

    //The above is the push_heap function of stl_heap.h

    //The priority queue is an ordered vector of wrappers in the priority queue order
    //When you initialise an element it needs to know its index in the vector

    //Swap always uses a copy constructor and an = operator
    //Each wrapper should know it's own location
    //We change the copy constructor as this is called by swap.

    //std::swap
