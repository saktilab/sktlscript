#!/usr/bin/env python3

import sys
def each_slice(size, iterable):
    """ Chunks the iterable into size elements at a time, each yielded as a list.

    Example:
      for chunk in each_slice(2, [1,2,3,4,5]):
          print(chunk)

      # output:
      [1, 2]
      [3, 4]
      [5]
    """
    current_slice = []
    for item in iterable:
        current_slice.append(item)
        if len(current_slice) >= size:
            yield current_slice
            current_slice = []
    if current_slice:
        yield current_slice


with open(sys.argv[1], 'r') as f:
    nat = int(f.readline().split()[0])

subs  = list(each_slice(10, [1]*nat))

for sub in subs:
    print (' '.join(map(str, sub)))

print ()
