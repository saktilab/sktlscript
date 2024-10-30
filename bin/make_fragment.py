#!/usr/bin/env python3
import sys
import os


# def split_file(filepath, lines_per_file=100):
#     """splits file at `filepath` into sub-files of length `lines_per_file`
#     """
#     lpf = lines_per_file
#     path, filename = os.path.split(filepath)
#     with open(filepath, 'r') as r:
#         name, ext = os.path.splitext(filename)
#         try:
#             w = open(os.path.join(path, '{}_{}{}'.format(name, 0, ext)), 'w')
            
#             for i, line in enumerate(r):
#                 if not i % lpf:
#                     #possible enhancement: don't check modulo lpf on each pass
#                     #keep a counter variable, and reset on each checkpoint lpf.
#                     w.close()
#                     filename = os.path.join(path,
#                                             '{}_{}{}'.format(name, i, ext))
#                     w = open(filename, 'w')
#                 w.write(line)
#         finally:
#             w.close()
def split_file(filepath, lines=100):
    """Split a file based on a number of lines."""
    path, filename = os.path.split(filepath)
    # filename.split('.') would not work for filenames with more than one .
    basename, ext = os.path.splitext(filename)
    # open input file
    with open(filepath, 'r') as f_in:
        try:
            # open the first output file
            f_out = open(os.path.join(path, '{}_{}{}'.format(basename, 0, ext)), 'w')
            # loop over all lines in the input file, and number them
            count = 0
            for i, line in enumerate(f_in):
                # every time the current line number can be divided by the
                # wanted number of lines, close the output file and open a
                # new one
                if i % lines == 0:
                    f_out.close()
                    f_out = open(os.path.join(path, '{}_{}{}'.format(basename, count, ext)), 'w')
                # write the line to the output file
                    count += 1    
                f_out.write(line)
                
        finally:
            # close the last output file
            f_out.close()


bigfile = sys.argv[1]
lines = sys.argv[2]
split_file(sys.argv[1],int(sys.argv[2]))

