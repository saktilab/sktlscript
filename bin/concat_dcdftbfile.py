#!/usr/bin/env python3

import sys

print(sys.argv)
folders = sys.argv[1:]

def ordered_set(in_list):
    out_list = []
    added = set()
    for val in in_list:
        if not val in added:
            out_list.append(val)
            added.add(val)
    return out_list

# job_ids = ordered_set(list(map(int, sys.argv[1:])))


# mulliken
print('Merging mulliken...')
stepnos = set()
with open('merged_nac', 'w') as fout:
    for i, jobid in enumerate(folders):
 #       try:
            with open('{}/mulliken'.format(jobid), 'r') as fin:
                for line in fin:
                    arr = line.split()
                    norb = int(arr[0])
                    nat = int(arr[1])
                    title = next(fin)
                    arr = title.split()
                    stepno = int(arr[9])
                    if (stepno in stepnos):
                        for j in range(norb):
                            next(fin)
                        continue
                    else:
                        stepnos.add(stepno)
                        print(nat, file=fout)
                        print(title, end='', file=fout)
                        nacs = {}
                        for j in range(norb):
                            line = next(fin)
                            arr = line.split()
                            iat = int(arr[0])
                            if (iat in nacs):
                                nacs[iat] += float(arr[3])
                            else:
                                nacs[iat] = float(arr[3])

                        for j in range(nat):
                            print('{} {:.6f}'.format(j+1, nacs[j+1]), file=fout)
#        except: 
#            print('Error:' , 'dftb.mul.{}'.format(jobid))
#            break

file_list = [ 'traject', 'velocity']
name_list = [ 'traject', 'merged_velocity']
for file, name in zip(file_list, name_list):
    print('Merging {}...'.format(file))
    stepnos = set()
    # print(file, name)
    with open(name, 'w') as fout:
        for i, jobid in enumerate(folders):
            try:
                with open('{}/{}'.format(jobid, file), 'r') as fin:
                    for line in fin:
                        nat = int(line)
                        dummy = next(fin)
                        arr = dummy.split()
                        stepno = int(arr[9])
                        if (stepno in stepnos):
                            for j in range(nat):
                                next(fin)
                            continue
                        else:
                            stepnos.add(stepno)
                            print(nat, file=fout)
                            print(dummy, end='', file=fout)
                            for j in range(nat):
                                print(next(fin), end='', file=fout)
            except:
                print('Error:' , name, '{}/{}'.format(jobid, file ))
                break
	      
    

print('Merging dftb.out...')
stepnos = set()
with open('merged_out', 'w') as fout:
  print('#{:>12s} {:>20s} {:>20s} {:>20s} {:>20s}'.format('Time', 'Temperature', 'Poten. Energy', 'Kinetic Energy', 'MD Energy'), file=fout)
  for i, jobid in enumerate(folders):
    try:
      with open('{}/{}'.format(jobid, 'dftb.out'), 'r') as fin: 
        for line in fin:
          if "FSEC, THIS RUN'S STEP NO.=" in line:
            arr = line.split()
            stepno = int(arr[9])
            time = float(arr[3])
            if (stepno in stepnos):
              continue
            else:
              stepnos.add(stepno)
              next(fin)
              line = next(fin)
              arr = line.split()
              temp = float(arr[2])
              line = next(fin)
              arr = line.split()
              ke = float(arr[3]) 
              line = next(fin)
              arr = line.split()
              mdene = float(arr[4]) 
              pe = mdene - ke 
              print(' {:12.4f} {:20.12f} {:20.12f} {:20.12f} {:20.12f}'.format(time, temp, pe, ke, mdene), file=fout)
    except:
      print('Warning:', 'dftb.out ', jobid)
      continue
