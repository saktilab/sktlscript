#!/usr/bin/env python3

import sys
import time

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from termcolor import colored

from subprocess import Popen, PIPE

process = Popen(["qstat", "-xml", *sys.argv[1:]], stdout=PIPE)
(output, err) = process.communicate()
exit_code = process.wait()


tree = ET.ElementTree(ET.fromstring(output))
root = tree.getroot()


def find_path(number):
    process = Popen(["qstat", "-j", number], stdout=PIPE)
    (output, err) = process.communicate()
    exit_code = process.wait()
    lines = str(output).split('\\n')
    for line in lines:
        if 'sge_o_workdir' in line:
            arr = line.split(':')
            return arr[1].strip()
    return ''

def print_job(job):
    if ('JB_submission_time' in job):
        sgetime = job['JB_submission_time']    
    elif ('JAT_start_time'):
        sgetime = job['JAT_start_time']
    time_obj = time.strptime(sgetime, "%Y-%m-%dT%H:%M:%S")   #2017-10-14T21:43:22
    str_time = time.strftime("%y/%m/%d %H:%M:%S", time_obj)
    print(' {} {} {:8s} {:40s} {:4s} {} {:2s} {}'.format(job['JB_job_number'], job['JAT_prio']
    , job['JB_owner'], job['JB_name'], job['state'], str_time, job['slots'], job['queue_name']
    
    ))

def print_job_path(job):
    print(' {} {:40s} {}'.format(job['JB_job_number'], 
    job['JB_name'], job['path']
    
    ))




jobs = []
for elem in tree.iter(tag='job_list'):
    
    res = {}
    res['status'] = elem.attrib['state']

    for item in elem.iter():
        # print (item)
        res[item.tag] = item.text

    jobs.append(res)

print('')
print('Number of jobs: ', len(jobs))
print (' JOB Status')
print("-"*100)    
for item in jobs:
    print_job (item)
print()
print(' JOB Path')
print("-"*100)    
for item in jobs:
    path = find_path(item['JB_job_number'])
    item['path'] = path
    print_job_path (item)
print()

