#!/usr/bin/env python
'''Usage:
send_softexit.py [options] walltime &

This is a version of watchdog which will work on computers without a
queueing system or using a queueing system other than Torque.

send_softexit must run in the same working directory as the NECI calculation.

Care must be taken that the send_softexit is terminated when the job finishes
(rather than waiting for the send_softexit to finish!).  send_softexit
does, however, listen out for the interrupt signal.  Recommended use in a
script is:
    
send_softexit.py [options] walltime &
[Job commands]
killall -2 send_softexit.py

or

send_softexit.py [options] walltime &
send_softexit_ps=$!
[Job commands]
kill -2 $send_softexit_ps

The latter should be used if multiple calculations are run on one computer.
'''

__author__ = 'James Spencer'

from optparse import OptionParser, OptionValueError
import signal
import sys
import time

def signal_handler(signal, frame):
    '''Capture signal and leave quietly.'''
    print('Signal has been caught.  Bye!')
    sys.exit()

def job_cleanup():
    '''Perform specified operations to perform when the job is approaching the end of its walltime.'''
    f = open('FCIQMC.COMM','w')
    f.write('softexit')
    f.close()


def parse_timer(t):
    '''Parse a time string and return the corresponding (integer) number of seconds.
    
t is a string eiher containing the number of seconds or in the format hh:mm:ss.'''
    if ':' in t:
        # assume hhmmss format
        t = [int(s) for s in t.split(':')]
        return (t[0]*60+t[1])*60+t[2]
    else:
        # already in seconds: just need to convert to integer.
        return int(t)

def parse_options(my_args):
    '''Parse command line options.  Return the amount of sleep time.'''
    parser = OptionParser(usage='''send_softexit.py [options] walltime

Monitor a running job and write SOFTEXIT to CHANGEVARS in the current directory
when the elapsed time gets to within a specified amount of the walltime allowed
for the job.

The walltime and grace period can be given either in seconds or in the format
hh:mm:ss.''')
    parser.add_option('-g','--grace',default='0',help='Amount of time before the walltime expires that SOFTEXIT is sent.  Default: %default.')
    (options,args) = parser.parse_args(my_args)
    if len(args) != 1:
        if len(args) == 0:
            print('Must specify walltime.')
        else:
            print('Do not understand options specified: %s.' % (' '.join(args)))
        parser.print_help()
        sys.exit(1)
    else:
        walltime = parse_timer(args[0])
        grace = parse_timer(options.grace)
        sleep_time = walltime - grace
    return sleep_time

def main(sleep_time):
    print('send_softexit sleeping for %is.' % (sleep_time))
    time.sleep(sleep_time)
    job_cleanup()
    sys.exit()

if __name__ == '__main__':
    sleep_time = parse_options(sys.argv[1:])
    signal.signal(signal.SIGINT,signal_handler) # Listen out for Ctrl-C.
    main(sleep_time)
