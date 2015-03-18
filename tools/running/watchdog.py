#!/usr/bin/env python
'''Monitor a running job and run a cleanup function when the elapsed time gets to within a specified amount of the walltime allowed for the job.

Usage:
watchdog.py [options] job_id &

This should be called from your PBS script, so that watchdog is automatically running in the correct working directory.

job_id needs to be the full job id, as used internally by torque (e.g. 47345.master.beowulf.cluster on tardis).

Care must be taken that the watchdog is terminated when the job finishes (rather than PBS waiting for the watchdog to finish!).  
watchdog does, however, listen out for the interrupt signal.  Recommended use in a PBS script is:
    
watchdog.py [options] $PBS_JOBID &
[Job commands]
killall -2 watchdog.py '''

__author__ = 'James Spencer'

import signal,sys,time
from optparse import OptionParser,OptionValueError
import PBSQuery

exit_time = 900

def signal_handler(signal,frame):
    '''Capture interupt signal and leave quietly.'''
    print('Interrupt signal has been caught.  Bye!')
    sys.exit()

def job_cleanup():
    '''Perform specified operations to perform when the job is approaching the end of its walltime.'''
    f = open('HANDE.COMM','w')
    f.write('softexit')
    f.close()

def hhmmss_to_seconds(hhmmss):
    '''Convert time in the string format of hh:mm:ss to seconds.'''
    t = [int(s) for s in t.split(':')]
    return (t[0]*60+t[1])*60+t[2]

def parse_options(my_args):
    '''Parse command line options.'''
    parser = OptionParser(usage=__doc__)
    parser.add_option('-e','--exit-time',type='float',default=exit_time,help='Amount of time (in seconds) before the walltime runs out at which the job_cleanup function is called to terminate the job. Default=%defaults.')
    (options,args) = parser.parse_args(my_args)
    if len(args) != 1:
        if len(args) == 0:
            print('Must specify the job id.')
        else:
            print('Do not understand options specified: %s.' % (' '.join(args)))
        parser.print_help()
        sys.exit()
    else:
        job_id = args[0]
    return (options,job_id)

def main(job_id):
    '''Monitor a job_id and run job_cleanup when the elapsed_time is within exit_time seconds of the wall_time.'''
    p = PBSQuery.PBSQuery()
    if not p.getjob(job_id):
        raise Exception('invalid job id %s.' % job_id)
    job = p.getjob(job_id)
    wall_time = hhmmss_to_seconds(job[job_id]['Resource_List.walltime'])
    try:
        sleep_time = wall_time-options.exit_time
    except NameError:
        # Try to use the default in case watchdog is being used as a module.
        sleep_time = wall_time-exit_time
    print('Watchdog sleeping for %i' % (sleep_time))
    time.sleep(sleep_time)
    job_cleanup()
    sys.exit

if __name__ == '__main__':
    (options,job_id) = parse_options(sys.argv[1:])
    signal.signal(signal.SIGINT,signal_handler) # Listen out for Ctrl-C.
    main(job_id)
