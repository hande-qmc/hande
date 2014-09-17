#!/usr/bin/env python

import argparse
import os
import sys
import socket
try:
    import socketserver
except:
    import SocketServer as socketserver

class DataServer(socketserver.BaseRequestHandler):
    def handle(self):
        filename = self.request.recv(1024)
        output = main(filename.decode('utf-8'))
        self.request.sendall(output.encode('utf-8'))

def main(filename):

    import pyhande

    (metadata, qmc_data, other_calcs) = pyhande.extract.extract_data(filename)

    output = []
    if not qmc_data.empty:
        qmc_data.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
        # Compare every 1/4 of the calculation...
        indxs = [0] + [int((i*len(qmc_data))/4)-1 for i in range(1,5)]
        test_data = qmc_data.ix[indxs]
        output.append(
                (metadata[(0, 'calc_type')], test_data.to_string(index=False, index_names=False))
                )
    for calc in other_calcs:
        if 'FCI' in calc.name:
            # Compare at most the first 5 eigenvalues (just to be quick/minimise
            # test output/handle different orders of states within Lanczos):
            select_data = calc[:min(len(calc),5)].to_frame().T
            select_data.rename(columns=lambda x: 'eigv_%s' % (x,), inplace=True)
        elif 'Hilbert space' in calc.name:
            select_data = calc.to_frame().T
            select_data.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
        output.append((calc.name, select_data.to_string(index=False)))

    return '\n\n'.join('\n'.join(calc_out) for calc_out in output)

def parse_args(args):

    parser = argparse.ArgumentParser(description='Extract selected data from HANDE output files for validation tests.',
            epilog='If --socket is specifed and a file is not provided, a data server is launched using the desired port.  '
                   'The server can then be used by other instances to avoid the overhead of module imports.')
    parser.add_argument('-s', '--socket', metavar='PORT', nargs='?', type=int, const=50007, default=-1,
            help='Use/launch a data server.  Not used by default.  Default port: %(const)s.')
    parser.add_argument('file', nargs='?', default='', help='File to be processed.')
    opts = parser.parse_args(args)
    if opts.socket <= 0 and not opts.file:
        parser.print_usage()
        print('Must specify a file if not launching a data server.')
        sys.exit()
    return (opts.socket, opts.file)

def launch_server(port, host=''):

    socketserver.TCPServer.allow_reuse_address = True
    server = socketserver.TCPServer((host, port), DataServer)
    server.serve_forever()

def use_server(filename, port, host=''):

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((host, port))
    filename = os.path.join(os.getcwd(), filename)
    s.sendall(filename.encode('utf-8'))
    output = s.recv(1024).decode('utf-8')
    s.close()
    return output

if __name__ == '__main__':

    (port, filename) = parse_args(sys.argv[1:])
    host = ''
    if filename:
        if port > 0:
            try:
                output = use_server(filename, port, host)
            except ConnectionRefusedError:
                output = main(filename)
        else:
            output = main(filename)
        print(output)
    elif port > 0:
        launch_server(port, host)
