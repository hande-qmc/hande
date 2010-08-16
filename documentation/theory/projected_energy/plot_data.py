#!/usr/bin/python

import numpy
import pylab

def plot_graph(data,x,y,bars,fci,title,plotfile):
    pylab.cla()
    pylab.errorbar(data[x],data[y],yerr=data[bars])
    pylab.axhline(fci, color='r')
    pylab.xlabel(r'$\langle N_w \rangle$')
    pylab.ylabel(r'$E$')
    pylab.title(title)
    pylab.savefig(plotfile)

if __name__ == '__main__':
    data = numpy.transpose(numpy.loadtxt('blocking.data', usecols=(2,3,5,6,8,9,11,12,14,15)))

    fci = data[-1][0]
    plot_graph(data,0,1,2,fci,'Shift','shift.pdf')
    plot_graph(data,0,3,4,fci,r'proj.e $\langle N_j/N_0\rangle$','proje1.pdf')
    plot_graph(data,0,5,6,fci,r'proj.e $\langle N_j^\prime/N_0^\prime\rangle$','proje2.pdf')
    plot_graph(data,0,7,8,fci,r'proj.e $\langle N_j^\prime\rangle/\langle N_0^\prime\rangle$','proje3.pdf')
