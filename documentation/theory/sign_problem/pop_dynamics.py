#!/usr/bin/python

# Plot the solution to 
#  n'(t) = p e^{2 b t} + q n(t) + r n(t)^2
# Setting f = r n / q and s = q t shows that there are actually only two
# meaningful parameters:
#  f'(s) = a e^{2 c s} + f(s) + f(s)^2
# We choose to solve this instead, where
# a = p r / q^2 and c = b / q.

# The differential equation is a Riccati differential equation with the
# solution given by
#  f(z) = -2 c z u'(z) / u(z)
# where
#  z(s) = - a e^{2 c s} / 4 c^2
#  u(z) = k _0F_1(;1-1/(2c);z) + z^{1/(2c)} _0F_1(;1+1/(2c);z)
#  u'(z) = 2kc/(2c-1) _0F_1(;2-1/(2c);z) + z^{1/(2c)-1}/(2c) _0F_1(;1+1/(2c);z)
#         + 2c z^{1/(2c)}/(2c+1) _0F_1(;2+1/(2c);z)
# and k is the constant of integration which we shall set such that n(t=0) = 0.

# WARNING: this assumes that 1/2c is not an integer.

import scipy
import scipy.special

import pylab
import matplotlib.widgets

def z(a, c, s):
    return - a*scipy.exp(2*c*s) / (4*c**2)

def u(k, c, z):
    ic2 = 1.0/(2*c)
    f01 = scipy.special.hyp0f1
    return k*f01(1-ic2,z) + z**(ic2)*f01(1+ic2,z)

def dudz(k, c, z):
    c2 = 2.0*c
    ic2 = 1.0/c2
    f01 = scipy.special.hyp0f1
    return k*c2/(c2-1)*f01(2-ic2,z) + ic2*z**(ic2-1)*f01(1+ic2,z) + c2*z**ic2/(c2+1)*f01(2+ic2,z)

def set_boundary_condition(a, c):
    c2 = 2.0*c
    ic2 = 1.0/c2
    f01 = scipy.special.hyp0f1
    z0 = z(a, c, 0.0)
    return -(ic2*z0**(ic2-1)*f01(1+ic2,z0) + c2*z0**ic2/(c2+1)*f01(2+ic2,z0))*(c2-1) / (c2*f01(2-ic2,z0))

def n(a, c, s):
    k = set_boundary_condition(a, c)
    z_values = z(a, c, s)
    return -2*c*z_values*dudz(k, c, z_values) / u(k, c, z_values)

def plot_dynamics(a = -0.006, c = 0.03, smin = 0, smax = 100, sstep = 0.01):

    s = scipy.arange(smin, smax, sstep)

    # main plot area
    ax = pylab.subplot(111)
    pylab.subplots_adjust(left=0.10, bottom=0.30, right=0.90)
    ax.set_xlabel('$t$')
    ax.set_ylabel('$n(t)$')
    ax.set_title("$n^{'}(t) = a e^{2 c t} + n(t) + n(t)^2$")

    # plot
    # Flip if a<0
    if a < 0:
        (line,) = pylab.plot(s, -n(a,c,s))
    else:
        (line,) = pylab.plot(s, n(a,c,s))

    # sliders
    axcolor = 'lightgoldenrodyellow'
    aax = pylab.axes([0.1, 0.1, 0.65, 0.03], axisbg=axcolor)
    cax = pylab.axes([0.1, 0.15, 0.65, 0.03], axisbg=axcolor)
    aslider = matplotlib.widgets.Slider(aax, '$a$', -0.04, 0, valinit=a, valfmt="%1.3f")
    cslider = matplotlib.widgets.Slider(cax, '$c$', 0, 0.3, valinit=c, valfmt="%1.3f")
    # button
    resetax = pylab.axes([0.8, 0.025, 0.1, 0.04])
    button = matplotlib.widgets.Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

    # callback functions
    def update(val):
        if aslider.val < 0:
            line.set_ydata(-n(aslider.val, cslider.val, s))
        else:
            line.set_ydata(n(aslider.val, cslider.val, s))
        ax.set_ylim([0.0,abs(n(aslider.val, cslider.val, smax))])
        pylab.draw()
        return None
    def reset(event):
        aslider.reset()
        cslider.reset()
        return None

    # call callbacks on event
    aslider.on_changed(update)
    cslider.on_changed(update)
    button.on_clicked(reset)

    pylab.show()

    return None

if __name__ == '__main__':

    plot_dynamics()
