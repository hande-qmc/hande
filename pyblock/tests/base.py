import unittest
import numpy

def corr_data(N, L, seed):
    '''

See: http://stackoverflow.com/questions/20357854/how-can-i-generate-gaussian-random-process-using-matlab.
    '''
    # Use identical seed so tests are reproducible.
    rand = numpy.random.RandomState(seed=seed)
    return numpy.convolve(rand.randn(2**N), numpy.ones(2**L)/10, 'same')

# produced from reblocking corr_data(10, 6, 5)
reblock_1D = [
        (0, 1024, numpy.array(0.20848489810011175), numpy.array(0.5059733357706058), numpy.array(0.022228688348921586), numpy.array(0.0004914292726987686)),
        (1, 512, numpy.array(0.20848489810011178), numpy.array(0.5015874592634239), numpy.array(0.0312995687889446), numpy.array(0.0009790681131713157)),
        (2, 256, numpy.array(0.20848489810011184), numpy.array(0.4956240227095556), numpy.array(0.04400035612025432), numpy.array(0.001948368517458123)),
        (3, 128, numpy.array(0.20848489810011173), numpy.array(0.48304884454343255), numpy.array(0.06143141784132584), numpy.array(0.003854549974240774)),
        (4, 64, numpy.array(0.2084848981001117), numpy.array(0.46303627895395055), numpy.array(0.08505846141716576), numpy.array(0.007577610011170585)),
        (5, 32, numpy.array(0.20848489810011178), numpy.array(0.40641824428237183), numpy.array(0.11269680622725792), numpy.array(0.014312508703377615)),
        (6, 16, numpy.array(0.20848489810011178), numpy.array(0.3388777084790293), numpy.array(0.14553300924511708), numpy.array(0.026570570675052838)),
        (7, 8, numpy.array(0.20848489810011178), numpy.array(0.228782900550399), numpy.array(0.16910902568698064), numpy.array(0.04519628822370252)),
        (8, 4, numpy.array(0.20848489810011178), numpy.array(0.21072809492758346), numpy.array(0.229525649398702), numpy.array(0.09370345398462808)),
        (9, 2, numpy.array(0.20848489810011175), numpy.array(0.011498881312406542), numpy.array(0.07582506614704185), numpy.array(0.05361641845649181))
    ]
reblock_1D_opt = [9]

# produced from reblocking [corr_data(10, 6, 5), corr_data(10, 6, 7)]
reblock_2D = [
        (0, 1024, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.50597334,  0.02776897], [ 0.02776897,  0.61260927]]), numpy.array([ 0.02222869,  0.02445917]), numpy.array([ 0.00049143,  0.00054074])),
        (1, 512, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.50158746,  0.02837632], [ 0.02837632,  0.60851945]]), numpy.array([ 0.03129957,  0.03447484]), numpy.array([ 0.00097907,  0.00107839])),
        (2, 256, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.49562402,  0.02892238], [ 0.02892238,  0.60270839]]), numpy.array([ 0.04400036,  0.04852143]), numpy.array([ 0.00194837,  0.00214857])),
        (3, 128, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.48304884,  0.02940551], [ 0.02940551,  0.59443823]]), numpy.array([ 0.06143142,  0.06814726]), numpy.array([ 0.00385455,  0.00427594])),
        (4, 64, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.46303628,  0.02918235], [ 0.02918235,  0.57765284]]), numpy.array([ 0.08505846,  0.09500434]), numpy.array([ 0.00757761,  0.00846366])),
        (5, 32, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.40641824,  0.03257275], [ 0.03257275,  0.5392826 ]]), numpy.array([ 0.11269681,  0.12981749]), numpy.array([ 0.01431251,  0.01648684])),
        (6, 16, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.33887771, -0.00386415], [-0.00386415,  0.45358458]]), numpy.array([ 0.14553301,  0.16837172]), numpy.array([ 0.02657057,  0.03074033])),
        (7, 8, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.2287829 ,  0.08372261], [ 0.08372261,  0.35810369]]), numpy.array([ 0.16910903,  0.21157259]), numpy.array([ 0.04519629,  0.05654515])),
        (8, 4, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.21072809,  0.12820696], [ 0.12820696,  0.20183304]]), numpy.array([ 0.22952565,  0.22462916]), numpy.array([ 0.09370345,  0.09170447])),
        (9, 2, numpy.array([ 0.2084849, -0.1791132]), numpy.array([[ 0.01149888,  0.04086769], [ 0.04086769,  0.14524615]]), numpy.array([ 0.07582507,  0.26948669]), numpy.array([ 0.05361642,  0.19055586]))
    ]
reblock_2D_opt = [9,8]
