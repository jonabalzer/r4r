import unittest
import bsplines as bs
import numpy as np

class TestSplineCurve(unittest.TestCase):

    def setUp(self):

        # FIXME: Approximation of circle better done with NURBS. Replace this
        # with something else. 
        knots = bs.knot_vector(2,True)
        knots.make_uniform(0,2*np.pi,150)
        self.circle = bs.curve(knots,2)
        
        t = np.linspace(0,2*np.pi,200)
        x = np.cos(t)+10
        y = np.sin(t)+10
        self.circle.interpolate(t,[x,y])
        
        knots = bs.knot_vector(3,False)
        knots.make_uniform(0,1.0,10)
        self.line = bs.curve(knots,2)
        
        t = np.linspace(0,1.0,20)
        self.line.interpolate(t,[t,t])

    def test_area(self):
        
        A = self.circle.area(2)
        print A
        self.assertTrue(np.abs(A-np.pi)<0.2)

    def test_barycenter(self):
        x,y,A = self.circle.barycenter(2)
        print x,y,A
        self.assertTrue(np.abs(x-10)<0.1 and np.abs(y-10)<0.1)
        
        
    def test_integration(self):

        f = lambda t:t                 
        res = self.line.integrate(f,2)
        self.assertTrue(np.abs(res-np.sqrt(2)*0.5)<0.01)
        

if __name__ == '__main__':
    try:
        unittest.main()
    except SystemExit:
        print "Finished testing."