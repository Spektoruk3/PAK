import unittest
import FluidCube as fc
import numpy as np

class TestFluidCube(unittest.TestCase):
    def testFluid(self):

        model = fc.FluidCube(0, 1e-7, 1e-2)
        self.AssertStep(model)
        density = model.density
        self.NanInf(density)

        model = fc.FluidCube(100000000, 1e-7, 1e-2)
        self.AssertStep(model)
        density = model.density
        self.NanInf(density)

        model = fc.FluidCube(1e-9 / 2, 1e-7, 1e-2)
        model.density[1, 1] = 1000
        model.density[2, 2] = 0
        self.AssertStep(model)
        density = model.density
        self.NanInf(density)

        self.assertLess(model.density[1, 1], 255)
        self.assertLess(model.density[2, 2], 255)
        return True

    def NanInf(self, density):
        error1 = np.all(not list(np.isnan(density)))
        error2 = np.all(not list(np.isinf(density)))
        self.assertEqual(error1, False)
        self.assertEqual(error2, False)
        return True

    def AssertStep(self, model):
        try:
            model.FluidCubeStep()
        except:
            print("An error occurred during the step function!")
            return False
        return True

def run_tests():
     unittest.main()

if __name__ == '__main__':
    run_tests()