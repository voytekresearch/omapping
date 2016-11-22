"""Create fake data for testing OM"""

from __future__ import print_function, division

import sys
import numpy as np

# Import FOOF (use sys to add location to path, then import)
sys.path.append('/Users/thomasdonoghue/Documents/GitCode/omegamappin/')
from om.gen import save_foof_pickle
from om.tests.helper_test_funcs import TestDB


def make_fake_foof_dat_1():
    """   """

    db = TestDB()

    v1 = (1.0, np.array([5, 10]), np.array([1.0e-22, 1.0e-22]), np.array([1.0, 1.0]))
    v2 = (1.0, np.array([5, 10, 15]), np.array([1.0e-22, 1.0e-22, 1.0e-22]), np.array([1.0, 1.0, 1.0]))

    fake_foof_dat = [v1, v2]

    save_foof_pickle(fake_foof_dat, db.foof_path, 'test1')

def make_fake_foof_dat_2():
    """   """

    db = TestDB()

    v1 = (0.0, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))
    v2 = (0.5, np.array([10, 11, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))
    v3 = (1.0, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))
    v4 = (1.5, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))
    v5 = (2.0, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))

    fake_foof_dat = [v1, v2, v3, v4, v5]

    save_foof_pickle(fake_foof_dat, db.foof_path, 'test2')


if __name__ == "__main__":
    make_fake_foof_dat_1()
    make_fake_foof_dat_2()
    print("Testing data created.")