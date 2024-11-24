import numpy as np

def test2(a, b):
    print(a, b, sep="|")

def test(*args):
    test2(*args)

np.array([1])[0]
test("aboba", "test")