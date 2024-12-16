import numpy as np
import math

class ProbabilityGenerator:
    def calculateProbabilityList(self, listLen):
        print("Base class: nothing")
        return []

class LinearDecrease(ProbabilityGenerator):
    def calculateProbabilityList(self, listLen):
        total = sum(range(1, listLen + 1))
        probList = [(listLen - i) / total for i in range(listLen)]
        return probList

class LinearIncrease(ProbabilityGenerator):
    def calculateProbabilityList(self, listLen):
        total = sum(range(1, listLen + 1))
        probList = [(i + 1) / total for i in range(listLen)]
        return probList

class Uniform(ProbabilityGenerator):
    def calculateProbabilityList(self, listLen):
        probList = [1/listLen for _ in range(listLen)]
        return probList

class GaussianDistribution(ProbabilityGenerator):
    def calculateProbabilityList(self, listLen):
        mean, std_dev = self.getParams()["mean"], self.getParams()["stdDev"]
        values = [(1/(std_dev * np.sqrt(2*np.pi)))*np.exp(-0.5*np.power((x-mean)/std_dev, 2)) for x in range(listLen)]
        total = np.sum(values)
        probList = values / total
        return probList

    def getParams(self):
        return {"mean": 0, "stdDev": 1}

class GaussianDistributionM0SD1(GaussianDistribution):
    def getParams(self):
        return {"mean": 500, "stdDev": 1000}

class ParetoDistribution(ProbabilityGenerator):
    def calculateProbabilityList(self, listLen):
        xm, k = self.getParams(listLen)["xm"], self.getParams(listLen)["k"]
        values = []
        for x in range(1, listLen+1):
            P = 1-np.power(xm/x, k)
            if P <= 0:
                P = 0
            values.append(P)
        total = np.sum(values)
        probList = values / total
        return probList

    def getParams(self, len):
        return {"xm": len*0.2, "k": 2}

class ParetoDistributionXm20K3(ParetoDistribution):
    def getParams(self, len):
        return {"xm": len*0.2, "k": 3}

Methods = {
    'linear-decrease': LinearDecrease(),
    'linear-increase': LinearIncrease(),
    'uniform': Uniform(),
    'gaussia mean=500 std_dev=1000': GaussianDistributionM0SD1(),
    'pareto xm=20% k=3': ParetoDistributionXm20K3(),
}

def GetProbabilities(num, method='linear-decrease'):
    generator = Methods[method]
    return np.array(generator.calculateProbabilityList(num))

def GetAllMethods():
    return Methods.keys()
