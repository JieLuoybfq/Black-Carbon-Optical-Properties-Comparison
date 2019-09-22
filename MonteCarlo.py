from ConfigParserModule import logging
import numpy as np
from decimal import Decimal


class MonteCarloMethods:

    def __init__(self):
        pass

    def createRandomNormalArr(self, center, width, size, decimalOutput=True, digit=3):
        try:
            A = np.random.normal(center, width, int(size))
            length = len(A)
            B = []
            if decimalOutput:
                for i in range(length):
                    B.append(round(Decimal(A[i]), digit))
            else:
                for i in range(length):
                    B.append(round(A[i], digit))
            return B

        except Exception as e:
            logging.exception(e)
            raise

    def getRandomFromArr(self, Array, Number):
        try:
            # uniform choice
            A = np.random.choice(Array, size=int(Number), replace=False)
            return A

        except Exception as e:
            logging.exception(e)
            raise
