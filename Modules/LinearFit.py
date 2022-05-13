import matplotlib.pyplot as plt
import numpy as np
import scipy.odr as odr
from scipy.optimize import least_squares, curve_fit

def fitFunc(p, x):
    '''
    Fit function
    '''
    f= p[0]*x+p[1]

    return f

def fitFuncDiff(p, x):
    '''
    Differential of fit function
    '''
    df= p[0]
    return df

def calcChiSq(p, x, y, xerr, yerr):
    '''
    Error function for fit
    '''
    e = (y - fitFunc(p, x))/(np.sqrt(yerr**2 + fitFuncDiff(p, x)**2*xerr**2))
    return e

def fitStdError(jacMatrix):

    # Compute covariance
    jMat2 = np.dot(jacMatrix.T, jacMatrix)
    detJmat2 = np.linalg.det(jMat2)

    # Prepare output
    output = np.zeros(jMat2.shape[0])
    if detJmat2 < 1E-32:
        print("Value of determinat detJmat2",detJmat2)
        print("Matrix singular, error calculation failed.")
        return output
    else:
        covar = np.linalg.inv(jMat2)
        for i in range(len(output)):
            output[i] = np.sqrt(covar[i, i])

        return output
        
class LinearFit:
    def __init__(self, xdata, ydata, yerror , xerror, m, c):
        self.xdata = xdata
        self.ydata = ydata
        self.xerror = xerror
        self.yerror = yerror
        self.m = m
        self.c = c
    

    
    def fit(self):
        initParams = np.array([self.m,self.c])
        pInit = initParams
        nPoints = len(self.xdata)
        nPars = len(initParams)

        # Run fit
        output = least_squares(calcChiSq, pInit, args = (self.xdata, self.ydata, self.xerror, self.yerror))
        # Get least_squares output, stored in array output.x[]
        self.grad = output.x[0]
        self.yintercept = output.x[1]
        # Get errors from our fits using fitStdError(), defined above
        pErrors = fitStdError(output.jac)
        self.gradErr = pErrors[0]
        self.yinterceptErr = pErrors[1]

        # Calculate fitted y-values using our fit parameters and the original fit function
        xPlot = np.linspace(np.min(self.xdata), np.max(self.xdata), 300)
        fitData = fitFunc(output.x, xPlot)

        # Calculate chis**2 per point, summed chi**2 and chi**2/NDF
        chiarr = calcChiSq(output.x, self.xdata, self.ydata, self.xerror, self.yerror)**2
        chisq = np.sum(chiarr)
        NDF = nPoints - nPars
        chisqndf = chisq/NDF
        self.fitData = fitData
        self.xPlot = xPlot

        ##### needs a lot of work ####
    def get_xPlot(self):
        return self.xPlot
    
    def get_fitData(self):
        return self.fitData

    def get_grad(self):
        return self.grad
    
    