import numpy as np
import matplotlib.pyplot as plt

class Histogram:

    def __init__(self, numberOfBins, binBot, binTop, name = "histo"):
        self._numberOfBins = int(numberOfBins)
        self._binBot = binBot
        self._binTop = binTop
        self._name = name

        self._binWidth = (self._binTop - self._binBot) / self._numberOfBins
        self._sumWeights = np.zeros(self._numberOfBins)
        self._sumWeights2 = np.zeros(self._numberOfBins)
        self._underflow = 0
        self._overflow = 0
        self._nfilled = 0

        # calculate and save the x coordinate of the centre of each bin
        self._binCenters = np.array([self._binBot + (i+0.5)*self._binWidth for i in range(self._numberOfBins)])
        
    def fill(self, value, weight = 1.0):
        self._value = value
        if value < self._binBot:
            self._underflow += 1
        elif value >= self._binTop:
            self._overflow += 1
        else:
            # add weight to the correct bin
            ibin = int( (value - self._binBot)/self._binWidth)
            self._sumWeights[ibin] +=  weight
            self._sumWeights2[ibin] += weight**2
            
        self._nfilled += 1

    def get_bin_errors(self):
        # returns all bin errors
        return np.sqrt(self._sumWeights2)

    def get_bin_content(self, nbin):
        # returns the contents on bin 'nbin'
        return self._sumWeights[nbin]

    def get_bin_error(self, nbin):
        # returns the error on bin 'nbin'
        return self.get_bin_errors()[nbin]

    def __str__(self):
        output = "Histogram " + self._name + "\n"
        for i in range (self._numberOfBins):
            output += f"Bin {i} = {self.get_bin_content(i)}  +- {self.get_bin_error(i)}\n"
        output += f"The number of fills was {self._nfilled}\n"
        output += f"Underflows = {self._underflow}, Overflows = {self._overflow}\n"
        return output

    def plot(self, xlabel = 'x', ylabel = 'N', color = 'r', label = 'Data', marker = '+', filename = 'unnamedhistogram.png'):
        plt.figure(figsize = (8, 6))
        plt.title(self._name)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.errorbar(self._binCenters, self._sumWeights, 
                     yerr = self.get_bin_errors(),
                     xerr = [self._binWidth/2]*self._numberOfBins,
                     color = color, marker = marker, linestyle = '', label = label)
        #plt.bar(self._binCenters, self._sumWeights, width=[self._binWidth]*self._numberOfBins)
        plt.legend()
        plt.savefig(filename)
        return plt
