# Author: Yang Long <longyang_123@yeah.net>
#
# License: LGPL-2.1
import matplotlib.pyplot as plt
import numpy as np

class Viewer:
    def __init__(self,bpm):
        self.bpm = bpm
    
    def viewintensity(self):
        field = self.bpm.output()
        plt.imshow(np.abs(field),cmap=plt.cm.hot)
        plt.colorbar()
        plt.show()

    def viewphase(self):
        field = self.bpm.output()
        plt.imshow(np.angle(field))
        plt.colorbar()
        plt.show()
    
    def viewrecord(self,t):
        record = self.bpm.info()
        plt.imshow(np.abs(record[t]),cmap=plt.cm.hot)
        plt.colorbar()
        plt.show()

        
