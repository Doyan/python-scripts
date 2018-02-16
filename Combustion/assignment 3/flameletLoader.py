# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:41:00 2017

@author: anygren
"""
import numpy as np
import cantera as ct
class loadFlamelet:
    def __init__(self, flameletLib):
        # Constructor - Loads a fluent flamelet library into python and sorts each flamelet into a python dictionary
        self.flamelets = dict()
        flameletScalars = dict()
        self.scalarDissRates = np.array([])
            
        f = open(flameletLib, 'r')
        lines = f.readlines()
        f.close()
        gas = ct.Solution('gri30.xml')
        
        Y_names = gas.species_names
        gridpoints = []
        for i in range(len(Y_names)): 
            Y_names[i] = 'massfraction-' + Y_names[i].lower()
        flameletScalarNames = ['Z','TEMPERATURE'] + Y_names
        for line in lines:
            if 'GRIDPOINTS' in line:
                gpSTR = line[11:-1]
                gridpoints.append(int(gpSTR))                
            elif 'CHI' in line:
                chiSTR = line[4:-1]
                self.scalarDissRates = np.append(self.scalarDissRates, float(chiSTR))
                
        flameletNo = 0
        ScalarTemp = []
        for i in range(len(lines)):
            if lines[i] == 'HEADER\n':
                flameletNo += 1
                flameletScalars = dict()
            if lines[i].strip('\n') in flameletScalarNames:
                scalarName = lines[i]
                if gridpoints[flameletNo - 1] % 5 == 0:
                    ScalarTemp = np.loadtxt(lines[(i + 1):(i + 1 + int(gridpoints[flameletNo - 1]/5))])
                    ScalarTemp = ScalarTemp.flatten()
                    flameletScalars[scalarName.strip('\n')] = ScalarTemp
                    self.flamelets[str(flameletNo)] = flameletScalars
                else:
                    ScalarTemp = np.loadtxt(lines[(i + 1):(i + 1 + int(gridpoints[flameletNo - 1]/5))])
                    ScalarTemp = np.append(ScalarTemp, np.loadtxt(lines[(i + 1 + int(gridpoints[flameletNo - 1]/5)):i + 2 + int(gridpoints[flameletNo - 1]/5)]))
                    ScalarTemp = ScalarTemp.flatten()
                    flameletScalars[scalarName.strip('\n')] = ScalarTemp         
                    self.flamelets[str(flameletNo)] = flameletScalars
    

        self.nFlamelets = flameletNo
        self.nGridPoints = gridpoints

          
    def getFlameletScalar(self, flameletNo, flameletScalar):
	# Access a Scalar in mixture fraction Space for a certain flamelet
        flamelet = self.flamelets[str(flameletNo)]
        scalar = flamelet[flameletScalar]
        return scalar
        
    def setFlameletScalar(self, flameletNo, flameletScalarName, flameletScalarArray):
	# Add a computed scalar to the flamelet library (For Instance Density)
        flamelet = self.flamelets[str(flameletNo)]
        flamelet[flameletScalarName] = flameletScalarArray
        return 0

    def getFlameletChi(self, flameletNo):
        Chi = self.scalarDissRates[flameletNo - 1]
        return Chi
    def interpolateFlamelet(self, chi, scalar, nInterpPoints = 50):
        # Takes any Value for chi and Linearly interpolate a flamelet between 
        # two existing flamelets in the library for any flameletscalar
        # Optional - Set the number of coordinates to be interpolated
        
        # Find closest chi in i flamelet to the Chi we want to interpolate
        chiIndex1 = (np.abs(self.scalarDissRates - chi)).argmin()
            
        if chiIndex1 == self.scalarDissRates[-1]:
            chiIndex2 = chiIndex1 - 1
        else:
            chiIndex2 = chiIndex1 + 1
            
        # First Interpolate Inside each flamelet due to different gridsize for each flamelet
        Zgrid = np.linspace(0, 1, nInterpPoints)
        scalarGrid1 = np.interp(Zgrid, self.getFlameletScalar(chiIndex1 + 1, 'Z'), self.getFlameletScalar(chiIndex1 + 1, scalar))            
        scalarGrid2 = np.interp(Zgrid, self.getFlameletScalar(chiIndex2 + 1, 'Z'), self.getFlameletScalar(chiIndex2 + 1, scalar))            
        
        # Then Interpolate between the two flamelets
        interpolatedScalar = np.zeros(np.shape(Zgrid))        
        for i in range(len(Zgrid)):
            interpolatedScalar[i] = np.interp(chi, np.array([self.scalarDissRates[chiIndex1], self.scalarDissRates[chiIndex2]]), np.array([scalarGrid1[i], scalarGrid2[i]]))
                       
        return interpolatedScalar
