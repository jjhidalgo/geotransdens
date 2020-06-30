"""
Functions to call PyGimli from geotransdens
authors: AAP, JHG
06-2020 
"""

import numpy as np
import pygimli as pg 
import pybert as pb
import time
import csv
from pygimli.utils import gmat2numpy


global meshERT, meshTD
meshERT = pg.load('ERTMESH.bms')
start = time.time()
meshTD = pg.load('HYDROMESH.vtk')
end = time.time()
#print("Time to load the VTK file and create TD mesh: "+str(end - start)+"seconds")

def interpolant():
    ''' Creates interpolant between ERT mesh and TD mesh'''
    start = time.time()
    interpMat = meshERT.interpolationMatrix(meshTD.positions())
    end = time.time() 
    print("Time to create interpolant  " + str(end - start) + " seconds")
    return interpMat

def interpolate_jacobian(inVector, interpolant=None):
    ''' Interpolates a vector defined in the ERT mesh in the TD mesh nodes.'''    
    if interpolant is None:
        outdata = np.array(meshERT.interpolationMatrix(meshTD.positions()) * pg.core.cellDataToPointData(meshERT, inVector))
    else:
        outdata = np.array(interpolant *  pg.core.cellDataToPointData(meshERT, inVector))
    return outdata

def findERTscheme(time):
    ''' Selects the ERT scheme for the given time from the auxiliary file ERTSCHEMES.dat''' 
    initTime = [] 
    endTime = []
    fileName = []
    with open('ERTSCHEMES.dat') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
                continue
            else:
                initTime.append(int(row[0]))
                endTime.append(int(row[1]))
                fileName.append(row[2])  
                line_count += 1
    
    for i in range(len(initTime)):
        if initTime[i] <= time <= endTime[i]:
            timeFound = i
        else:
            continue
    
    return fileName[timeFound]
    
    
    
def RUN_ERT(concentration, porosity, dc_dp, sim_time):
    ''' Computes resisitivity (rhoa) and the derivative of the resisitivity to the
        calibrated parameters (drhoa_dp)
        drhoa_dp = drhoa_dc * dc_dp'''    

    start_total = time.time()

    cell_conc = pg.meshtools.nodeDataToCellData(meshTD, concentration)
    
    #FROM SALT MASS FRACTION TO SALT CONCENTRATION TO WATER CONDUCTIVITY
    csw = 5.5 #[S/m]
    cfw = 0.005 #[mS/m]
    #Normalized salt_mass fraction for conversion to water EC. 
    salt_mass = cell_conc/max(cell_conc)
    conductivity = (salt_mass * csw) + (1-salt_mass)*cfw


    #PETROPHYSICS: FROM WATER CONDUCTIVITY TO BULK RESISTIVITY
    m = 1.3
    resistivity_bulk = 1./conductivity * porosity**(-m)
    rhob = pg.core.RVector(meshERT.cellCount())
    pg.interpolate(meshTD, resistivity_bulk, meshERT.cellCenters(), rhob)
    rhob = pg.meshtools.fillEmptyToCellArray(meshERT, rhob)
    
    #SIMULATION OF ERT ACQUISITION USING BULK RESISTIVITY MODEL
    print('Acquisition...')
    start_acq = time.time()
    smFileName = findERTscheme(sim_time)
    sm = pb.DataContainerERT(smFileName)
    fERT = pb.DCSRMultiElectrodeModelling(meshERT, sm)
    RM = fERT.regionManager()
    RM.region(0).setBackground(True)
    rhoa = fERT.response(rhob)
    end_acq = time.time()
    print("Time to create acquisition: " + str(end_acq - start_acq) + " seconds")
    
    start_createJ = time.time()  
    fERT.createJacobian(rhob)
    end_createJ = time.time()
    print("Time to create the jacobian: " + str(end_createJ - start_createJ) + " seconds")
    
    start_convertJ = time.time()  
    jacobian = gmat2numpy(fERT.jacobian())
    end_convertJ = time.time()
    print("Time to convert jacobian to numpy array " + str(end_convertJ - start_convertJ) + " seconds")
    
    interpMat = interpolant()
    
    start_interp = time.time()
    outdata = np.apply_along_axis(interpolate_jacobian, 1, jacobian,interpolant=interpMat)
    end_interp = time.time()
    print("Time to interpolate jacobian by axes (interpolant as arg): "+str(end_interp - start_interp)+"seconds")
    
    drhoa_drhob = np.array(outdata)
    drhoa_drhob = drhoa_drhob.reshape(sm.size(), meshTD.nodeCount())  
    
    #ANALYTICAL COMPUTATION TO ACCOUNT FOR PETROPHYSICAL TRANSFORMATION. 
    node_por = pg.meshtools.cellDataToNodeData(meshTD, np.array(porosity))

    for i in range(sm.size()):
        drhoa_drhob[i,:] = drhoa_drhob[i,:] * (1/(np.power(node_por, -(1)*m))) 
    for i in range(sm.size()):
        drhoa_drhob[i,:] = drhoa_drhob[i,:] * (csw - cfw)
    
    #TO DO: check if this below works better
    #drhoa_drhob = np.multiply(drhoa_drhob, (csw - cfw)/np.power(node_por, -m))

    drhoa_dp = drhoa_drhob.dot(dc_dp)
    
    end_total = time.time()
    print("Total simulation time: " + str((end_total-start_total)/60) + ' minutes')

    return np.array(rhoa), drhoa_dp
    
 
