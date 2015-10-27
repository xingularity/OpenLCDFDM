#########################################################################
## Copyright (C) 2015 Zong-han, Xie <icbm0926@gmail.com>.
## All rights reserved.
##
## You may use this file under the terms of the BSD license as follows:
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##   * Neither the name of OpenLCDFDM nor the names of its contributors
##     may be used to endorse or promote products derived from this
##     software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
##
#########################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lcd1d import *

def readLightSourceSpectrum(fname=''):
    if (len(fname) == 0):
        return{}
    temp = pd.read_csv(fname, sep=',',header=None).values
    answer={}
    for i in range(temp.shape[0]):
        answer[temp[i,0]]=temp[i,1]
    return answer

def readIsotropicSpectrum(fname=''):
    if (len(fname) == 0):
        return{}
    temp = pd.read_csv(fname, sep=',',header=None).values
    answer={}
    for i in range(temp.shape[0]):
        answer[temp[i,0]]=temp[i,1]+1j*temp[i,2]
    return answer

def writeTransmissions(fname, data, inAngles):
    f = open(fname, 'w')
    if (len(data.shape) != 2):
        raise Exception("data has wrong dimension numbers.")
    if (len(inAngles.shape) != 3):
        raise Exception("inAngles has wrong dimension numbers.")
    if (data.shape != inAngles.shape[0:2]):
        raise Exception("Shapes of data and inAngles dooen't match.")
    for i in range(inAngles.shape[0]):
        for j in range(inAngles.shape[1]):
            f.write(str(inAngles[i, j, 0]*180.0/np.pi)+", "+str(inAngles[i, j, 1]*180.0/np.pi)+", "+str(data[i, j])+"\n")
    f.close()

def writeNormalTrans(fname, data):
    f = open(fname, 'w')
    if (len(data.shape) != 1):
        raise Exception("data has wrong dimension numbers.")
    for i in range(data.shape[0]):
        f.write(str(data[i])+"\n")
    f.close()

def getAnglesForPolarPlot(inAngles):
    azimuths = set()
    zeniths = set()
    for i in range(inAngles.shape[0]):
        for j in range(inAngles.shape[1]):
            zeniths.add(inAngles[i, j, 0]*180.0/np.pi)
            azimuths.add(inAngles[i, j, 1]*180.0/np.pi)
    azimuths = np.radians(list(azimuths))
    zeniths = np.array(list(zeniths))
    r, theta = np.meshgrid(zeniths, azimuths)
    return (theta, r)

def glass_only():
    lcd1dstaticmain = pyLCD1DStaticMain()
    nk = readIsotropicSpectrum('TestGlassNKSpectrum.csv')
    lcd1dstaticmain.addOpticalGlassLayer(20.0, nk)
    lcd1dstaticmain.setOMPThreadNum(1)
    lcd1dstaticmain.setOpticalIncidentAngles(1,1)
    lcd1dstaticmain.setOpticalWavelength(0.55)
    lcd1dstaticmain.createExtendedJones()
    lcd1dstaticmain.calculate()
    transmissions = np.array(lcd1dstaticmain.getTransmissions())
    transmissions = transmissions[0]
    inAngles = np.array(lcd1dstaticmain.getIncidentAngles())
    writeTransmissions("GlassOnly_550nm.csv", transmissions, inAngles)
    lcd1dstaticmain.setOpticalWavelength(0.38, 0.78, 0.01)
    lcd1dstaticmain.setOpticalSourceSpectrum(readLightSourceSpectrum('TestLightSrc.csv'))
    lcd1dstaticmain.createExtendedJones()
    lcd1dstaticmain.calculate()
    transmissions = np.array(lcd1dstaticmain.getTransmissions())
    transmissions = transmissions[0]
    writeTransmissions("GlassOnly_MutiWavelength.csv", transmissions, inAngles)
    lcd1dstaticmain.useOptical2X2Lambertian();
    lcd1dstaticmain.createExtendedJones()
    lcd1dstaticmain.calculate()
    transmissions = np.array(lcd1dstaticmain.getTransmissions())
    transmissions = transmissions[0]
    #plot transmissions
    theta, r = getAnglesForPolarPlot(inAngles)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    #fig, ax = plt.subplots()
    ax.set_xlabel("phi(degree)")
    ax.set_ylabel("theta(degree)")
    cs = ax.contourf(theta, r, transmissions.T, 255)
    plt.colorbar(cs)
    plt.show()


def main():
    glass_only()


if __name__ == '__main__':
    main()
