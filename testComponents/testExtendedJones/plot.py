"""
Copyright (C) 2015 Zong-han, Xie <icbm0926@gmail.com>.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt

myfile="CrossPolarizer_MultiWaveLength_TestLightSrc_Lambertian.csv"
answerfile="Answer_CrossPolarizer_MultiWaveLength_TestLightSrc_Lambertian.csv"

myfiledata = np.loadtxt(myfile,delimiter=',')
answerfiledata = np.loadtxt(answerfile,delimiter=',')

plotphidegree = 60.0
plotxdata = np.array(range(0,81))

myfileplotdata = np.zeros(81)
index = 0
for i in myfiledata:
    if (abs(i[1] - plotphidegree) < 1.0e-10):
        myfileplotdata[index] = i[2]
        index = index + 1

index = 0
answerfileplotdata = np.zeros(81)
for i in answerfiledata:
    if (abs(i[1] - plotphidegree) < 1.0e-10):
        answerfileplotdata[index] = i[2]
        index = index + 1

plt.plot(plotxdata, myfileplotdata, 'b')
plt.plot(plotxdata, answerfileplotdata, 'r')
plt.show()
