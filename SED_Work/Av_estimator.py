import numpy as np
import matplotlib.pyplot as plt

def Linenum(infile):
    """Calculates the number
    of lines in infile and return it"""
    with open(infile,'r') as inputFile:
        for i,line in enumerate(inputFile):
            pass
    nLines = i
#    print 'Numero de lineas:',nLines
    return nLines

path = 'path_where_my_catalogs_and_tables_are'

infile_stars = path+'my_catalog.txt'
infile_KH = path+'KH95_tableA5.txt'

### Calculates the number of lines in files
Lines1 = Linenum(infile_stars)
Lines2 = Linenum(infile_KH)

### define variables and load values from tables
Spt_KH,Teff_KH,B_V_KH = [],[],[]

with open(infile_KH,'r') as inputFile:
    for j in range(Lines2+1):
            
        line = inputFile.readline()
        ## I put a # symbol at the begining of my tables
        ## so the columns labels are ignored
        if line.find('#') < 0:
            
            linelist = line.split()
            Spt_KH.append(linelist[1].strip())
            Teff_KH.append(float(linelist[2].strip()))
            B_V_KH.append(float(linelist[5].strip()))

### Do the same for your catalog

Spt,Teff,B_V = [],[],[]

with open(infile_stars,'r') as inputFile:
    for j in range(Lines1+1):
            
        line = inputFile.readline()
        if line.find('#') < 0:
            
            linelist = line.split()
            Spt.append(linelist[1].strip())
            Teff.append(float(linelist[2].strip()))
            B_V.append(float(linelist[5].strip()))

## creates a vector of zeros, then I can fill it up
## with the values I want
Av = np.zeros(len(Spt))

### estimate the Av of each star
for i in range(len(Spt)):
    for j in range(len(Spt_KH)):
        if Spt[i] == Spt_KH[j]:
            Av[i] = 3.1*(B_V[i]-B_V_KH[j])

### In the case you have to interpolate the B_V_KH color
### on the table you need to:

"""each spt on the KH table have a number, for example
B0 is the number 1, B1 number 2, B2 number 3 and so on..
so you have a vector like this x=[1,2,3,...N].
in this way I can refer to a certain Spt by a number. Now if I
want to know the B_V color for a Spt of B1.5, I have to interpolate
the B_V_KH color between spts B1 and B2. In orden to do that
I use the numpy linear interpolation command like this:
np.interp(1.5,x,y), where x = [1,2,3,...N] and
y = [B_V_KH[1],B_V_KH[2],B_V_KH[3],...], e.g is your B_V_KH vector.
The 1.5 is just the x-coordinates of the interpolated values.
Like this I will have the B_V_KH color for a
Spt between B1 and B2 --> B1.5. The vector x must be increasing.
Note that you can use and array to define the x-coor of the
interpolated values. So you can obtain a vector of
interpolated values of B_V colors and then appended to you original vector
of B_V_KH"""



