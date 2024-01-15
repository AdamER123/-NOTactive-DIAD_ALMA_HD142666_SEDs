#program for getting HAEBE spectral type from SIMBAD

import urllib, urllib2

i = 0


#resetting file of spectral types for the HAEBE stars
with open('HAEBE_SpT.txt', 'w') as f:
	f.write('Star' + '	 ' + 'SpT')

table = open('HAEBEinput.txt', 'r') #for names from catalogs of HAEBE stars, edited to include as many as possible
lines = table.readlines()

#beginning loop that reads through name inputs from HAEBEinput.txt
for obj in lines:

	url = 'http://simbad.u-strasbg.fr/simbad/sim-id?Ident=' #URL that is queried to get data
	obj = obj.replace(" ", "+")

	#encoding search to get data from SIMBAD, put into variable raw
	p = {}
	p['Ident'] = obj
	'''p['NbIdent'] = 1
	p['Radius'] = 2    
	p['Radius.unit'] = 'arcmin'
	p['submit'] = 'submit+id'''

	#p = '&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id'
	query = urllib.urlencode(p)
	get_url = url + query
	#get_url = url + obj + query
	handler = urllib2.urlopen(get_url)
	raw = handler.readlines()

	print raw[135]

	#appending spectral type data for each HAEBE star
	with open('HAEBE_SpT.txt', 'a') as f:
		#try:
		f.write(obj + '		' + raw[430]) #appends line with data, second to last line, from the SIMBAD's information about each star (used HTML source to figure out line)
		#except IndexError:
			#print obj

	i +=1
	print i

table.close()
