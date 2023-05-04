import os
dirname, filename = os.path.split(os.path.abspath(__file__))
mytemplates = dirname + '/templates'
os.environ['PROTOCOL2_TEMPLATE']=mytemplates+'/protocol2.html'
os.environ['ANCIENT_TEMPLATE']=mytemplates+'/ancient.html'
os.environ['DHDS_TEMPLATE']=mytemplates+'/dhds.html'
os.environ['CLASSCOMPARE_TEMPLATE']=mytemplates+'/classcompare.html'
