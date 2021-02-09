"""
Created on Fri Dec 15 21:54:40 2017

@author: abhijit
"""
class Reader:
    
    def __init__(self,f):
         self.data = self.__read(f)
        
    @staticmethod
    def __read(f):
        """Returns tuple of tuple of read data
        """
        return tuple(tuple([line[0:6].strip(),line[6:11].strip(),\
                            line[12:16].strip(),line[16:17].strip(),\
                            line[17:20].strip(),line[21:22].strip(),\
                            line[22:26].strip(),\
                            line[26:27].strip(),line[30:38].strip(),\
                            line[38:46].strip(),line[46:54].strip(),\
                            line[54:60].strip(),line[60:66].strip(),\
                            line[76:78].strip(),line[78:80].strip()])\
                            for line in open(f,'r') if ('ATOM' or 'HETATM') in line[0:6].strip())

    
class Writer:
    
    def __init__(self,*args):
        self.out=self.__write(args)
    
    @staticmethod
    
    def __write(args):
        pdb=dict(zip([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],args))
        return "{PDB[1]:6s}{PDB[2]:5d} {PDB[3]:^4s}{PDB[4]:1s}{PDB[5]:3s} {PDB[6]:1s}{PDB[7]:4d}{PDB[8]:1s}   {PDB[9]:8.3f}{PDB[10]:8.3f}{PDB[11]:8.3f}{PDB[12]:6.2f}{PDB[13]:6.2f}          {PDB[14]:>2s}{PDB[15]:2s}".format(PDB=pdb)
    
    def write(self):
        return self.out