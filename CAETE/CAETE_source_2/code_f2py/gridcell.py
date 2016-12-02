import numpy as np
 
class gridcell:
    
    def __init__(self, x, y, cell_id):
        
        self.x = np.int32(x)
        
        self.y = np.int32(y)
        
        self.cell_id = cell_id


        self.shp = (self.x, self.y)

        self.attrs = dict()

        
    def __str__(self):
        
        return "gridcell at x = %d; y=%d" %(self.x, self.y)
    
    def getPos(self):
        return self.shp
