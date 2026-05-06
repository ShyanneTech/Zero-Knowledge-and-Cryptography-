# 41 X103 = 4223 to R1CS

import numpy as np

O= np.matrix([[0,1,0,0]])
L= np.matrix([[0,0,1,0]])
R= np.matrix([[0,0,0,1]])

# Witness vector
a = np.array([1,4223,41,103])


# element wise for each constraint 
result = np.matmul(O,a) == np.matmul(L,a) * np.matmul(R,a) 

print(result.all(), "result contains an inequality")
assert result.all(), "result contains an inequality"