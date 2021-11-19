import matplotlib 
print(matplotlib.get_backend())

import matplotlib.pyplot as plt 

x = range(1,10)

plt.plot(x,x)
plt.show()