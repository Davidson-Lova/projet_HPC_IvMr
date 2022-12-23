import numpy as np

res = np.loadtxt('res.txt',dtype="float64")

bla = int(np.sqrt(res.shape[0]))
bla_2 = bla//2
bla_l = res[1:1,:]
for i in range(bla) :
    bla_t = np.concatenate((res[i*bla_2:(i+1)*bla_2,:], res[i*bla_2 + bla_2*bla: (i+1)*bla_2 + bla_2*bla,:]))
    bla_l = np.concatenate((bla_l, bla_t))

content = ""
count = 0
for i in range(bla_l.shape[0]) :
    content += "  "
    for j in range(bla_l.shape[1]) :
        content += "{} ".format(bla_l[i,j])
    content += "\n"
    if (count == bla - 1) :
        content += "  \n"
    if(count == bla) :
        count = 0
    count += 1

file = open("res.txt", "w+")
file.write(content)
file.close()