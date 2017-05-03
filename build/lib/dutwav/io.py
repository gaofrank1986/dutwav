# @func: print node dictionary [ndic] to file at [filepath]
def print_nodedic(filepath,ndic):
    f = open(filepath,"w")
    tmp =[]
    for i in ndic:
        line = "{:10d}  {:12.6f} {:12.6f} {:12.6f}".format(i,ndic[i][0],ndic[i][1],ndic[i][2])
        #for j in range(ndic[i][1]):
            #line += " {:12.6f}".format(ndic[i][2][j])
        line += "\n"
        tmp.append(line)
    f.writelines(tmp)
    f.close()

# @func: print elem dict [ndict] at [filepath] with string prefix [line]
# @param: [line] is prefix put at beginning
def print_elemdic(filepath,ndic,line):
    f = open(filepath,"w")
    if (line[-1]!='\n'):
        line+='\n'
    tmp =[line]
    for i in ndic:
        line = "{:5d}  {:5d}".format(i,ndic[i][1])
        for j in range(ndic[i][1]):
            line += " {:5d}".format(ndic[i][2][j])
        line += "\n"
        tmp.append(line)
    f.writelines(tmp)
    f.close()

