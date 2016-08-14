# @func: generate r1,r2,r3
# @param: [path] = file path
# @desrcipt: read r1,r2 from each line of file, calc r3, return (r1,r2,r3)
def parse_result(path):
    with open(path,"rb") as f:
        r = f.readlines()
    r1={}
    r2={}
    r3={}
    for i in range(len(r)):
        tmp=[float(j) for j in r[i].split()]
        r1[i+1]=tmp[1]
        r2[i+1]=tmp[2]
        r3[i+1]=tmp[2]-tmp[1]
    return (r1,r2,r3) 

def display_value(fi,fo,mesh):
    import dutwav.mesh
    assert(isinstance(mesh,dutwav.mesh.Mesh))
    with open(fi,"rb") as f:
        r = f.readlines()
        r1={}
        for i in range(len(r)):
            tmp=[float(j) for j in r[i].split()]
            #tmp[0] is node id
            r1[i+1]=tmp[1]
    # print len(r1),len(mesh.nodes)
    mesh.tecplt_surface(fo,[r1],1,kind=2)  


# @param path : file path
# @param mesh : mesh object
# @param ubound : upper bounder of time counting from time 000
# @param prefix : 
def create_animation(path,mesh,ubound,prefix="./fort.7"):
    import dutwav.mesh
    assert(isinstance(mesh,dutwav.mesh.Mesh))
    for i in range(ubound):
        infilename=prefix+'{0:0>3d}'.format(i)
        outfilename="./tecplt_animation_"+'{0:0>4d}'.format(i)+'.dat'
        (r1,r2,r3)=parse_result(infilename)
        mesh.tecplt_surface(outfilename,[r1,r2,r3],i)

