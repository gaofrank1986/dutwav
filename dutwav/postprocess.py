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
def create_animation(mesh,prefix,postfix='.out',loc="./",kind=1):
    import dutwav.mesh
    from os.path import isfile,isdir
    from os import mkdir,rename
    from shutil import rmtree
    assert(isinstance(mesh,dutwav.mesh.Mesh))
    fd=loc+prefix+'_anim'
    fd2=loc+prefix
    if(not(isdir(fd))):
        mkdir(fd)
    else:
        rmtree(fd)
        mkdir(fd)
    if(not(isdir(fd2))): 
        mkdir(fd2)

    i=0
    while(True):
        if(kind==1):
            infile=loc+prefix+'.'+'{0:0>7d}'.format(i)+postfix
            mov_dest=loc+prefix+'/'+prefix+'.'+'{0:0>7d}'.format(i)+postfix
        else:
            infile=loc+prefix+'/'+prefix+'.'+'{0:0>7d}'.format(i)+postfix

        if(not(isfile(infile))):
            break
        
        print infile
        outfile=fd+'/'+prefix+'.'+'{0:0>7d}'.format(i)+'.dat'
        (r1,r2,r3)=parse_result(infile)
        mesh.tecplt_surface(outfile,[r1,r2,r3],i)
        # print mov_dest
        if(kind==1):
            rename(infile,mov_dest)
        i=i+1
