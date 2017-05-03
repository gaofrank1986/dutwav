# @func: return a function define the damping zone characteristics using [r],[alpha],[L]
# @param: [L] define the length of damping zone
# @param: [r] define the beginning position 
# @param: [alpha] is coefficient
def def_damp_func(r,alpha,L):
  r = float(r)
  alpha = float(alpha)
  L = float(L)
  def f(x):
     if x<r:
        return 0.
     else:
        value = (alpha*((x-r)/L)**2)
        if value>1:
            return 1.
        return value
  return f

# @func: generate sippem kind of mesh for one element
# @var : [nid] is node number for source point
# @var : [eid] element id for target elem
def output_ss(path,mesh,eid,nid):
    with open(path,"wb") as f:
        f.write(' 3      8      1      8    3.   8\n')
        info = mesh.elems[eid][2]
        j=0
        for i in [1,3,5,7,2,4,6,8]:
            j+=1 
            str1='{0:<5d}'.format(j)
            str1+=('    '.join('{0:<14.8f}'.format(i) for i in mesh.nodes[info[i-1]]))
            if (info[i-1]==nid):
                jth=i
            str1+='\n'
            f.write(str1)
        f.write('1     1     2    3     4    5   6   7   8   -3\n')
        str2=('    '.join('{0:<14.8f}'.format(i) for i in mesh.nodes[nid]))
        str2+='\n'
        f.write(str2)
        if j==1: xi=[-1.,-1.]
        if j==2: xi=[1.,-1.]
        if j==3: xi=[1.,1.]
        if j==4: xi=[-1.,1.]
        if j==5: xi=[0.,-1.]
        if j==6: xi=[1.,0.]
        if j==7: xi=[0.,1.]
        if j==8: xi=[-1.,0.]
        f.write('{0:<14.8f}  {0:<14.8f}'.format(xi[0],xi[1]))
        
