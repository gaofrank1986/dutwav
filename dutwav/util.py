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


def get_node_linking(edic):
    node_elem_rela ={}
    for i in edic:
        for j in range(edic[i][1]):
            node = edic[i][2][j]
            if node in node_elem_rela:
                node_elem_rela[node]['elem'].append(i)
                node_elem_rela[node]['pos'].append(j+1)
            else:
                node_elem_rela[node] ={}
                node_elem_rela[node]['elem']=[i]
                node_elem_rela[node]['pos']=[j+1]
    for i in node_elem_rela:
        node_elem_rela[i]['num_links'] = len(node_elem_rela[i]['elem'])
        node_elem_rela[i]['elem'] = tuple(node_elem_rela[i]['elem'])
        node_elem_rela[i]['pos'] = tuple(node_elem_rela[i]['pos'])
    return node_elem_rela


def get_nrml_on_node(edic):
    nelem = len(edic.keys())
    node_nrml_rela={}
    print "{:d} elemnts".format(nelem)
    for i in range(nelem):
        etype = edic[i+1][1]
        for k in range(etype):#not all elem has 8 node
            nid = edic[i+1][2][k]
            mid = edic[i+1][3][k]
            if nid in node_nrml_rela:
                #print nid
                node_nrml_rela[nid].add(mid)
            else:
                node_nrml_rela[nid]=set()
                node_nrml_rela[nid].add(mid)
    return node_nrml_rela

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

def set_aspect_equal_3d(ax):
    """Fix equal aspect bug for 3D plots."""

    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)

    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])

    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])


def get_sufrace_potential(m,num,amp,omeg,depth,wk):
    import numpy as np
    g=9.801
    res={}
    for i in range(num):
        x=m.nodes[i+1][0]
        y=m.nodes[i+1][1]
        z=m.nodes[i+1][2]
        theta = np.arctan2(y,x)
        r=np.sqrt(x**2+y**2)
        phi=-1j*g*amp/omeg
        phi=phi*np.cosh(wk*(z+depth))/np.cosh(wk*depth)
        phi=phi*np.exp(1j*wk*r*np.cos(theta))
        res[i+1]=phi
    return res

def get_ind_potential(mesh,num,amp,omeg,depth,wk): 
    import numpy as np
    from dutwav.analytical import bj
    from scipy.special import jv
    mesh.nodes[1]
    g=9.801
    res={}
    for i in range(num):
        x=mesh.nodes[i+1][0]
        y=mesh.nodes[i+1][1]
        z=mesh.nodes[i+1][2]
        theta = np.arctan2(y,x)
        r=np.sqrt(x**2+y**2)
        tmp=0.
        for m in range(10):
            eps=2.
            if m==0:
                eps=1.
            tp=-1j*g*amp/omeg
            tp=tp*np.cosh(wk*(z+depth))/np.cosh(wk*depth)
            tp=tp*eps*(1j)**m*bj(wk*r,m)*np.cos(m*theta)
            # tp=tp*eps*(1j)**m*jv(m,wk*r)*np.cos(m*theta)
            tmp+=tp
        res[i+1]=tmp
    return res
