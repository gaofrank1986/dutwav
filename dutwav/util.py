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

