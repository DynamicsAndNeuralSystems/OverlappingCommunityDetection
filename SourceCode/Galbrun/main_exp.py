#!/usr/bin/python
from experimental import *
import pdb


########################################################
################# MINING ACTIONS
### Denset subgraph then choose labels, possibly stealing from previous cliques (cf. Cohen and Katzir)
def findDensestCK(graph, nodes_to_labels, data_parameters):        
    ticE = datetime.datetime.now()
    min_supp = data_parameters["min_supp"] #*(data_parameters["min_supp"]-1)
    ks = data_parameters["ks"]
    external = "nlgd"
    print "------------METHOD=%s K=%s" % (external, ks)
    
    current_graph = dict([(n, set(v.keys())) for n,v in graph.items()])

    labels_to_nodes = getReversedLabels(nodes_to_labels, graph, data_parameters)
    unassigned = len(graph)+10
    edges_nbassi = {}
    assE = {}
    for n in graph.keys():
        for nei in graph[n].keys():
            assE[(nei,n)] = (unassigned, None)
            if nei < n:
                edges_nbassi[(nei,n)] = 1
    nb_n = float(len(graph))
    nb_e = float(len(edges_nbassi))
    all_edges = set(edges_nbassi.keys())

    stop = False    
    ti = -1
    while not stop:
        history = []
        n = None
        lc = float(len(current_graph))
        if lc > min_supp:
            dn = dict([(x, sum([assE[(nei,x)][0] > lc for nei in vs])) for x, vs in current_graph.items()])
            ss = sum(dn.values())
        else:
            ss = 0
        while ss > 0:
            history.append((lc/ss, n))
            n = sorted(dn.keys(), key=lambda x: (dn[x], x))[0]
            for neigh in current_graph[n]:
                current_graph[neigh].remove(n)
            del current_graph[n]
            lc = float(len(current_graph))
            if lc > min_supp:
                dn = dict([(x, sum([assE[(nei,x)][0] > lc for nei in vs])) for x, vs in current_graph.items()])
                ss = sum(dn.values())
            else:
                ss = 0
            
        history.append((float("Inf"), n))
        for n in current_graph.keys():
            history.append((float("Inf"), n))
        t = min([(v[0],-i) for i,v in enumerate(history)])

        if t[0] == float("Inf"):
            stop = True
        else:
            ti += 1
            # print "[%s]" % ", ".join([h[1] for h in history[1:-t[1]+1]])
            dense = set([h[1] for h in history[-t[1]+1:]])
            ld = len(dense)
            edges_dense = set()
            for node in dense:
                for nei in dense.intersection(graph[node].keys()):
                    if nei < node and assE[(nei, node)]  > ld:
                        assE[(nei,node)] = (ld, ti)
                        assE[(node,nei)] = (ld, ti)
                        edges_dense.add((nei, node))  
            current_graph = dict([(n, set(v.keys())) for n,v in graph.items()])

        if ti+1 >= max(ks):
            stop = True

        if ti+1 in ks or stop:
            edgesets = [set() for i in range(ti+1)]
            for edge, v in assE.items():
                if edge[0] > edge[1] and v[1] != None:
                    edgesets[v[1]].add(edge)

            nodesets = []
            for edgeset in edgesets:
                nodesets.append(getNodeset(edgeset))

            idsort = sorted(range(len(nodesets)), key=lambda x: (len(nodesets[x])))

            covered = {}
            coveredn = {}

            fo, fe = initResOut(data_parameters, external, ti+1)
            
            for ti, tt in enumerate(idsort):
                labelset, lspec = chooseLabels(labels_to_nodes, nodesets[tt], data_parameters)
                tmp, ne, nn = addGroup(ti, score(len(edgesets[tt]), len(nodesets[tt])), lspec, edgesets[tt], nodesets[tt], labelset=labelset, covered=covered, coveredn=coveredn, nb_e=nb_e, nb_n=nb_n, all_edges=all_edges, graph=graph)
                fo.write(tmp+"\n")
                if fe is not None:
                    # fe.write(" ".join(sorted([getNodeId(n, graph, data_parameters) for n in nodesets[tt]])) + "\n")
                    fe.write("%d\t%d\t%s\t%s\n" % ( len(nodesets[tt]), len(labelset),
                                                   " ".join(sorted(nodesets[tt])),
                                                   " ".join(sorted(labelset))))

            fo.write("Done in %s\n" % (datetime.datetime.now() - ticE))
            if fe is not None:
                fo.close()
                fe.close()

### Denset subgraph then choose labels, no stealing from previous cliques 
def findDensest(graph, nodes_to_labels, data_parameters):
    min_supp = data_parameters["min_supp"] #*(data_parameters["min_supp"]-1)
    ks = data_parameters["ks"]
    external = "dense"
    print "------------METHOD=%s K=%d" % (external, ks)
    
    kmax = 0
    if 0 not in ks:
        kmax = max(ks)

    fo, fe = initResOut(data_parameters, external, kmax)
    ticE = datetime.datetime.now()

    current_graph = dict([(n, set(v.keys())) for n,v in graph.items()])
    uncovered_graph = dict([(n, set(v.keys())) for n,v in graph.items()])

    covered = {}
    coveredn = {}
    edges_nbassi = set()
    for n in graph.keys():
        edges_nbassi.update([tuple(sorted([nei, n])) for nei in graph[n].keys()])
    edges_nbassi = dict([(e, 1) for e in edges_nbassi])
    nb_n = float(len(graph))
    nb_e = float(len(edges_nbassi))
    all_edges = set(edges_nbassi.keys())

    labels_to_nodes = getReversedLabels(nodes_to_labels, graph, data_parameters)
    stop = False    
    ti = -1
    while not stop:
        ti += 1
        history = []
        n = None
        while len(current_graph) > min_supp and sum([len(v) for v in current_graph.values()]) > 0:
            history.append((len(current_graph)*2.0/sum([len(v) for v in current_graph.values()]), n))
            n = sorted(current_graph.keys(), key=lambda x:(len(current_graph[x]), x))[0]
            for neigh in current_graph[n]:
                current_graph[neigh].remove(n)
            del current_graph[n]
            
        history.append((float("Inf"), n))
        for n in current_graph.keys():
            history.append((float("Inf"), n))
        t = min([(v[0],-i) for i,v in enumerate(history)])

        if t[0] == float("Inf"):
            stop = True
        else:
            # print "[%s]" % ", ".join([h[1] for h in history[1:-t[1]+1]])
            dense = set([h[1] for h in history[-t[1]+1:]])

            edges_dense = set()
            for node in dense:
                edges_dense.update([tuple(sorted([nei, node])) for nei in dense.intersection(graph[node].keys())])
                uncovered_graph[node] -= dense

            labelset, lspec = chooseLabels(labels_to_nodes, dense, data_parameters)
            tmp, ne, nn = addGroup(ti, t[0], lspec, edges_dense, nodeset=dense, labelset=labelset, covered=covered, coveredn=coveredn, nb_e=nb_e, nb_n=nb_n, all_edges=all_edges, graph=graph)
            fo.write("Here in %s\n" % (datetime.datetime.now() - ticE))
            fo.write(tmp+"\n")

            if fe is not None:
                #fe.write(" ".join(sorted([getNodeId(n, graph, data_parameters) for n in dense])) + "\n")
                fe.write("%d\t%d\t%s\t%s\n" % ( len(dense), len(labelset),
                                               " ".join(sorted(dense)),
                                               " ".join(sorted(labelset))))


            current_graph = dict([(n, set(v)) for n,v in uncovered_graph.items()])

        if 0 not in ks and ti+1 >= max(ks):
            stop = True
            
    fo.write("Done in %s\n" % (datetime.datetime.now() - ticE))
    if fe is not None:
        fo.close()
        fe.close()

### Denset subgraph using labels, with/without spreading labels at mining/evaluation time 
### this allows for stealing edges from previously found cliques (Cohen Katzir)..??
def mineGreedyCK(graph, nodes_to_labels, data_parameters):
    spread_mine=data_parameters["spread_mine"]
    spread_eval=data_parameters["spread_eval"]
    min_supp = data_parameters["min_supp"] #*(data_parameters["min_supp"]-1)
    ks = data_parameters["ks"]
    external = "lgdck%d%d" % (spread_mine, spread_eval)
    
    kmax = 0
    if 0 not in ks:
        kmax = max(ks)
    
    nodes = set(graph.keys())
    if spread_mine or spread_eval:
        ## TIME ticT = datetime.datetime.now()
        tmp_dict = spreadLabels(graph, data_parameters["ndlb_dict"], data_parameters)
        nodes_to_labels_spread = getNodesToLabels(graph, tmp_dict, data_parameters)
        labels_to_nodes_spread = getReversedLabels(nodes_to_labels_spread, graph, data_parameters, store=False)
        labels_to_edges_spread, edges_to_labels_spread = getEdgesCands(nodes_to_labels_spread, graph, data_parameters, store=False)
        ## TIME print "Spreading T=%s" % (datetime.datetime.now() - ticT)
        
    if not spread_mine or not spread_eval:
        ## TIME ticT = datetime.datetime.now()
        nodes_to_labels_unspread = nodes_to_labels
        labels_to_nodes_unspread = getReversedLabels(nodes_to_labels_unspread, graph, data_parameters)
        labels_to_edges_unspread, edges_to_labels_unspread = getEdgesCands(nodes_to_labels_unspread, graph, data_parameters)
        ## TIME print "Reversing T=%s" % (datetime.datetime.now() - ticT)

    ## TIME nn = float(len(nodes))
    ## TIME print "NODES TO LABELS: av=%f ap=%f" % ( numpy.sum([len(v) for v in nodes_to_labels_unspread.values()])/nn, numpy.sum([len(v) for v in nodes_to_labels_spread.values()])/nn)
    ## TIME print "LABELS TO NODES: av=%f ap=%f" % ( numpy.mean([len(v) for v in labels_to_nodes_unspread.values()]), numpy.mean([len(v) for v in labels_to_nodes_spread.values()]))

    
    if spread_mine:
        labels_to_nodes, nodes_to_labels = (labels_to_nodes_spread, nodes_to_labels_spread)
        labels_to_edges, edges_to_labels = (labels_to_edges_spread, edges_to_labels_spread)
    else:
        labels_to_nodes, nodes_to_labels = (labels_to_nodes_unspread, nodes_to_labels_unspread)
        labels_to_edges, edges_to_labels = (labels_to_edges_unspread, edges_to_labels_unspread)

    if spread_eval:
        labels_to_nodes_eval, nodes_to_labels_eval = (labels_to_nodes_spread, nodes_to_labels_spread)
        labels_to_edges_eval, edges_to_labels_eval = (labels_to_edges_spread, edges_to_labels_spread)
    else:
        labels_to_nodes_eval, nodes_to_labels_eval = (labels_to_nodes_unspread, nodes_to_labels_unspread)
        labels_to_edges_eval, edges_to_labels_eval = (labels_to_edges_unspread, edges_to_labels_unspread)

    if not spread_eval and spread_mine:
        type_support = "maj"
    else:
        type_support = "int"
    
    ticE = datetime.datetime.now()

    ## TIME ticT = datetime.datetime.now()
    candidates = sorted([l for l,v in labels_to_nodes.items() if len(v) > min_supp and labels_to_edges.has_key(l)])
    ## TIME print "CANDIDATES", len(candidates)
    ## TIME print "Sorting T=%s" % (datetime.datetime.now() - ticT)
        
    ### reduce labels on edges to frequent ones
    # for k in labels.keys():
    #     labels[k] &= candidates
    nb_e = float(len(edges_to_labels))
    nb_n = float(len(nodes))
    covered = {} #set()
    coveredn = {} #set()
    all_edges = set(edges_to_labels.keys())
    sizes = {None:len(nodes)+1}
    selected = []

    top = (1, 1, 1, 1)
    ti = -1
    stop = False
    while top[1] is not None and not stop:

        ## TIME ticT = datetime.datetime.now()
        ti += 1
        keep = []
        top = (1, None, None, None)
        new_top = (float("Inf"), None, None, None)
        for k in candidates:
            ucov = labels_to_edges[k]
            esupp = [e for e in labels_to_edges[k] if sizes[covered.get(e, None)] > len(labels_to_nodes[k])]
            r = score(len(esupp), len(labels_to_nodes[k]))
            if r < new_top[0]:
                new_top = (r, [k], set(labels_to_nodes[k]), ucov)

        while new_top[1] is not None:
            top = new_top
            keep.append(new_top[0])
            new_top = (float("Inf"), None, None, None)
            todoe = sorted([k for k,v in labels_to_nodes.items() if k not in top[1] and len(v & top[2]) > min_supp and labels_to_edges.has_key(k)])
            for l in todoe:
                if True: #len(labels_to_edges[l].intersection(top[3])) > min_supp:
                
                    esupp = [e for e in labels_to_edges[l].intersection(top[3]) if sizes[covered.get(e, None)] > len(labels_to_nodes[l].intersection(top[2]))]
                    r = score(len(esupp), len(labels_to_nodes[l].intersection(top[2])))
                    if r < new_top[0]: # and r != 1:
                        new_top = (r, [l]+top[1], labels_to_nodes[l].intersection(top[2]), labels_to_edges[l].intersection(top[3]))
        ## TIME print "Candidate T=%s" % (datetime.datetime.now() - ticT)

        if top[1] is not None:
            tscore, pos = min([(v,-i) for (i,v) in enumerate(keep)])
            labelset = top[1][pos-1:]
            ## TIME ticT = datetime.datetime.now()
            cnodes, cedges = getSupport(type_support, labelset, labels_to_nodes_eval, labels_to_edges_eval, nodes_to_labels=nodes_to_labels_eval, edges_to_labels=edges_to_labels_eval)
            ## TIME print "Support T=%s" % (datetime.datetime.now() - ticT)
            if len(cnodes) > 1 and len(cedges) > 0:
                for e in cedges:
                    if sizes[covered.get(e, None)] > len(cnodes):
                        covered[e] = ti
                    sizes[ti] = len(cnodes)
                selected.append({"nodeset": cnodes, "labelset": labelset})
            else:
                ti -= 1
            if spread_mine and not spread_eval:
                ## TIME ticT = datetime.datetime.now()
                updateSpreads(graph, data_parameters["ndlb_dict"], labelset, cnodes, labels_to_nodes, labels_to_edges, labels_to_nodes_eval, data_parameters)
                ## TIME print "Update T=%s" % (datetime.datetime.now() - ticT)
            if 0 not in ks and ti+1 >= max(ks):
                stop = True

            if ti+1 in ks or stop:
                for ci in range(len(selected)):
                    selected[ci]["edgeset"] = set()
                for e, ci in covered.items():
                    selected[ci]["edgeset"].add(e) 
                printSelected(selected, graph, data_parameters, external, ticE, nb_n, nb_e, all_edges)

    if len(selected) > 0  and ti+1 < min(ks):
        for ci in range(len(selected)):
            selected[ci]["edgeset"] = set()
        for e, ci in covered.items():
            selected[ci]["edgeset"].add(e) 
        printSelected(selected, graph, data_parameters, external, ticE, nb_n, nb_e, all_edges)


### Denset subgraph using labels, with/without spreading labels at mining/evaluation time 
def mineGreedy(graph, nodes_to_labels, data_parameters):
    spread_mine=data_parameters["spread_mine"]
    spread_eval=data_parameters["spread_eval"]
    min_supp = data_parameters["min_supp"] #*(data_parameters["min_supp"]-1)
    ks = data_parameters["ks"]
    external = "lgd%d%d" % (spread_mine, spread_eval)
    
    kmax = 0
    if 0 not in ks:
        kmax = max(ks)
    fo, fe = initResOut(data_parameters, external, kmax)
    
    nodes = set(graph.keys())
    if spread_mine or spread_eval:
        ## TIME ticT = datetime.datetime.now()
        tmp_dict = spreadLabels(graph, data_parameters["ndlb_dict"], data_parameters)
        nodes_to_labels_spread = getNodesToLabels(graph, tmp_dict, data_parameters)
        labels_to_nodes_spread = getReversedLabels(nodes_to_labels_spread, graph, data_parameters, store=False)
        labels_to_edges_spread, edges_to_labels_spread = getEdgesCands(nodes_to_labels_spread, graph, data_parameters, store=False)
        ## TIME print "Spreading T=%s" % (datetime.datetime.now() - ticT)
        
    if not spread_mine or not spread_eval:
        ## TIME ticT = datetime.datetime.now()
        nodes_to_labels_unspread = nodes_to_labels
        labels_to_nodes_unspread = getReversedLabels(nodes_to_labels_unspread, graph, data_parameters)
        labels_to_edges_unspread, edges_to_labels_unspread = getEdgesCands(nodes_to_labels_unspread, graph, data_parameters)
        ## TIME print "Reversing T=%s" % (datetime.datetime.now() - ticT)

    ## TIME nn = float(len(nodes))
    ## TIME print "NODES TO LABELS: av=%f ap=%f" % ( numpy.sum([len(v) for v in nodes_to_labels_unspread.values()])/nn, numpy.sum([len(v) for v in nodes_to_labels_spread.values()])/nn)
    ## TIME print "LABELS TO NODES: av=%f ap=%f" % ( numpy.mean([len(v) for v in labels_to_nodes_unspread.values()]), numpy.mean([len(v) for v in labels_to_nodes_spread.values()]))

    
    if spread_mine:
        labels_to_nodes, nodes_to_labels = (labels_to_nodes_spread, nodes_to_labels_spread)
        labels_to_edges, edges_to_labels = (labels_to_edges_spread, edges_to_labels_spread)
    else:
        labels_to_nodes, nodes_to_labels = (labels_to_nodes_unspread, nodes_to_labels_unspread)
        labels_to_edges, edges_to_labels = (labels_to_edges_unspread, edges_to_labels_unspread)

    if spread_eval:
        labels_to_nodes_eval, nodes_to_labels_eval = (labels_to_nodes_spread, nodes_to_labels_spread)
        labels_to_edges_eval, edges_to_labels_eval = (labels_to_edges_spread, edges_to_labels_spread)
    else:
        labels_to_nodes_eval, nodes_to_labels_eval = (labels_to_nodes_unspread, nodes_to_labels_unspread)
        labels_to_edges_eval, edges_to_labels_eval = (labels_to_edges_unspread, edges_to_labels_unspread)

    if not spread_eval and spread_mine:
        type_support = "maj"
    else:
        type_support = "int"
    
    ticE = datetime.datetime.now()

    ## TIME ticT = datetime.datetime.now()
    candidates = sorted([l for l,v in labels_to_nodes.items() if len(v) > min_supp and labels_to_edges.has_key(l)])
    ## TIME print "CANDIDATES", len(candidates)
    ## TIME print "Sorting T=%s" % (datetime.datetime.now() - ticT)
        
    ### reduce labels on edges to frequent ones
    # for k in labels.keys():
    #     labels[k] &= candidates
    nb_e = float(len(edges_to_labels))
    nb_n = float(len(nodes))
    covered = {} #set()
    coveredn = {} #set()
    all_edges = set(edges_to_labels.keys())

    top = (1, 1, 1, 1)
    ti = -1
    stop = False
    while top[1] is not None and not stop:

        ## TIME ticT = datetime.datetime.now()
        ti += 1
        if type(covered) == dict:
            tco = set(covered.keys())
        else:
            tco = covered
        keep = []
        top = (1, None, None, None)
        new_top = (float("Inf"), None, None, None)
        for k in candidates:
            ucov = labels_to_edges[k] - tco
            r = score(len(ucov), len(labels_to_nodes[k]))
            if r < new_top[0]:
                new_top = (r, [k], set(labels_to_nodes[k]), ucov)

        while new_top[1] is not None:
            top = new_top
            keep.append(new_top[0])
            new_top = (float("Inf"), None, None, None)
            todoe = sorted([k for k,v in labels_to_nodes.items() if k not in top[1] and len(v & top[2]) > min_supp and labels_to_edges.has_key(k)])
            for l in todoe:
                if True: #len(labels_to_edges[l].intersection(top[3])) > min_supp:
                    r = score(len(labels_to_edges[l].intersection(top[3])), len(labels_to_nodes[l].intersection(top[2])))
                    if r < new_top[0]: # and r != 1:
                        new_top = (r, [l]+top[1], labels_to_nodes[l].intersection(top[2]), labels_to_edges[l].intersection(top[3]))
        ## TIME print "Candidate T=%s" % (datetime.datetime.now() - ticT)

        if top[1] is not None:
            tscore, pos = min([(v,-i) for (i,v) in enumerate(keep)])
            labelset = top[1][pos-1:]
            ## TIME ticT = datetime.datetime.now()
            cnodes, cedges = getSupport(type_support, labelset, labels_to_nodes_eval, labels_to_edges_eval, nodes_to_labels=nodes_to_labels_eval, edges_to_labels=edges_to_labels_eval)
            ## TIME print "Support T=%s" % (datetime.datetime.now() - ticT)
            if type(covered) == dict:
                tco = set(covered.keys())
            else:
                tco = covered
            nr = score(len(cedges-tco), len(cnodes))

            tmp, ne, nn = addGroup(ti, nr, 1, cedges, cnodes, labelset=labelset, covered=covered, coveredn=coveredn, nb_e=nb_e, nb_n=nb_n, all_edges=all_edges, graph=graph)

            if spread_mine and not spread_eval:
                ## TIME ticT = datetime.datetime.now()
                updateSpreads(graph, data_parameters["ndlb_dict"], labelset, cnodes, labels_to_nodes, labels_to_edges, labels_to_nodes_eval, data_parameters)
                ## TIME print "Update T=%s" % (datetime.datetime.now() - ticT)

            fo.write("Here in %s\n" % (datetime.datetime.now() - ticE))
            fo.write(tmp+"\n")
            
            if fe is not None:
                # fe.write(" ".join(sorted([getNodeId(n, graph, data_parameters) for n in cnodes])) + "\n")
                fe.write("%d\t%d\t%s\t%s\n" % ( len(cnodes), len(labelset),
                                                " ".join(sorted(cnodes)),
                                                " ".join(sorted(labelset))))

            if 0 not in ks and ti+1 >= max(ks):
                stop = True
    fo.write("Done in %s\n" % (datetime.datetime.now() - ticE))
    if fe is not None:
        fo.close()
        fe.close()
        

def mineGreedyMaj(graph, nodes_to_labels, data_parameters):
    min_supp = data_parameters["min_supp"] #*(data_parameters["min_supp"]-1)
    ks = data_parameters["ks"]
    kmax = 0
    if 0 not in ks:
        kmax = max(ks)
    external = "maj"

    fo, fe = initResOut(data_parameters, external, kmax)

    nodes = set(graph.keys())    
    labels_to_nodes = getReversedLabels(nodes_to_labels, graph, data_parameters)
    labels_to_edges, edges_to_labels = getEdgesCands(nodes_to_labels, graph, data_parameters)

    ticE = datetime.datetime.now()
    nb_e = float(len(edges_to_labels))
    nb_n = float(len(nodes))
    covered = {} #set()
    coveredn = {} #set()
    all_edges = set(edges_to_labels.keys())
    uncovered_edges = set(all_edges)

    ti = -1
    while True:
        ti += 1
        labelset_tmp, scores = densestMaj(nodes_to_labels, labels_to_nodes, uncovered_edges)
        if len(labelset_tmp) == 0:
            break

        tscore, pos = min([(v,-i) for (i,v) in enumerate(scores)])
        labelset = labelset_tmp[:-pos+1]
        cnodes, cedges = getSupport("maj", labelset, labels_to_nodes, labels_to_edges, nodes_to_labels=nodes_to_labels, edges_to_labels=edges_to_labels)
        tmp, ne, nn = addGroup(ti, tscore, 1, cedges, cnodes, labelset=labelset, covered=covered, coveredn=coveredn, nb_e=nb_e, nb_n=nb_n, all_edges=all_edges, graph=graph)
        uncovered_edges.difference_update(covered) 
        fo.write("Here in %s\n" % (datetime.datetime.now() - ticE))
        fo.write(tmp+"\n")

        if fe is not None:
            # fe.write(" ".join(sorted([getNodeId(n, graph, data_parameters) for n in cnodes])) + "\n")
            fe.write("%d\t%d\t%s\t%s\n" % ( len(cnodes), len(labelset),
                                               " ".join(sorted(cnodes)),
                                               " ".join(sorted(labelset))))


        if 0 not in ks and ti+1 >= max(ks):
            break
    fo.write("Done in %s\n" % (datetime.datetime.now() - ticE))
    if fe is not None:
        fo.close()
        fe.close()

def mineGreedyMajCK(graph, nodes_to_labels, data_parameters):
    min_supp = data_parameters["min_supp"] #*(data_parameters["min_supp"]-1)
    ks = data_parameters["ks"]
    kmax = 0
    if 0 not in ks:
        kmax = max(ks)
    external = "majck"

    nodes = set(graph.keys())    
    labels_to_nodes = getReversedLabels(nodes_to_labels, graph, data_parameters)
    labels_to_edges, edges_to_labels = getEdgesCands(nodes_to_labels, graph, data_parameters)

    ticE = datetime.datetime.now()
    nb_e = float(len(edges_to_labels))
    nb_n = float(len(nodes))
    covered = {} #set()
    coveredn = {} #set()
    all_edges = set(edges_to_labels.keys())
    sizes = {None:nb_n+1}
    selected = []
    
    stop = False
    ti = -1
    while not stop:
        ti += 1
        labelset_tmp, scores = densestMaj(nodes_to_labels, labels_to_nodes, all_edges, covered, sizes) #,
                                          # {"le":labels_to_edges, "el": edges_to_labels,
                                          #  "ln":labels_to_nodes, "nl": nodes_to_labels})
        if len(labelset_tmp) == 0:
            stop = True

        else:
            tscore, pos = min([(v,-i) for (i,v) in enumerate(scores)])
            labelset = labelset_tmp[:-pos+1]
            cnodes, cedges = getSupport("maj", labelset, labels_to_nodes, labels_to_edges, nodes_to_labels=nodes_to_labels, edges_to_labels=edges_to_labels)
            ce = 0
            for e in cedges:
                if sizes[covered.get(e, None)] > len(cnodes):
                    covered[e] = ti
                    ce +=1
            sizes[ti] = len(cnodes)
            selected.append({"nodeset": cnodes, "labelset": labelset})

        if 0 not in ks and ti+1 >= max(ks):
            stop = True

        if ti+1 in ks or stop:
            for ci in range(len(selected)):
                selected[ci]["edgeset"] = set()
            for e, ci in covered.items():
                selected[ci]["edgeset"].add(e) 
            printSelected(selected, graph, data_parameters, external, ticE, nb_n, nb_e, all_edges)

def mineTriangles(graph, nodes_to_labels, data_parameters):
    min_supp = data_parameters["min_supp"]
    external = "tri"
    ticE = datetime.datetime.now()
    labels_to_nodes = getReversedLabels(nodes_to_labels, graph, data_parameters)
    labels_to_edges, edges_to_labels = getEdgesCands(nodes_to_labels, graph, data_parameters)
    candidates = []
    seen_labels = set()
    for node in graph.keys():
        tmp = mineTriangleNode(node, graph, nodes_to_labels, edges_to_labels, labels_to_nodes, labels_to_edges, data_parameters, seen_labels)
        if tmp is not None and len(tmp["nodeset"]) > min_supp:
            candidates.append(tmp)
    selectCandidates(candidates, graph, data_parameters, external, ticE)


def mineLabelsOrdering(graph, nodes_to_labels, data_parameters):
    min_supp = data_parameters["min_supp"]
    external = "lo"
    ticE = datetime.datetime.now()
    labels_to_nodes = getReversedLabels(nodes_to_labels, graph, data_parameters)
    labels_to_edges, edges_to_labels = getEdgesCands(nodes_to_labels, graph, data_parameters)
    lo, ldists = getLabelsOrdering(graph, labels_to_edges, data_parameters)
    candidates = getCandidatesFromLO(graph, labels_to_nodes, labels_to_edges, lo, data_parameters,
                                     nodes_to_labels=nodes_to_labels, edges_to_labels=edges_to_labels)
    candidates = [c for c in candidates if len(c["nodeset"]) > min_supp]
    selectCandidates(candidates, graph, data_parameters, external, ticE)
    

############################################################
###             ACTIONS
############################################################

ACTIONS = {"ego" : extractSubs,
           "nlgd": findDensestCK,
           "dense": findDensest,
           "lgd": mineGreedy,
           "lgdck": mineGreedyCK,
           "maj": mineGreedyMaj,
           "majck": mineGreedyMajCK,
           "tri": mineTriangles,
           "lo": mineLabelsOrdering}


################ MAIN
def run(args):
    org_parameters = None
    if len(args) > 1:
        action = args[1]
        if len(args) > 2 and os.path.isfile(args[2]):
            fp = open(args[2])
            org_parameters = readParameters(fp)
            fp.close()
    if org_parameters is None:
        Warning("No parameters given")
        exit()

    data_parameters = dict(org_parameters)

    if not action in ACTIONS:
        print "Action %s does not exist..." % action
        exit()
        
    elif action in ["table", "tables", "project"]:
        ACTIONS[action](data_parameters)
        exit()

    series = [{}]

    bfs = []
    print "PARAMETERS:\n" + ("\n".join(["  * %s= %s" % (k,v) for (k,v) in data_parameters.items()]))
    #### Load data network from file 
    if action != "ego" and data_parameters.has_key("centers"):
        for centert in data_parameters["centers"].strip(" ;").split(";"):

            center = centert.strip()
            bfs.append(("wegonet_%s%d-%.1f.graph" % (center.replace(" ",""), data_parameters["ego"], data_parameters["min_weight"]),
                        "wegonet_%s%d-%.1f.labels" % (center.replace(" ",""), data_parameters["ego"], data_parameters["min_weight"]), center))

    elif data_parameters.has_key('graph_file') and data_parameters.has_key('labels_file'):
        if "*" in data_parameters['graph_file']:
            starind = data_parameters['graph_file'].index("*")
            fnames = glob.glob(data_parameters['data_rep']+data_parameters['graph_file'])
            fnames.sort(key=lambda x: os.path.getsize(x))
            for gname in fnames:
                ser_b = gname[len(data_parameters['data_rep'])+starind:-len(data_parameters['graph_file'])+starind+1]
                lname = data_parameters['labels_file'].replace("*", ser_b)
                gname = gname[len(data_parameters['data_rep'])::]
                bfs.append((gname, lname, None))
        else:
            bfs.append((data_parameters['graph_file'], data_parameters['labels_file'], None))

    for bf in bfs:
        data_parameters['basis'] = ".".join(bf[0].split(".")[:-1])
        print "--------- %s -------" % data_parameters["basis"]
        if  action == "nlgd" and os.path.isfile("%s%s_A.%s_%d.communities" % (data_parameters['results_rep'], data_parameters['basis'], "nlgd", 20)):
            print "Skipping, done already"
            continue

        cleanCache(data_parameters)
        graph = readGraph(data_parameters['data_rep']+ bf[0], center=bf[-1], min_weight=data_parameters['min_weight'])
        nodes_to_labels = None
        print "Read edges [%d]..." % len(graph)
        
        if action != "ego":
            ndlb_dict = readLabels(data_parameters['data_rep']+ bf[1], graph=graph, min_cl=data_parameters["min_count"])
            if ndlb_dict.has_key(bf[-1]):
                del ndlb_dict[bf[-1]]
            if data_parameters.has_key("min_frac"):
                normalizeLabelsCounts(ndlb_dict, data_parameters)
            data_parameters["ndlb_dict"] = ndlb_dict
            nodes_to_labels = getNodesToLabels(graph, ndlb_dict, data_parameters)

        print "Read labels [%d]..." % len(nodes_to_labels)
        # print sum([len(i) for i in nodes_to_labels.values()])/len(nodes_to_labels)
        # pdb.set_trace()

        store_parameters = data_parameters
        for serie in series:
            data_parameters = dict(store_parameters)
            data_parameters.update(serie)
            data_parameters['basis'] += "_%s" % data_parameters["serie_suff"]
            
            ACTIONS[action](graph, nodes_to_labels, data_parameters)

if __name__ == "__main__":
    run(sys.argv)
