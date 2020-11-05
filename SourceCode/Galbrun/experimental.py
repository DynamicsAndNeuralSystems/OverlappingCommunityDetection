#!/usr/bin/python
import os, os.path, sys, re, copy, glob, itertools, datetime, math, random, numpy
from scipy.sparse.linalg import eigsh
from scipy.sparse import coo_matrix, spdiags
import pdb

########################################################
################# READING
def readParameters(fp=None):
    tf_dict = {"True":True, "False":False, "true":True, "false": False}
    noparse = ["externals", "ks"]

    ### pickle: -1 disabled, 0 store, 1 load
    ### log: - stdout, + auto filename, * both, other filename
    ### normalize comp: flip comparison edges A>B to A<B
    ### try simply: try finding simpler clauses when remapping
    ### try V: which shapes to try 1 /\, 2 \/, 3 -- positive indices: even if simple edge has enough support
    org_parameters = {'verbosity': 3, 'logfile': "+", "ego": -1, "center":"", "min_supp": 3,
                      'results_rep': "../results/", "data_rep":"../data/", "min_count": 1.0, "action": "mine",
                      'spread_weighted': 1, 'spread_outdiv': 0, 'spread_outweight': 0, 'spread_indiv': 1, 'spread_thres': 1/4.0,
                      "externals":"", "ks":"0", "output": "+", "min_weight": None, "type_support":"int",
                      "tri_neighf":0.5, "spread_mine":0, "spread_eval": 0}

    # 'limit_paths_itemsets': 1000, 'limit_graphs_itemsets': 1000,
    # 'limit_unit_instances': 0, 'limit_abs_instances': 10**5,

    if fp is not None:
        read_parameters = {}
        for line in fp:
            if re.match("^ *#", line) is not None:
                continue
            parts = line.split()
            if len(parts) > 1:
                if parts[0] in noparse:
                    tmp = " ".join(parts[1:])
                    read_parameters[parts[0]] = tmp
                elif parts[1] in tf_dict:
                    read_parameters[parts[0]] = tf_dict[parts[1]]
                else:
                    try:
                        read_parameters[parts[0]] = int(parts[1])
                    except ValueError:
                        try:
                            read_parameters[parts[0]] = float(parts[1])
                        except ValueError:
                            if "~" in parts[1]:
                                tmp = re.sub("~", os.path.expanduser("~"), parts[1])
                            else:
                                tmp = " ".join(parts[1:])
                            read_parameters[parts[0]] = tmp
        org_parameters.update(read_parameters)

    poss_externals = set(["dbp", "conga", "cmod", "metis", "maxc", "blda"]) # "bmlpa"
    org_parameters["externals"] = list(poss_externals.intersection([s.strip() for s in org_parameters["externals"].strip().split(",")]))
    org_parameters["ks"] = [int(s) for s in org_parameters["ks"].strip().split(",")]
    return org_parameters

def readGraph(filename, min_weight=None, center=None):
    graph = {}
    # with codecs.open(filename, encoding='utf-8', mode='r') as f:
    sep = None
    with open(filename, mode='r') as f:
        for line in f:
            if sep is None:
                tmp = re.match("^(?P<from>[^\t]*)\t(?P<to>[^\t]*)(\t(?P<count>[0-9.]*) *)?$", line.strip())
                if tmp is None:
                    sep = " "
                else:
                    sep = "\t"
            tmp = re.match("^(?P<from>[^"+sep+"]*)"+sep+"(?P<to>[^"+sep+"]*)("+sep+"(?P<count>[0-9.]*) *)?$", line.strip())
            if tmp is not None and tmp.group("from") != tmp.group("to"):
                nfrom, nto, ncount = tmp.group("from"), tmp.group("to"), None
                if tmp.group("count") is not None:
                    ncount = float(tmp.group("count"))
                if ncount >= min_weight and center not in [nfrom, nto]:
                    for (f,t) in [(nfrom, nto), (nto, nfrom)]:
                        if not graph.has_key(f):
                            graph[f] = {t:ncount}
                        else:
                            graph[f][t] = ncount
    return graph

def readLabels(lfilename, graph, min_cl=0):
    if lfilename.split(".")[-1] == "feat":
        return readLabelsF(lfilename, graph, min_cl)
    if lfilename.split(".")[-1] == "featlabels":
        return readLabelsL(lfilename, graph, min_cl)
    else:
        return readLabelsW(lfilename, graph, min_cl)

def readLabelsF(lfilename, graph=None, min_cl=0):
    ndlb_dict = {}
    names = []
    with open(lfilename+"names", mode='r') as f:
        for line in f:
            names.append(line.strip())

    with open(lfilename, mode='r') as f:
        for line in f:
            tmp = re.match("^(?P<node>[^ ]*) (?P<feats>[01 \t]*)$", line)
            if tmp is not None and ( graph is None or graph.has_key(tmp.group("node"))):
                # ndlb_dict[tmp.group("node")] = dict([("%d" % i,1) for i,v in enumerate(tmp.group("feats").strip().split()) if v != "0"])
                ndlb_dict[tmp.group("node")] = dict([(names[i],1) for i,v in enumerate(tmp.group("feats").strip().split()) if v != "0"])
    return ndlb_dict

def readLabelsL(lfilename, graph=None, min_cl=0):
    ndlb_dict = {}
    with open(lfilename, mode='r') as f:
        for line in f:
            tmp = re.match("^(?P<node>[^ ]*) (?P<term>.*)$", line)
            if tmp is not None and ( graph is None or graph.has_key(tmp.group("node"))):
                if ndlb_dict.has_key(tmp.group("node")):
                    ndlb_dict[tmp.group("node")][tmp.group("term").strip()] = 1
                else:
                    ndlb_dict[tmp.group("node")] = {tmp.group("term").strip(): 1}
    return ndlb_dict

def readLabelsW(lfilename, graph=None, min_cl=0):
    ndlb_dict = {}
    cnode = None
    with open(lfilename, mode='r') as f:
        for line in f:
            tmp = re.match("^(?P<node>[^\t]*)\t(?P<term>[^\t]*)\t(?P<count>\d*)$", line)
            if tmp is not None and ( int(tmp.group("count")) >= min_cl ) and ( graph is None or graph.has_key(tmp.group("node"))):
                if ndlb_dict.has_key(tmp.group("node")):
                    ndlb_dict[tmp.group("node")][tmp.group("term")] = int(tmp.group("count"))
                else:
                    ndlb_dict[tmp.group("node")] = {tmp.group("term"): int(tmp.group("count"))}
    return ndlb_dict

def readLabelsS(lfilename, graph=None, min_cl=0):
    stop_words = "algorithm analysy applicat approach base comput control data design gener informat method model multy network study system time".split()
    ndlb_dict = {}
    cnode = None
    with open(lfilename, mode='r') as f:
        for line in f:
            tmp = re.match("^(?P<node>[^\t]*)\t(?P<term>[^\t]*)\t(?P<count>\d*)$", line)
            if tmp is not None and ( int(tmp.group("count")) >= min_cl ) and ( graph is None or graph.has_key(tmp.group("node"))):
                if tmp.group("node") != cnode and cnode is not None:
                    sm = float(sum(ndlb_dict[cnode].values()))
                    keep = set(ndlb_dict[cnode].keys()) - remove
                    ndlb_dict[cnode] = dict([(k, ndlb_dict[cnode][k]) for k in keep if ndlb_dict[cnode][k]/sm >= data_parameters["min_frac"] and ndlb_dict[cnode][k] >= data_parameters["min_count"]])
                    cnode = tmp.group("node")

                if tmp.group("term") not in stop_words:
                    if ndlb_dict.has_key(tmp.group("node")):
                        ndlb_dict[tmp.group("node")][tmp.group("term")] = int(tmp.group("count"))
                    else:
                        ndlb_dict[tmp.group("node")] = {tmp.group("term"): int(tmp.group("count"))}
    return ndlb_dict

def getNodesToLabels(graph, ndlb_dict, data_parameters):
    nodes_to_labels = {}
    for node, lbd in ndlb_dict.items():
        nodes_to_labels[node] = set(lbd.keys())
    return nodes_to_labels

def readCandidates(lfilename, lsorted):
    candidates = []
    with open(lfilename, mode='r') as f:
        for line in f:
            candidates.append([lsorted[int(p)] for p in line.strip().split()[:-1]])
    return candidates

def readClusters(lfilename, uls):
    clusters = []
    with open(lfilename, mode='r') as f:
        for line in f:
            try:
                clusters.append(set([uls[int(pi)-1] for pi in line.split()]))
            except IndexError:
                clusters.append(set())
    return clusters

def readCLabels(lfilename):
    labelsets = []
    with open(lfilename, mode='r') as f:
        for line in f:
            labelsets.append(sorted([l.strip() for l in line.split()]))
    return labelsets

def cleanCache(data_parameters):
    for v in ["lsdict", "labels_to_nodes", "labels_to_edges", "edges_to_labels"]:
        if data_parameters.has_key(v):
            del data_parameters[v]

########################################################
################# PRINTING
def openMakeP(filename, mode="w"):
    d = os.path.dirname(filename)
    if not os.path.isdir(d):
        os.makedirs(d)
    return open(filename, mode)

def initResOut(data_parameters, external, k):
    if data_parameters["output"] != "-":
        if data_parameters["output"] == "+":
            fbasis = "%s%s.%s_%d" % (data_parameters['results_rep'], data_parameters['basis'], external, k)
        else:
            fbasis = data_parameters["output"]
        fo = openMakeP(fbasis+".communities", "w")
        fe =None
        # fe = openMakeP(fbasis+".cclust", "w")
    else:
        fo = sys.stdout
        fe = None
    if fe is not None:
        print "------------METHOD=%s K=%d" % (external, k)
    fo.write("------------METHOD=%s K=%d\n" % (external, k))
    return fo, fe

def labelsOut(nodes_to_labels, lmap, fname_mtx):
    fo = openMakeP(fname_mtx, "w")
    for e in sorted(nodes_to_labels.keys()):
        fo.write(" ".join([lmap[l] for l in nodes_to_labels[e]])+"\n")
    fo.close()

def edgesOut(graph, fname_grl, fname_grd):
    uls = sorted([k for k,vs in graph.items() if len(vs) > 0])
    fo = openMakeP(fname_grl, "w")
    fo.write("\n".join(uls))
    fo.close()

    lsdict = dict([(v,i+1) for i,v in enumerate(uls)])

    fo = openMakeP(fname_grd, "w")
    for k, vs in graph.items():
        fo.write(("\n".join(["%d %d" % (ldict[k], ldict[t]) for t in vs.keys()]))+"\n")
    fo.close()
    return uls

def graphOut(graph, fname_grl, fname_grd, weighted=False):
    uls = sorted([k for k,vs in graph.items() if len(vs) > 0])
    fo = openMakeP(fname_grl, "w")
    fo.write("\n".join(uls))
    fo.close()

    lsdict = dict([(v,i+1) for i,v in enumerate(uls)])

    fo = openMakeP(fname_grd, "w")
    if weighted:
        fo.write("%d %d 001\n" % (len(lsdict), sum([len(k) for k in graph.values()])/2))
    else:
        fo.write("%d %d\n" % (len(lsdict), sum([len(k) for k in graph.values()])/2))
    for k, vs in graph.items():
        if weighted:
            tmp = sorted([(lsdict[v], 1000*w) for v, w in vs.items()])
            fo.write((" ".join(["%d %d" % t for t in tmp]))+"\n")
        else:
            tmp = sorted([lsdict[v] for v, w in vs.items()])
            fo.write((" ".join(["%d" % t for t in tmp]))+"\n")
    fo.close()
    return uls

def externalOut(graph, nodes_to_labels, basis):
    nls = getNodeOrds(graph)

    #### GRAPH EDGES LIST AND METIS FORMAT
    lsdict = dict([(v,i+1) for i,v in enumerate(nls)])
    fe = openMakeP(basis+".ledges", "w")
    fo = openMakeP(basis+".metis", "w")
    fo.write("%d %d\n" % (len(lsdict), sum([len(k) for k in graph.values()])/2))
    for k, vs in graph.items():
        tmp = sorted([lsdict[v] for v, w in vs.items()])
        fo.write((" ".join(["%d" % t for t in tmp]))+"\n")
        for t in tmp:
            if t > lsdict[k]:
                fe.write("%d %d\n" % (lsdict[k], t))
    fo.close()
    fe.close()

    ### LABEL MATRIX
    lbls = set()
    for vs in nodes_to_labels.values():
        lbls.update(vs)
    lbls = sorted(lbls)

    ### LABEL MAP
    fo = openMakeP(basis+".lmap", "w")
    fo.write("\n".join(lbls))
    fo.close()

    lbdict = dict([(v,"%d" % (i+1)) for i,v in enumerate(lbls)])
    fo = openMakeP(basis+".lmatrix", "w")
    fo.write("%d\n" % len(nls))
    fo.write("%d\n" % len(lbls))
    for u in nls:
        fo.write("%d " % len(nodes_to_labels.get(u, [])))
        fo.write((" ".join(sorted([lbdict[l] for l in nodes_to_labels.get(u, [])])))+"\n")
    fo.close()
    return nls

def addGroup(ti, scr, lspec, edgeset, nodeset, labelset=None, covered=None, coveredn=None, nb_e=-1, nb_n=-1, all_edges=set(), graph=None, sorted_nodes=False):
    sorted_nodes=True
    if covered is not None:
        if type(covered) == dict:
            newedges = edgeset.difference(covered.keys())
        else:
            newedges = edgeset.difference(covered)
    if coveredn is not None:
        if type(covered) == dict:
            newnodes = nodeset.difference(coveredn.keys())
        else:
            newnodes = nodeset.difference(coveredn)
    nn = len(nodeset)

    if type(covered) == dict:
        for edge in newedges:
            covered[edge] = ti
        for node in nodeset:
            if coveredn.has_key(node):
                coveredn[node].append(ti)
            else:
                coveredn[node]= [ti]
    elif covered is not None:
        covered |= newedges
        coveredn |= newnodes

    tmp = ("*(%d)\ts=%f sL=%f d=%d dR=%f E=%d N=%d L=%d" % (ti, scr, lspec, density(edgeset, nn), densityR(edgeset, nn), len(edgeset), nn, len(labelset)))
    if covered is not None and coveredn is not None:
        tmp += (" Cd=%d CdR=%f CE=%d CN=%d" % (density(newedges, nn), densityR(newedges, nn), len(newedges), len(newnodes)))
        tmp += (" (%d/%d=%f, %d/%d=%f)"  % (len(covered), nb_e, len(covered)/nb_e, len(coveredn), nb_n, len(coveredn)/nb_n))
        tmp += (" cut_edges=%d overlap_nodes=%f"  % (len([e for e in all_edges - set(covered.keys()) if coveredn.has_key(e[0]) and coveredn.has_key(e[1])]), numpy.mean([len(v) for v in coveredn.values()])))
    if labelset is not None:
        tmp += "\n[%s]" % ", ".join(labelset)

    # if graph is not None:
    #     # tmp += ("\n (%s)\n" % (", ".join(["%s [%.2f]" % nif for nif in sorted([(n, sum(graph[n].values())) for n in nodeset], key=lambda x: x[1], reverse=True)])))
    #     tmp += ("\n (%s)\n" % (", ".join(sorted(nodeset))))
    if sorted_nodes:
        scores = dict([(n,0) for n in nodeset])
        for e in edgeset:
            scores[e[0]] += 1
            scores[e[1]] += 1

        tmp += ("\n (%s)\n" % (", ".join(sorted([k for k in scores.keys()], key=lambda x: (-scores[x], x))[:10])))

    elif graph is not None and len(nodeset) > 5:
        tmp += ("\n (%s ...)\n" % (", ".join(sorted(random.sample(nodeset,5)))))
    else:
        tmp += ("\n (%s)\n" % (", ".join(sorted(nodeset))))
    return tmp, newedges, newnodes

########################################################
################# MISC
def reverseNdlb(ndlb_dict, nodes=None):
    labels_to_nodes = {}
    if nodes is None:
        nodes = ndlb_dict.keys()
    for n in nodes:
        for vi in ndlb_dict.get(n,{}).keys():
            if labels_to_nodes.has_key(vi):
                labels_to_nodes[vi].add(n)
            else:
                labels_to_nodes[vi] = set([n])
    return labels_to_nodes

def reverseNdlbSets(ndlb_dict, nodes=None):
    labels_to_nodes = {}
    if nodes is None:
        nodes = ndlb_dict.keys()
    for n in nodes:
        for vi in ndlb_dict.get(n,set()):
            if labels_to_nodes.has_key(vi):
                labels_to_nodes[vi].add(n)
            else:
                labels_to_nodes[vi] = set([n])
    return labels_to_nodes


def reverseLabels(nodes_to_labels, nodes=None):
    labels_to_nodes = {}
    if nodes is None:
        nodes = nodes_to_labels.keys()
    for n in nodes:
        for vi in nodes_to_labels.get(n,[]):
            if labels_to_nodes.has_key(vi):
                labels_to_nodes[vi].add(n)
            else:
                labels_to_nodes[vi] = set([n])
    return labels_to_nodes

def getReversedLabels(nodes_to_labels, graph, data_parameters, store=True):
    if not data_parameters.has_key("labels_to_nodes") or not store:
        labels_to_nodes = reverseLabels(nodes_to_labels, graph.keys())
        if store:
            data_parameters["labels_to_nodes"] = labels_to_nodes
    else:
        labels_to_nodes = data_parameters["labels_to_nodes"]
    return labels_to_nodes

def getNodeset(edgeset):
    return set(itertools.chain.from_iterable(zip(*edgeset)))

def getNodeOrds(graph):
    return sorted([k for k,vs in graph.items() if len(vs) > 0]) + sorted([k for k,vs in graph.items() if len(vs) == 0])

# def getNodeId(node, graph, data_parameters):
#     if not data_parameters.has_key("lsdict"):
#         lsdict = dict([(v,"%d" % (i+1)) for i,v in enumerate(getNodeOrds(graph))])
#         data_parameters["lsdict"] = lsdict
#     else:
#         lsdict = data_parameters["lsdict"]
#     return lsdict[node]

### to print the node names rather than num id
def getNodeId(node, graph, data_parameters):
    return node

def projectGraph(nodes, graph, min_weight=None):
    subgraph = {}
    for node in nodes:
        subgraph[node] = dict([(nt, nc) for nt, nc in graph.get(node, {}).items() if nc >= min_weight])
        if len(subgraph[node]) == 0:
            del subgraph[node]
    return subgraph

def makeEdgesCands(graph, nodes_to_labels):
    labels_to_edges = {}
    edges_to_labels = {}
    for nf, nts in graph.items():
        for nt, c in nts.items():
            if nf == min(nf, nt):
                tlbs = nodes_to_labels.get(nf, set()) & nodes_to_labels.get(nt, set())
                edges_to_labels[(nf,nt)] = tlbs
                for l in tlbs:
                    if labels_to_edges.has_key(l):
                        labels_to_edges[l].add((nf,nt))
                    else:
                        labels_to_edges[l] = set([(nf,nt)])
    return labels_to_edges, edges_to_labels

def getEdgesCands(nodes_to_labels, graph, data_parameters, store=True):
    if not data_parameters.has_key("labels_to_edges") or not store:
        labels_to_edges, edges_to_labels = makeEdgesCands(graph, nodes_to_labels)
        if store:
            data_parameters["labels_to_edges"], data_parameters["edges_to_labels"] = labels_to_edges, edges_to_labels
    else:
        labels_to_edges, edges_to_labels = data_parameters["labels_to_edges"], data_parameters["edges_to_labels"]
    return labels_to_edges, edges_to_labels

def cleanLabels(ndlb_dict, data_parameters):
    for a in ndlb_dict.keys():
        sm = float(sum(ndlb_dict[a].values()))
        ndlb_dict[a] = dict([k for k in ndlb_dict[a].items() if k[1]/sm >= data_parameters["min_frac"] and k[1] >= data_parameters["min_count"]])

def normalizeLabelsCounts(ndlb_dict, data_parameters):
    remove = set()
    if data_parameters.has_key("stop_frac"):
        nbstop = len(ndlb_dict)*data_parameters["stop_frac"]
        labels_to_nodes = reverseNdlb(ndlb_dict)
        remove = set([l for l,v in labels_to_nodes.items() if len(v) >= nbstop])

    for a in ndlb_dict.keys():
        sm = float(sum(ndlb_dict[a].values()))
        keep = set(ndlb_dict[a].keys()) - remove
        ndlb_dict[a] = dict([(k, ndlb_dict[a][k]) for k in keep if ndlb_dict[a][k]/sm >= data_parameters["min_frac"] and ndlb_dict[a][k] >= data_parameters["min_count"]])


def updateSpreads(graph, ndlb_dict, labelset, cnodes, labels_to_nodes, labels_to_edges, labels_to_nodes_unspread, data_parameters):
#    print labelset
#    pdb.set_trace()
    for label in labelset:
        up_nodes = set([node for node in labels_to_nodes[label] if (node not in cnodes) and len(labels_to_nodes_unspread[label].intersection(graph.get(node, {}).keys()).difference(cnodes)) > len(graph.get(node, {}).keys())/data_parameters["spread_thres"]])
        up_edges = set([edge for edge in labels_to_edges[label] if (edge[0] in up_nodes) and (edge[1] in up_nodes)])
        labels_to_nodes[label] = up_nodes
        labels_to_edges[label] = up_edges

def spreadLabels(graph, ndlb_dict, data_parameters):
    all_nodes = set(graph.keys()).union(ndlb_dict.keys())
    spread_ndlb_dict = {}
    for node in all_nodes:
        neighs = graph.get(node, {})
        spread_ndlb_dict[node] = dict(ndlb_dict.get(node, {}))
        cand_labels = {}
        for neigh, nw in neighs.items():
            for label, lw in ndlb_dict.get(neigh, {}).items():
                if data_parameters["spread_weighted"] == 0:
                    lw = 1.0
                if data_parameters["spread_outdiv"] == 1:
                    lw = lw / float(len(graph.get(neigh, {})))
                if data_parameters["spread_outweight"] == 1:
                    lw = lw / float(sum(ndlb_dict[neigh].values()))
                ### options for lw: 1, 1.0/len(graph.get(node, {}))
                if cand_labels.has_key(label):
                    cand_labels[label] += lw
                else:
                    cand_labels[label] = lw

        for label, lw in cand_labels.items():
            if data_parameters["spread_indiv"] == 1:
                lw = lw / float(len(graph.get(node, {})))

        # lbls = set(cand_labels.keys()).difference(spread_ndlb_dict[node])
        # for label in lbls:
        #     if re.search("c:[0-9]", node) and re.search("C[0-9]", label):
        #         print node, label
        #         pdb.set_trace()

            if lw >= data_parameters["spread_thres"]:
                if label not in spread_ndlb_dict[node]:
                    spread_ndlb_dict[node][label] = 1
    return spread_ndlb_dict

def extractSubs(graph, nodes_to_labels, data_parameters):
    min_weight = data_parameters["min_weight"]
    min_cl = 0
    nodes = set()
    if data_parameters["ego"] > 0:
        ego = data_parameters["ego"]
    # for ego in [2]:
        for centert in data_parameters["centers"].split(";"):
            center = centert.strip()

            if ego > 0 and center is not None and graph.has_key(center):
                print "Filtering %s %d" % (center, ego)
                ### GATHER THE EGO NET
                last_nodes = set([center])
                new_nodes = set()
                for i in range(ego+1):
                    for node in last_nodes:
                        if graph.has_key(node):
                            new_nodes.update([t for t,c in graph[node].items() if c >= min_weight])
                    nodes |= last_nodes
                    last_nodes = new_nodes - nodes
                    new_nodes = set()

            if "ego" in data_parameters["action"]:
                fo = openMakeP(data_parameters['data_rep'] + ("wegonet_%s%d-%.1f.graph" % (center.replace(" ",""), ego, min_weight)), "w")
                for node in nodes:
                    for nt, nc in graph.get(node, {}).items():
                        if nt in nodes and nc >= min_weight:
                            fo.write("%s\t%s\t%f\n" % (node, nt, nc))
                fo.close()


                if data_parameters.has_key('labels_file'):
                    fo = openMakeP(data_parameters['data_rep'] + ("wegonet_%s%d-%.1f.labels" % (center.replace(" ",""), ego, min_weight)), "w")
                    with open(data_parameters['data_rep'] + data_parameters['labels_file'], mode='r') as f:
                        for line in f:
                            tmp = re.match("^(?P<node>[^\t]*)\t(?P<term>[^\t]*)\t(?P<count>\d*)$", line)
                            if tmp is not None and ( int(tmp.group("count")) >= min_cl ) and ( nodes is None or tmp.group("node") in nodes):
                                fo.write(line)
                    fo.close()
                nodes = set()

########################################################
################# SCORING
def score(nb_edges, nb_nodes=None, new_edges=None):
    ### new_edges is not used in this score, only here for compatibility
    if nb_nodes is None: ## in this case provided the actual edges and not just their number
        nb_nodes = len(getNodeset(nb_edges))
    if nb_edges > 0:
        return nb_nodes/(2.0*nb_edges)
    return float("Inf")

def score1(nb_edges, nb_nodes=None, new_edges=None):
    if nb_nodes is None: ## in this case provided the actual edges and not just their number
        nb_nodes = len(getNodeset(nb_edges))
    if type(nb_edges) is not int: ## if actual edges provided turn to number
        nb_edges = len(nb_edges)
    if new_edges is None:
        new_edges = nb_edges
    if nb_edges > 0:
        return ((nb_nodes*(nb_nodes-1))/2.0-nb_edges)/new_edges
    return float("Inf")

def density(edgeset, nbN=None):
    if nbN is None:
        nbN = len(getNodeset(edgeset))
    r = (nbN*(nbN-1))-2*len(edgeset)
    return r

def densityR(edgeset, nbN=None):
    if nbN is None:
        nbN = len(getNodeset(edgeset))
    if nbN > 1:
        r = 2.0*len(edgeset)/(nbN*(nbN-1))
    else:
        r = float("Inf")
    return r

def ratio(edgeset, nbN=None):
    if nbN is None:
        nbN = len(getNodeSet(edgeset))
    r = 2.0*len(edgeset)/(nbN*(nbN-1)- 2.0*len(edgeset))
    return r

def computeLabelsSpec(labels_to_nodes, nodeset, labelset, data_parameters):
    if len(labelset) == 0:
        return 0.0
    supp, suppe = getSupportInter(labelset, labels_to_nodes)
    return len(nodeset.intersection(supp))/(1.0*len(nodeset.union(supp)))

def chooseLabels(labels_to_nodes, nodeset, data_parameters):
    candidates = sorted([l for l, supp in labels_to_nodes.items() if len(supp & nodeset) > 0])
    best = (0, None)
    for l in candidates:
        if min(len(nodeset), len(labels_to_nodes[l]))/float(max(len(nodeset), len(labels_to_nodes[l]))) > best[0]:
            t = len(nodeset.intersection(labels_to_nodes[l]))/(1.0*len(nodeset.union(labels_to_nodes[l])))
            if t > best[0]:
                best = (t, l)
    if best[1] is None:
        return [], 0
    top = (best[0], [best[1]], set(labels_to_nodes[best[1]]))
    while top[1] is not None:
        best = top
        top = (best[0], None, None)
        for l in candidates:
            supp = labels_to_nodes[l].intersection(best[2])
            if min(len(nodeset), len(supp))/float(max(len(nodeset), len(supp))) > top[0]:
                t = len(nodeset.intersection(supp))/(1.0*len(nodeset.union(supp)))
                if t > top[0]:
                    top = (t, [l]+best[1], supp)
    return best[1], best[0]


def chooseLabelsTheta(labels_to_nodes, nodeset, theta, feats, data_parameters):
    candidates = [feats[k] for (k, v) in sorted(enumerate(theta[1:]), key=lambda x:x[1]) if v > 0]
    best = [0, [], set(nodeset)]
    while len(candidates) > 0:
        l = candidates.pop()
        suppL = labels_to_nodes.get(l,[])
        if min(len(best[2]), len(suppL))/float(max(len(best[2]), len(suppL))) > best[0]:
            t = len(best[2].intersection(suppL))/(1.0*len(best[2].union(suppL)))
            if t > best[0]:
                best[1].append(l)
                best[0] = t
                best[2].intersection_update(suppL)
            else:
                candidates = []
        else:
            candidates = []
    pdb.set_trace()
    print best
    return best[1], best[0]


def getSupportInter(labelset, labels_to_nodes, labels_to_edges=None):
    cnodes = set(labels_to_nodes.get(labelset[0], set()))
    for l in labelset[1:]:
        cnodes &= labels_to_nodes.get(l, set())

    if labels_to_edges is not None:
        cedges = set(labels_to_edges.get(labelset[0], set()))
        for l in labelset[1:]:
            cedges &= labels_to_edges.get(l, set())
    else:
        cedges = None
    return cnodes, cedges

def getSupportMajority(labelset, nodes_to_labels, edges_to_labels=None, all_edges=None):
    if all_edges is None and edges_to_labels is not None:
        all_edges = set(edges_to_labels.keys())
    lss = set(labelset)
    mini = len(labelset)/2+1
    try:
        cnodes = set([node for (node, ls) in nodes_to_labels.items() if len(lss.intersection(ls)) >= mini])
    except AttributeError:
        pdb.set_trace()
        print "A"
        cnodes = set()
    if all_edges is not None:
        cedges = set([edge for edge in all_edges if edge[0] in cnodes and edge[1] in cnodes])
    else:
        cedges = None
    return cnodes, cedges

def getSupport(stype, labelset, labels_to_nodes=None, labels_to_edges=None, nodes_to_labels=None, edges_to_labels=None, all_edges=None):
    if stype == "maj" and (labels_to_edges is None or len(labelset) > 2):
        return getSupportMajority(labelset, nodes_to_labels, edges_to_labels, all_edges)
    else:
        return getSupportInter(labelset, labels_to_nodes, labels_to_edges)

def getSizeLInter(node, labelset, nodes_to_labels):
    return len(nodes_to_labels.get(node, set()).intersection(labelset))

def improveSubgraph(graph, edgeset, covered, scoref=densityR, f=1):
    cdR = scoref(edgeset)
    odR = cdR
    nodes = getNodeset(edgeset)
    neighs = set()
    for node in nodes:
        neighs.update(graph[node].keys())
    added = []
    nedges = set()
    neighs = sorted([(n, nodes & set(graph[n].keys())) for n in neighs - nodes], key= lambda x:len(x[1]))
    while len(neighs) > 0 and len(added) < 0.1*len(nodes):
        n, nds = neighs.pop()
        eds = set([tuple(sorted([n, o])) for o in nds])
        tdR = scoref(edgeset | nedges | eds)
        if f*tdR > f*cdR:
            cdR = tdR
            nedges.update(eds)
            added.append(n)
    #print "Improved with %d nodes from %f to %f" % (len(added), odR, cdR)
    return nedges, added

########################################################
################# MINING
def densestMaj(nodes_to_labels, labels_to_nodes, all_edges, assigned=None, sizes=None, debug_add=None):
    distsN = [set(nodes_to_labels.keys())]
    onemissE = [[], [], list(all_edges)]
    uncoveredE, coveredE = [], []
    minC = 1
    scores = []
    labelset = []
    cand_labels = labels_to_nodes.keys()
    while True:
        if debug_add is not None:
            debug_add["lbls"] = labelset
        sc, label = bestLabelMaj(cand_labels, labels_to_nodes, minC, distsN, onemissE, coveredE, assigned, sizes, debug_add=debug_add)
        if label is None:
            break
        cand_labels.remove(label)
        labelset.append(label)
        scores.append(sc)
        old_minC = minC
        minC = (len(labelset)+3)/2
        if debug_add is not None: ## DEBUG ##
            cN = sum([len(c) for c in distsN[old_minC:]])
            cE = len(coveredE)
            N, E = NandELabelMaj(labels_to_nodes[label], old_minC, distsN, onemissE)
            cnodes, cedges = getSupport("maj", labelset, nodes_to_labels=nodes_to_labels, all_edges=all_edges)
            if cN+N != len(cnodes) or cE+E != len(cedges):
                print "TOP %s N:%d+%d=%d E:%d+%d=%d" % (labelset, cN, N, cN+N, cE, E, cE+E)
                print "VS.", len(cnodes), len(cedges)
                pdb.set_trace()
        distsN, onemissE, coveredE, uncoveredE = updateLabelMaj(labelset, labels_to_nodes, nodes_to_labels,
                                                                old_minC, minC, distsN, onemissE, coveredE, uncoveredE)

    return labelset, scores


def updateLabelMaj(labelset, labels_to_nodes, nodes_to_labels, old_minC, minC, distsN, onemissE, coveredE, uncoveredE):
    ### update list of nodes by intersection size with label set
    up_distsN = [set()]
    for distN in distsN:
        up_distsN[-1].update(distN.difference(labels_to_nodes[labelset[-1]]))
        up_distsN.append(distN.intersection(labels_to_nodes[labelset[-1]]))

    ### update list of edges with one label missing on either or both ends of the edge
    ## minC = old_minC + 0 or 1
    if minC == old_minC:
        up_coveredE = list(coveredE)
        up_uncoveredE = []
        up_onemissE = [[], [], []]
        retest = uncoveredE
        for i in [0,1]:
            for edge in onemissE[i]:
                if edge[i] in labels_to_nodes[labelset[-1]]:
                    up_coveredE.append(edge)
                else:
                    up_onemissE[i].append(edge)

        for i in [2]:
            for edge in onemissE[i]:
                if edge[0] in labels_to_nodes[labelset[-1]] and edge[1] in labels_to_nodes[labelset[-1]]:
                    up_coveredE.append(edge)
                elif edge[0] in labels_to_nodes[labelset[-1]]:
                    up_onemissE[1].append(edge)
                elif edge[1] in labels_to_nodes[labelset[-1]]:
                    up_onemissE[0].append(edge)
                else:
                    up_onemissE[2].append(edge)

    else:
        up_coveredE = []
        up_uncoveredE = []
        up_onemissE = [[], [], []]
        retest = uncoveredE + coveredE
        for cE in onemissE:
            retest.extend(cE)

    for edge in retest:
        c0, c1 = getSizeLInter(edge[0], labelset, nodes_to_labels), getSizeLInter(edge[1], labelset, nodes_to_labels)
        if c0 < minC-1 or c1 < minC-1:
            up_uncoveredE.append(edge)
        elif c0 >= minC and c1 >= minC:
            up_coveredE.append(edge)
        elif c0 == minC-1:
            if c1 == minC-1:
                up_onemissE[2].append(edge)
            else:
                up_onemissE[0].append(edge)
        else:
            up_onemissE[1].append(edge)

    return up_distsN, up_onemissE, up_coveredE, up_uncoveredE

def bestLabelMaj(cand_labels, labels_to_nodes, minC, distsN, onemissE, coveredE, assigned=None, sizes=None, debug_add=None):
    top = (float("Inf"), None)
    cN = sum([len(c) for c in distsN[minC:]])
    acE = len(coveredE)
    for label in cand_labels:
        N, E = NandELabelMaj(labels_to_nodes[label], minC, distsN, onemissE, assigned, sizes, cN)
        if assigned is not None:
            cE = len([c for c in coveredE if sizes[assigned.get(c, None)] > N+cN])
        else:
            cE = acE
        if cE+E > 0:
            sc = score(E+cE,N+cN)

            if debug_add is not None:
                cnodes, cedges = getSupport("maj", [label]+debug_add["lbls"], debug_add["ln"], debug_add["le"], nodes_to_labels=debug_add["nl"], edges_to_labels=debug_add["el"])
                ce = 0
                for e in cedges:
                    if sizes[assigned.get(e, None)] > len(cnodes):
                        ce +=1
                if score(ce, len(cnodes)) != sc:
                    N, E = NandELabelMaj(labels_to_nodes[label], minC, distsN, onemissE, assigned, sizes, cN)
                print "COMP SCORES", score(ce, len(cnodes)), sc

            if sc < top[0]:
                top = (sc, label)
    return top

def NandELabelMaj(label_nsupp, minC, distsN, onemissE, assigned=None, sizes=None, cN=0):
    N = len(label_nsupp & distsN[minC-1])
    if assigned is None:
        E  = sum([v_0 in label_nsupp for (v_0, v_1) in onemissE[0]])
        E += sum([v_1 in label_nsupp for (v_0, v_1) in onemissE[1]])
        E += sum([v_0 in label_nsupp and v_1 in label_nsupp for (v_0, v_1) in onemissE[2]])
    else:
        E  = sum([v_0 in label_nsupp for (v_0, v_1) in onemissE[0] if sizes[assigned.get((v_0, v_1), None)] > N+cN])
        E += sum([v_1 in label_nsupp for (v_0, v_1) in onemissE[1] if sizes[assigned.get((v_0, v_1), None)] > N+cN])
        E += sum([v_0 in label_nsupp and v_1 in label_nsupp for (v_0, v_1) in onemissE[2] if sizes[assigned.get((v_0, v_1), None)] > N+cN])
    return N, E


########################################################
################# MINING TOOLS

def parseClustersL(graph, nodes_to_labels, candidates, data_parameters, labelsets=None, fo=None, fc=None):
    ##### PARSE CLUSTERS BASED ONLY ON LABELS RECONSTRUCT SUPPORT!!
    min_supp = data_parameters["min_supp"]
    if fo is None:
        fo = sys.stdout

    nodes = set(graph.keys())
    labels_to_nodes = getReversedLabels(nodes_to_labels, graph, data_parameters)
    labels_to_edges, edges_to_labels = getEdgesCands(nodes_to_labels, graph, data_parameters)

    ticE = datetime.datetime.now()
    supps, suppsn, ncands = [], [], []
    while len(labelsets) > 0:
        cand = labelsets.pop()
        ndc, edgc = getSupport(data_parameters["type_support"], cand, labels_to_nodes, labels_to_edges, nodes_to_labels=nodes_to_labels, edges_to_labels=edges_to_labels)
        if len(ndc) > min_supp:
            supps.append(edgc)
            suppsn.append(ndc)
            ncands.append(cand)
    candidates = ncands

    cis = range(len(candidates))
    cis.sort(key= lambda x:(len(suppsn[x]),candidates[x]))
    covered = {}
    coveredn = {}
    nb_n = float(len(graph))
    nb_e = float(len(edges_to_labels))
    all_edges = set(edges_to_labels.keys())
    for ti, cid in enumerate(cis):
        nodeset = suppsn[cid]
        edgeset = supps[cid]
        labelset = candidates[cid]
        lspec = computeLabelsSpec(labels_to_nodes, nodeset, labelset, data_parameters)
        tmp, ne, nn = addGroup(ti, score(len(edgeset), len(nodeset)), lspec, edgeset, nodeset, labelset=labelset, covered=covered, coveredn=coveredn, nb_e=nb_e, nb_n=nb_n, all_edges=all_edges, graph=graph)
        fo.write(tmp+"\n")
        if fc is not None:
            fc.write("%d\t%d\t%s\t%s\n" % ( len(nodeset), len(labelset),
                                            ", ".join(sorted(nodeset)),
                                            ", ".join(sorted(labelset))))

def selectGreedy(graph, nodes_to_labels, candidates, data_parameters, fo=None, fc=None):
    min_supp = data_parameters["min_supp"]
    if fo is None:
        fo = sys.stdout

    nodes = set(graph.keys())
    labels_to_nodes = getReversedLabels(nodes_to_labels, graph, data_parameters)
    labels_to_edges, edges_to_labels = getEdgesCands(nodes_to_labels, graph, data_parameters)

    ticE = datetime.datetime.now()
    if type(candidates) == dict:
        candidates, supps = zip(*candidates.items())
        suppsn = []
        while len(suppsn) < len(candidates):
            nds = getNodeset(supps[len(suppsn)])
            if len(nds) > min_supp:
                suppsn.append(nds)
            else:
                supps.pop(len(suppsn))
                candidates.pop(len(suppsn))

    else:     ### generate supports for candidates if not provided
        supps, suppsn, ncands = [], [], []
        while len(candidates) > 0:
            cand = candidates.pop()
            ndc, edgc = getSupport(data_parameters["type_support"], cand, labels_to_nodes, labels_to_edges, nodes_to_labels=nodes_to_labels, edges_to_labels=edges_to_labels)

            if len(ndc) > min_supp:
                supps.append(edgc)
                suppsn.append(ndc)
                ncands.append(cand)
        candidates = ncands

    ### start greedy selection
    cis = range(len(candidates))
    cis.sort(key= lambda x:(len(suppsn[x]),candidates[x]))
    nb_e = float(len(edges_to_labels))
    nb_n = float(len(nodes))
    all_edges = set(edges_to_labels.keys())
    covered = {}
    coveredn = {}
    resultsc = 0
    top = (1, 1, 1, 1)
    while top[1] is not None:
        top = (float("Inf"), None, None)
        for ci in cis:
            if len(suppsn[ci].difference(coveredn.keys())) > min_supp: # and r < 1:
                r = score(len(supps[ci].difference(covered)), len(suppsn[ci]), len(supps[ci]))
                if r < top[0]: # and r != 1:
                    top = (r, ci)

        if top[1] is not None:
            lspec = computeLabelsSpec(labels_to_nodes, suppsn[top[1]], candidates[top[1]], data_parameters)
            tmp, ne, nn = addGroup(resultsc, top[0], lspec, supps[top[1]], suppsn[top[1]], labelset=candidates[top[1]], covered=covered, coveredn=coveredn, nb_e=nb_e, nb_n=nb_n, all_edges=all_edges, graph=graph)
            resultsc +=1
            fo.write(tmp+"\n")
            if fc is not None:
                fc.write("%d\t%d\t%s\t%s\n" % ( len(suppsn[top[1]]), len(candidates[top[1]]),
                                                ", ".join(sorted(suppsn[top[1]])),
                                                ", ".join(sorted(candidates[top[1]]))))
            cis.remove(ci)
    print "Done in %s" % (datetime.datetime.now() - ticE)

def getCandidatesFromLO(graph, labels_to_nodes, labels_to_edges, lo, data_parameters,
                        nodes_to_labels=None, edges_to_labels=None):
    start = 0
    mul = 0
    candidates = []
    while start < len(lo):
        lset = set([lo[start]])
        nsupp = labels_to_nodes[lo[start]]
        esupp = labels_to_edges[lo[start]]
        #if score(len(esupp), len(nsupp)) < 3:
        candidates.append({"nodeset":set(nsupp), "edgeset": set(esupp), "labelset": set(lset)})
        end = start+1
        while end < len(lo):
            esupp &= labels_to_edges[lo[end]]
            if len(esupp) < 2:
                end = len(lo)
            else:
                nsupp &= labels_to_nodes[lo[end]]
                lset.add(lo[end])
                #if score(len(esupp), len(nsupp)) < 3:
                mul += 1
                if data_parameters["type_support"] != "int":
                    nodesetB, edgesetB = getSupport(data_parameters["type_support"], list(lset), labels_to_nodes, labels_to_edges,
                                                    nodes_to_labels=nodes_to_labels, edges_to_labels=edges_to_labels)
                    nsupp = nodesetB
                    esupp = edgesetB
                # if nodesetB != nsupp or edgesetB != esupp:
                #     pdb.set_trace()
                #     print len(nsupp)
                candidates.append({"nodeset":set(nsupp), "edgeset": set(esupp), "labelset": set(lset)})
                # candidates.append({"nodeset":set(nodesetB), "edgeset": set(edgesetB), "labelset": set(lset)})
            end += 1
        start +=1
    # print "Returning %d candidates, %d multiple, out of %d labels" % (len(candidates), mul, len(lo))
    return candidates

def getLabelsOrdering(graph, labels_to_nodes, data_parameters):
    dM, labels = makeLabelDistMatrix(labels_to_nodes, data_parameters)
    if dM is not None:
        try:
            fiedo = getFiedlerOrd(dM)
        except Exception:
            fiedo = range(len(labels))
        return [labels[i] for i in fiedo], [-dM[fiedo[i],fiedo[i+1]] for i in range(len(fiedo)-1)]
    return labels, [0 for i in range(len(labels)-1)]

def makeLabelDistMatrix(labels_to_items, data_parameters):
    labels = sorted(labels_to_items.keys())
    vij = []
    for lAi, lA in enumerate(labels):
        vij += [(len(labels_to_items[lA] & labels_to_items[lB])/float(len(labels_to_items[lA] | labels_to_items[lB])), lAi, lBi)
               for lBi, lB in enumerate(labels[:lAi]) if len(labels_to_items[lA] & labels_to_items[lB]) > 0]
    if len(vij) > 0:
        V,I,J = zip(*vij)
        dM = coo_matrix((-numpy.array(V+V), (I+J,J+I)), shape=(len(labels),len(labels)) ).tocsr()
        return dM + spdiags(-dM.sum(0), [0], len(labels), len(labels), "csr"), labels
    else:
        return None, labels

def getFiedlerOrd(dM):
    eval, evect = eigsh(dM, 2)
    return numpy.argsort(evect[:,1])

def mineTriangleNode(node, graph, nodes_to_labels, edges_to_labels, labels_to_nodes, labels_to_edges, data_parameters, seen_labels=None):
    neighs_one = set(graph[node].keys())
    if len(neighs_one) == 1:
        return
    neighs_close = set([n for n in set(itertools.chain(*[graph[n].keys() for n in neighs_one])).difference(neighs_one)
                        if len(neighs_one.intersection(graph[n].keys())) > len(graph[n].keys())*data_parameters["tri_neighf"]])
    neighs_close.update(neighs_one)
    # inter_labels = set(nodes_to_labels[node]).intersection(*[nodes_to_labels[n] for n in neighs_close])
    neighs_close.add(node)
    labelset, lspec = chooseLabels(labels_to_nodes, neighs_close, data_parameters)
    labelset = tuple(labelset)
    # print node, len(neighs_close), inter_labels, labelset
    # pdb.set_trace()
    if len(labelset) > 0 and seen_labels is not None and labelset not in seen_labels:
        seen_labels.add(labelset)
        nodesetB, edgesetB = getSupport(data_parameters["type_support"], labelset, labels_to_nodes, labels_to_edges,
                                        nodes_to_labels=nodes_to_labels, edges_to_labels=edges_to_labels)
        if len(edgesetB) > 0:
            return {"labelset": labelset, "nodeset": nodesetB, "edgeset": edgesetB}

### greedy selection of candidates {"labelset", "nodeset", "edgeset"}, possibly stealing edges (CK)
def selectCandidates(candidates, graph, data_parameters, external, ticE):
    ks = data_parameters["ks"]
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
    selected = []
    ti = -1
    while len(candidates) > 0:
        select = sorted([(score(len([e for e in candidate["edgeset"] if assE[e][0] > len(candidate["nodeset"])]),
                                len(candidate["nodeset"])), ci) for ci, candidate in enumerate(candidates)])[0]
        if select[0] < float('inf'):
            ti += 1
            selected.append(candidates.pop(select[1]))
            for e in selected[-1]["edgeset"]:
                if assE[e][0] > len(selected[-1]["nodeset"]):
                    assE[e] = (len(selected[-1]["nodeset"]), len(selected)-1)
        else:
            candidates = []

        if ti+1 >= max(ks):
            candidates = []

        if ti+1 in ks or stop:
            for i in range(ti+1):
                selected[i]["edgeset"] = set()
            for edge, v in assE.items():
                if v[1] != None:
                    selected[v[1]]["edgeset"].add(edge)
            printSelected(selected, graph, data_parameters, external, ticE, nb_n, nb_e, all_edges)

def printSelected(selected, graph, data_parameters, external, ticE, nb_n, nb_e, all_edges):
    idsort = sorted(range(len(selected)), key=lambda x: len(selected[x]["nodeset"]))

    covered = {}
    coveredn = {}

    fo, fe = initResOut(data_parameters, external, len(selected))

    for ti, tt in enumerate(idsort):
        tmp, ne, nn = addGroup(ti, score(len(selected[tt]["edgeset"]), len(selected[tt]["nodeset"])), 1, selected[tt]["edgeset"], selected[tt]["nodeset"], labelset=selected[tt]["labelset"], covered=covered, coveredn=coveredn, nb_e=nb_e, nb_n=nb_n, all_edges=all_edges, graph=graph)
        fo.write(tmp+"\n")
        if fe is not None:
            # fe.write(" ".join(sorted([getNodeId(n, graph, data_parameters) for n in nodesets[tt]])) + "\n")
            fe.write("%d\t%d\t%s\t%s\n" % ( len(selected[tt]["nodeset"]), len(selected[tt]["labelset"]),
                                           " ".join(sorted(selected[tt]["nodeset"])),
                                           " ".join(sorted(selected[tt]["labelset"]))))

    fo.write("Done in %s\n" % (datetime.datetime.now() - ticE))
    if fe is not None:
        fo.close()
        fe.close()

def getEdgeset(nodes, graph):
    edgeset = set()
    for nodeA in nodes:
        for nodeB in graph.get(nodeA, {}).keys():
            if nodeA < nodeB and nodeB in nodes:
                edgeset.add((nodeA, nodeB))
    return edgeset
