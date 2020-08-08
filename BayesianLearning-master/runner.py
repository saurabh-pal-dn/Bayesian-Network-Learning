import BayesianNetwork as bn
import numpy as np
import pandas as pd
import time


def setupNetwork(alarm, record):
    print("0. Reading the network")
    Alarm = bn.readNetwork(alarm)
    print("1. Setting up Markov Blanket for Network")
    Alarm.setMb()
    # print("2. Printing the network")
    # bn.printNetwork(Alarm)
    print("3. Getting data of records")
    df = bn.getData(record)
    print("4. Initializing parameters....")
    init_parameters(df, Alarm)
    print("5. Getting missing Data")
    missingIndex = getMissingData(df)
    return Alarm, df, missingIndex


def getMissingData(df):
    miss = []
    for index in range(0, df.shape[0]):
        row = df.loc[index, :]
        for j in range(0, df.shape[1]):
            if row[j] == -71:
                miss.append(j)
    return miss


def init_parameters(df, bn):
    #using laplace smoothing for count
    N = df.shape[0]
    currentIteration = 0
    for node in bn.PresGraph.values():
        parents = bn.getParents(node)
        nParents = len(parents)
        currentIteration += 1
        if nParents == 0:
            v0 = []
            counts = []
            for p0 in range(0, node.nvalues):
                a = df[node.NodeName] == p0
                counts.append(pd.DataFrame(df[a]).shape[0] + 1)
                v0.append(p0)
            cptDf = pd.DataFrame({
                node.NodeName: v0,
                'p': np.ones(len(counts)) * (-1),
                'counts': counts
            })
            node.setCPTData(cptDf)

        elif nParents == 1:
            v0 = []
            v1 = []
            counts = []
            for p0 in range(0, parents[0].nvalues):
                a = df[parents[0].NodeName] == p0
                for p1 in range(0, node.nvalues):
                    b = df[node.NodeName] == p1
                    counts.append(pd.DataFrame(df[a & b]).shape[0] + 1)
                    v0.append(p0)
                    v1.append(p1)
            cptDf = pd.DataFrame({
                node.NodeName: v1,
                parents[0].NodeName: v0,
                'p': np.ones(len(counts)) * (-1),
                'counts': counts
            })
            node.setCPTData(cptDf)

        elif nParents == 2:
            v0 = []
            v1 = []
            v2 = []
            counts = []
            for p0 in range(0, parents[0].nvalues):
                a = df[parents[0].NodeName] == p0
                for p1 in range(0, parents[1].nvalues):
                    b = df[parents[1].NodeName] == p1
                    for p2 in range(0, node.nvalues):
                        c = df[node.NodeName] == p2
                        counts.append(pd.DataFrame(df[a & b & c]).shape[0] + 1)
                        v0.append(p0)
                        v1.append(p1)
                        v2.append(p2)
            cptDf = pd.DataFrame({
                node.NodeName: v2,
                parents[1].NodeName: v1,
                parents[0].NodeName: v0,
                'p': np.ones(len(counts)) * (-1),
                'counts': counts
            })
            # print(node.NodeName, cptDf)
            node.setCPTData(cptDf)

        elif nParents == 3:
            v0 = []
            v1 = []
            v2 = []
            v3 = []
            counts = []
            for p0 in range(0, parents[0].nvalues):
                a = df[parents[0].NodeName] == p0
                for p1 in range(0, parents[1].nvalues):
                    b = df[parents[1].NodeName] == p1
                    for p2 in range(0, parents[2].nvalues):
                        c = df[parents[2].NodeName] == p2
                        for p3 in range(0, node.nvalues):
                            d = df[node.NodeName] == p3
                            counts.append(
                                pd.DataFrame(df[a & b & c & d]).shape[0] + 1)
                            v0.append(p0)
                            v1.append(p1)
                            v2.append(p2)
                            v3.append(p3)
            cptDf = pd.DataFrame({
                node.NodeName: v3,
                parents[2].NodeName: v2,
                parents[1].NodeName: v1,
                parents[0].NodeName: v0,
                'p': np.ones(len(counts)) * (-1),
                'counts': counts
            })
            node.setCPTData(cptDf)

        elif nParents == 4:
            v0 = []
            v1 = []
            v2 = []
            v3 = []
            v4 = []
            counts = []
            for p0 in range(0, parents[0].nvalues):
                a = df[parents[0].NodeName] == p0
                for p1 in range(0, parents[1].nvalues):
                    b = df[parents[1].NodeName] == p1
                    for p2 in range(0, parents[2].nvalues):
                        c = df[parents[2].NodeName] == p2
                        for p3 in range(0, parents[3].nvalues):
                            d = df[parents[3].NodeName] == p3
                            for p4 in range(0, node.nvalues):
                                e = df[node.NodeName] == p4
                                counts.append(
                                    pd.DataFrame(df[a & b & c & d
                                                    & e]).shape[0] + 1)
                                v0.append(p0)
                                v1.append(p1)
                                v2.append(p2)
                                v3.append(p3)
                                v4.append(p4)
            cptDf = pd.DataFrame({
                node.NodeName: v4,
                parents[3].NodeName: v3,
                parents[2].NodeName: v2,
                parents[1].NodeName: v1,
                parents[0].NodeName: v0,
                'p': np.ones(len(counts)) * (-1),
                'counts': counts
            })
            node.setCPTData(cptDf)

    for X in bn.PresGraph.keys():
        bn.normalizeCpt(X)


def normaliseArray(vals):
    denom = np.sum(vals)
    normVals = []
    for val in vals:
        normVals.append(val / float(denom))
    return normVals


def getAssignmentFor(factor, E, nval):
    currFactor = factor
    for key, value in E.items():
        if key in list(factor.columns):
            condition = currFactor[key] == value
            currFactor = currFactor[condition]
        if currFactor.shape[0] == nval:
            return currFactor
    return currFactor


def markovBlanketsampling(X, E, bn):
    distX = []
    children = bn.PresGraph[X].getChildren()
    parents = bn.PresGraph[X].getParents()
    XCpt = bn.PresGraph[X].CPTData
    facX = getAssignmentFor(XCpt, E, bn.PresGraph[X].nvalues)
    facC = np.log(np.asarray(facX['p']))

    for child in children:
        cCPT = bn.PresGraph[list(bn.PresGraph.keys())[child]].CPTData
        temp = getAssignmentFor(cCPT, E, bn.PresGraph[X].nvalues)
        facC += np.log(np.asarray(temp['p']))
    return normaliseArray(np.exp(facC))


def Expectation(bn, df, missingIndex):
    """
    Return
        newdf: each missing value in a row replaced by the possible values it can take
        newWeight: array of weigh assigned for each data point
    """
    newWeights = []
    mydict = df.to_dict(orient='records')
    newDf = pd.DataFrame()
    for i in range(0, df.shape[0]):
        row = pd.DataFrame(df.loc[i, ]).T
        X = list(bn.PresGraph.keys())[missingIndex[i]]
        mbX = bn.MarkovBlanket[X]
        E = {
            key: value
            for key, value in mydict[i].items() if (key != X and key in mbX)
        }
        distX = markovBlanketsampling(X, E, bn)
        for n in range(bn.PresGraph[X].nvalues):
            row.iloc[0, list(bn.PresGraph.keys()).index(X)] = n
            newWeights.append(distX[n])
            newDf = pd.concat([newDf, row])
    return newWeights, newDf


def Maximization(df, net, weights):
    df['wts'] = weights
    N = df.shape(0)
    currentIter = 0
    for X, node in net.PresGraph.items():
        parents = net.getParents(node)
        nParents = len(parents)
        if nParents == 0:
            counts = []
            for p0 in range(0, node.nvalues):
                a = df[node.NodeName] == p0
                count = float(pd.DataFrame(df[a])['wts'].sum())
                if count != 0:
                    counts.append(float(count))
                else:
                    counts.append(0.000005)
            node.CPTData['counts'] = counts
            net.normalizeCpt(X)

        elif nParents == 1:
            counts = []
            for p0 in range(0, parents[0].nvalues):
                a = df[parents[0].NodeName] == p0
                for p1 in range(0, node.nvalues):
                    b = df[node.NodeName] == p1
                    count = float(pd.DataFrame(df[a & b])['wts'].sum())
                    if count != 0:
                        counts.append(float(count))
                    else:
                        counts.append(0.000005)
            node.cpt_data['counts'] = counts
            net.normalise_cpt(X)

        elif nParents == 2:
            counts = []
            for p0 in range(0, parents[0].nvalues):
                a = df[parents[0].NodeName] == p0
                for p1 in range(0, parents[1].nvalues):
                    b = df[parents[1].NodeName] == p1
                    for p2 in range(0, node.nvalues):
                        c = df[node.NodeName] == p2
                        count = float(pd.DataFrame(df[a & b & c])['wts'].sum())
                        if count != 0:
                            counts.append(count)
                        else:
                            counts.append(0.000005)
            node.cpt_data['counts'] = counts
            net.normalise_cpt(X)

        elif nParents == 3:
            counts = []
            for p0 in range(0, parents[0].nvalues):
                a = df[parents[0].NodeName] == p0
                for p1 in range(0, parents[1].nvalues):
                    b = df[parents[1].NodeName] == p1
                    for p2 in range(0, parents[2].nvalues):
                        c = df[parents[2].NodeName] == p2
                        for p3 in range(0, node.nvalues):
                            d = df[node.NodeName] == p3
                            count = float(
                                pd.DataFrame(df[a & b & c & d])['wts'].sum())
                            if count != 0:
                                counts.append(float(count))
                            else:
                                counts.append(0.000005)
            node.cpt_data['counts'] = counts
            net.normalise_cpt(X)

        elif nParents == 4:
            counts = []
            for p0 in range(0, parents[0].nvalues):
                a = df[parents[0].NodeName] == p0
                for p1 in range(0, parents[1].nvalues):
                    b = df[parents[1].NodeName] == p1
                    for p2 in range(0, parents[2].nvalues):
                        c = df[parents[2].NodeName] == p2
                        for p3 in range(0, parents[3].nvalues):
                            d = df[parents[3].NodeName] == p3
                            for p4 in range(0, node.nvalues):
                                e = df[node.NodeName] == p4
                                count = float(
                                    pd.DataFrame(df[a & b & c & d
                                                    & e])['wts'].sum())
                                if count != 0:
                                    counts.append(float(count))
                                else:
                                    counts.append(0.000005)
            node.cpt_data['counts'] = counts
            net.normalise_cpt(X)

        currentIter += 1


def ExpMax(df, bn, missingIndex):
    currentIter = 0
    time_i = time.time()
    while True:
        print("ITERATION.......", currentIter)
        step0 = time.time()
        wts, newDf = Expectation(bn, df, missingIndex)
        prevCpts = []
        for X in bn.PresGraph.keys():
            prevCpts.append(np.array(list(bn.PresGraph[X].CPTData['p'])))
        Maximization(newDf, bn, wts)
        newCpts = []
        for X in bn.PresGraph.keys():
            newCpts.append(np.array(list(bn.PresGraph[X].CPTData['p'])))
        diffs = []
        for i in range(len(prevCpts)):
            maxdiff = max(abs(np.subtract(prevCpts[i], newCpts[i])))
            diffs.append(maxdiff)
        delta = max(diffs)
        print("Delta.....", delta)

        if delta <= 0.00005:
            break
        currentIter += 1
    print("Convergence reached in .......", currentIter)
    return bn

