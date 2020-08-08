from collections import OrderedDict
import pandas as pd
import numpy as np
# from pandas import dataFrame


class GraphNode():
    def __init__(self,name,n,vals):
        self.NodeName=name#variable Name
        self.nvalues=n#No. of categories this variable can take
        self.values=vals#possible values
        
        self.Children=[]#To store the children
        self.Parents=[]#to store the Parents
        self.CPT=[]#Coditional Probibilty table
        self.markovBlanket=[]#to store the markov blanket
        self.CPTData=pd.DataFrame()

    
    def getName(self):
        return self.NodeName
    
    def getChildren(self):
        return self.Children
    
    def getParents(self):
        return self.Parents
    
    def getCPT(self):
        return self.CPT
    
    def getNValue(self):
        return self.nvalues 

    def getValues(self):
        return self.values 

    def getMarkovBlanket(self):
        return self.markovBlanket

    def setCPT(self, new_CPT):
        del(self.CPT)
        self.CPT=new_CPT
    
    def setCPTData(self,newCPTData):
        self.CPTData.drop(columns=list(self.CPTData.columns))
        self.CPTData=newCPTData
    
    def setParents(self,new_Parents):
        self.Parents=new_Parents
    
    def setMB(self,MB):
        self.markovBlanket=MB

    def addChild(self, newChildIndex):
        if newChildIndex in self.Children:
                return 0
        else:
            self.Children.append(newChildIndex)
            return 1

    def printNode(self):
        print(self.NodeName)
        print(self.values)
        print(self.nvalues)
        print(self.CPT)
        print(self.Parents)
        print(self.Children)

class Network():
    """
    The whole graph is a Dict of Node. 
    PresGraph= key:Node name/Variable Name  Value:Node
    MB=key:Node name/Variable Name  Value:List of Nodes that are the MarkovBlanket of the Node
    """

    def __init__(self,PresGraph=OrderedDict(),MarkovBlanket=OrderedDict()):
        self.PresGraph=PresGraph
        self.MarkovBlanket=MarkovBlanket
    
    def addNode(self,node):
        self.PresGraph[node.getName()]=node

    def netSize(self):
        return len(self.PresGraph)
    
    def getindex(self,varName):
        if varName in self.PresGraph:
            return list(self.PresGraph.keys()).index(varName)
        else:
            return -1

    def getNthNode(self,n):
        return list(self.PresGraph.values())[n]

    def searchNode(self,varName):
        try:
            return self.PresGraph[varName]
        except:
            return None
    
    def getParents(self,node):
        parentNodes=[]
        parents=node.getParents()
        for p in parents:
            parentNodes.append(self.searchNode(p))
        return parentNodes

    def getChildren(self,valName):
        children=[]
        Children=self.PresGraph[valName].getChildren()
        for child in Children:
            children.append(list(self.PresGraph.keys())[child])
        return children

    
    def setMb(self):
        for varName in self.PresGraph.keys():
            self.MarkovBlanket[varName]=getMarkovBlanket(self,varName)

    def normalizeCpt(self,X):
        l=[X]+self.PresGraph[X].getParents()+['p','counts']
        cpt=self.PresGraph[X].CPTData
        nvals=self.PresGraph[X].getNValue()
        cardinality=cpt.shape[0]
        noOfGroups=int(cardinality/nvals)
        df=pd.DataFrame()
        i=0
        for n in range(noOfGroups):
            currDf=pd.DataFrame(cpt.iloc[i:i+nvals,:])
            currDf['p']=normalizeCounts(currDf['counts'])
            df=df.append(currDf)
            i+=nvals
        self.PresGraph[X].CPTData=df[l]



def getMarkovBlanket(network,varName):
    mb=[]
    mb.append(varName)
    Node=network.searchNode(varName)
    for parent in Node.getParents():
        mb.append(parent)
    for childNumber in Node.getChildren():
        child=network.getNthNode(childNumber)
        # childName=childGetter.keys()
        mb.append(child.getName())
        # now for parent of Child
        for parentofchild in child.getParents():
            if parentofchild not in mb:
                mb.append(parentofchild)
    
    return mb


def readNetwork(bifFile):
    Alarm=Network()
    with open(bifFile,"r") as file:
        while True:
            line=file.readline().strip()
            if line=="":
                break
            token=line.split()

            if token[0]=="variable":
                name=token[1]
                nextline=file.readline().strip()
                nextToken=nextline.split(" ")
                # size=str(nextToken[2])[9]
                values=[]
                for i in range(4,len(nextToken)-1,2):
                    values.append(nextToken[i])
                newNode=GraphNode(name=name,n=len(values),vals=values)
                Alarm.addNode(newNode)
                
            if token[0]=="probability":
                vals=[]
                tmpName=token[2]
                index=Alarm.getindex(tmpName)
                node=Alarm.searchNode(tmpName)
                i=3
                # print(line,tmpName)
                while token[i]!=")":
                    newNode=Alarm.searchNode(token[i])
                    newNode.addChild(index)
                    vals.append(token[i])
                    i+=1
                node.setParents(vals)
                cptline=file.readline().strip()
                cpttoken=cptline.split(" ")
                CPT=[]
                i=1
                while cpttoken[i]!=";":
                    CPT.append(cpttoken[i])
                    i+=1
                node.setCPT(CPT)
    file.close()
    return Alarm


def printNetwork(network):
    PresGraph=network.PresGraph
    MarkovBlanket=network.MarkovBlanket
    for Name,Node,mb in zip(PresGraph.keys(),PresGraph.values(),MarkovBlanket.values()):
        print(Name,network.getindex(Name))
        print(mb)
        print(Node.printNode())
        print(Node.CPTData)
        print("-----------")

def getData(records):
    with open(records,'r') as f:
        df=pd.DataFrame(l.rstrip().split() for l in f)
    df.columns = ['"Hypovolemia"','"StrokeVolume"','"LVFailure"','"LVEDVolume"','"PCWP"','"CVP"','"History"','"MinVolSet"','"VentMach"','"Disconnect"','"VentTube"','"KinkedTube"','"Press"','"ErrLowOutput"','"HRBP"','"ErrCauter"','"HREKG"','"HRSat"','"BP"','"CO"','"HR"','"TPR"','"Anaphylaxis"','"InsuffAnesth"','"PAP"','"PulmEmbolus"','"FiO2"','"Catechol"','"SaO2"','"Shunt"','"PVSat"','"MinVol"','"ExpCO2"','"ArtCO2"','"VentAlv"','"VentLung"','"Intubation"']
    features=list(df.columns)
    # float('nan')
    mapping1={'True':0,'False':1,'?':-71}
    mapping2={'Zero':0,'Low':1,'Normal':2,'High':3,'?':-71}
    mapping3={'Normal':0,'Esophaegeal':1,'OneSided':2,'?':-71}
    mapping4={'Low':0,'Normal':1,'High':2,'?':-71}
    mapping5={'Low':0,'Normal':1,'High':2,'?':-71}
    mapping6={'Normal':0,'High':1,'?':-71}
    OverallMapping={ '"Hypovolemia"':mapping1 , u'"StrokeVolume"':mapping4, u'"LVFailure"':mapping1,u'"LVEDVolume"':mapping4, u'"PCWP"':mapping4, u'"CVP"':mapping4,  u'"History"':mapping1, u'"MinVolSet"':mapping4, u'"VentMach"':mapping2, u'"Disconnect"':mapping1,u'"VentTube"':mapping2, u'"KinkedTube"':mapping1, u'"Press"':mapping2,u'"ErrLowOutput"':mapping1, u'"HRBP"':mapping4,u'"ErrCauter"':mapping1, u'"HREKG"':mapping4, u'"HRSat"':mapping4,u'"BP"':mapping4, u'"CO"':mapping4, u'"HR"':mapping4, u'"TPR"':mapping4,u'"Anaphylaxis"':mapping1, u'"InsuffAnesth"':mapping1, u'"PAP"':mapping4,u'"PulmEmbolus"':mapping1, u'"FiO2"':mapping5,u'"Catechol"':mapping6, u'"SaO2"':mapping4, u'"Shunt"':mapping6,u'"PVSat"':mapping4, u'"MinVol"':mapping2, u'"ExpCO2"':mapping2,u'"ArtCO2"':mapping4, u'"VentAlv"':mapping2, u'"VentLung"':mapping2, u'"Intubation"':mapping3}
    df=df.replace(OverallMapping)
    return df

def normalizeCounts(vals):
    vals[vals==0]=0.000005
    total=np.sum(vals)
    normalizedVals=[]
    for val in vals:
        normalizedVals.append(val/total)
    return normalizedVals
