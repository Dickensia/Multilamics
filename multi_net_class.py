# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 10:06:42 2019

@author: Bioinfo
"""
import networkx as nx
#import pymnet as pn
class Multilayer:
    def __init__(self):
        self.orinet = nx.Graph()
        self.laylst = []
        self.laydict = {}
    def add_node(self,nodename,Layer=0,Prize=1,**params):
        self.orinet.add_node(str(nodename)+'∈L'+str(Layer),
                             layer=Layer,prize=Prize,**params)
        if Layer not in self.laylst:
            self.laylst.append(Layer)
            self.laydict['Layer'+str(Layer)]=[nodename]
        elif (Layer in self.laylst) and (nodename not in self.laydict['Layer'+str(Layer)]):
            self.laydict['Layer'+str(Layer)].append(nodename)
        self.laylst.sort()
    def add_node_from(self,nbunch,Layer=0,Prize=1,**params):
        '''nbunch is a list of nodenames attched with prizes'''
        for (node,prize) in nbunch:
            self.add_node(node,Layer,Prize=prize,**params)###此处参数无需给self，思考为什么？？
            #self.orinet.add_node(str(node)+'∈L'+str(Layer
            #                     ),layer=Layer,**params)
#{n for n, d in B.nodes(data=True) if d['bipartite']==0}
    def nodes(self,default=True,Layer=0):
        dict = self.orinet.nodes(data=True)._nodes
        #print(dict)
        if default == True:
            for layer in self.laylst:
                print('Layer '+str(layer)+': \n')
                for n in dict:
                    if dict[n]['layer']==layer:
                        print(n.split('∈')[0]+'\n')
        elif default == False:
            print('Layer '+str(Layer)+': \n')
            for n in dict:
                if dict[n]['layer'] == Layer:
                    print(n.split('∈')[0]+'\n')
        
    def add_edge(self,edgetuple,Weight=1):
        '''edgetuple is like ((node1,Layer1),(node2,Layer2))'''
        for (n,l) in edgetuple:#此处也可用try-except写，但是注意不能合并判断分支，会报错
            if (str('Layer'+str(l)) not in self.laydict.keys()):
                self.add_node(n,Layer=l)
            elif (n not in self.laydict['Layer'+str(l)]):
                self.add_node(n,Layer=l)
        self.orinet.add_edge((str(edgetuple[0][0])+'∈L'+str(edgetuple[0][1]
        ),str(edgetuple[1][0])+'∈L'+str(edgetuple[1][1])),weight=Weight)
        
    def add_edge_from(self,ebunch):
        '''ebunch is a list of edgestuple'''
        for edgetuple in ebunch:
            self.add_edge(edgetuple)
    def visualize(self):
        pass#后续用pymnet来实现
        

G = Multilayer()        
G.add_node('x',Layer=1,attr='x')
G.add_node_from([1,2,3,4],attr='y')
G.nodes(default=False)
G.add_edge(((1,0),('x',2)),2)
ebunch = [((2,0),('z',1)),(('u',1),('w',9))]
G.add_edge_from(ebunch)


