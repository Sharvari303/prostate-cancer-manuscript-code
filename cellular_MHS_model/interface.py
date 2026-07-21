#!/usr/bin/env python

import os,json,re
import xml.etree.ElementTree as ET
#ET.register_namespace("",'{http://www.copasi.org/static/schema}')
from tools import asciitree, xml2tree, uniq, forn

#---CLASSES
#-------------------------------------------------------------------------------------------------------------

class CopasiState():

    """
    Class which loads, stores, manipulates, and writes COPASI CPS file representing the state of a COPASI
    simulation.
    """

    #---class constants
    namespace_tag = ['{http://www.copasi.org/static/schema}',''][1]

    #---lookup rules
    rules = [
            lambda elem,tag : forn(uniq(elem.findall(CopasiState.namespace_tag+tag))),
            lambda elem,tag : forn(uniq([a for a in elem if a.attrib['name']==tag])),
            lambda elem,tag : forn(uniq([a for a in elem if a.attrib['key']==tag])),
            lambda elem,tag : forn(uniq([a for a in elem if a.attrib['cn']==tag])),
            lambda elem,tag : forn(uniq([a for a in elem
                    if dict([r.split('=') for r in a.attrib['cn'].split(',')])['String']==tag])),
            #---very open ended rule, beware!
            lambda elem,tag : None if '@' not in tag else forn(uniq(
                    [i for i in elem if
                    set([l for l in [uniq(k) for k in [re.findall('\[([^\]]+)\]',j.split('=')[1])
                    for j in i.attrib['cn'].split(',')]] if l])==set(tag.split('@'))])),
            #---! highly specific rulemay need rewritten
            lambda elem,tag : forn(uniq([a for a in elem
                    if re.findall('Values\[([^\]]+)',dict([r.split('=')
                    for r in a.attrib['cn'].split(',')])['Vector']).pop()==tag])),
            lambda elem,tag : forn(uniq([a for a in elem if 'target' in a.attrib and a.attrib['target']==tag]))
            ]

    #---operate without namespaces
    def noname(self,tag): return re.findall('^\{[^\}]+\}(.+)',tag).pop()

    def __init__(self,init_cps,script=None):

        #---constants
        self.tag_prefix = '{http://www.copasi.org/static/schema}'
        #---variables
        #---! multiple ways to load the state or should it always be from a cps file?
        self.init_cps = init_cps
        self.tree = None
        self.root = None
        #---load the initial cps file
        self.read(init_cps)
        self.space = {}
        if script != None: self.source_script(script)

    def __str__(self):

        asciitree(xml2tree(self.root))
        return 'class is CopasiSimulation (see above for tree)'

    def source_script(self,script):

        #---script interpreter rules
        interp_rules = {
                'path':lambda x : dict([(key,{'key':val,'value':tuple(val.split(','))}[key])
                        for key,val in dict([(j,re.findall('^(\w+)\s*=\s*(.+)',x).pop()[jj])
                        for jj,j in enumerate(['key','value'])]).items()]),
                }

        #---interpret the script
        self.space,self.spacetypes = {},{}
        for line in [l for l in script.split('\n') if re.match('^(\w+)\s',l)]:
            name,code = re.findall('^(\w+)\s+(.+)',line).pop()
            convert = interp_rules[name](code)
            self.space[convert['key']] = convert['value']
            self.spacetypes[convert['key']] = name

    def read(self,fn):

        """
        Read an XML state from a CPS file.
        """

        with open(fn,'r') as fp: xmlstring = fp.read()
        xmlstring = re.sub(' xmlns="[^"]+"','',xmlstring,count=1)
        self.root = ET.fromstring(xmlstring)
        #---Newer COPASI (>=4.3x) writes MiriamAnnotation metadata elements
        #with no attributes throughout the tree (e.g. as a direct child of
        #ModelParameterSet). The attribute-based lookup rules in stepdown()
        #raise on these attribute-less elements and abort the traversal.
        #These elements are pure RDF metadata and are not used by the
        #simulation, so strip them after parsing to keep traversal robust.
        for parent in self.root.iter():
            for child in list(parent):
                if child.tag.split('}')[-1] == 'MiriamAnnotation':
                    parent.remove(child)

    def stepdown(self,elem,tag):

        """
        Move down a tree by one step according to a set of rules.
        """

        down,next = elem,None
        for rr,rule in enumerate(CopasiState.rules):
            if next == None:
                #---? is try neccessary?
                try: next = rule(down,tag)
                except: pass
                if next != None: return next
        #---! very useful to explain the lookup failures here
        raise Exception('[ERROR] lookup failure: %s,%s'%(str(elem),tag))

    def _traverse(self,*path,**kwargs):

        """
        Traverses the ElementTree according to a path and returns the root object or list of children.
        """

        #---start our traversal at the root node
        elem = self.root
        #---for each string in the the *path list, step down to the next level of the tree
        for p in path: elem = self.stepdown(elem,p)
        #---do nothing if kwargs lacks a key,value pair
        if not ('key' in kwargs and 'value' in kwargs): pass
        #---if key,value are in kwargs, then change the attributes dictionary
        else: elem.attrib[kwargs['key']] = kwargs['value']
        return elem

    def _mod(self,*path,**kwargs):

        """
        Modify an XML element of the state.
        """

        for key,value in kwargs.items():
            self._traverse(*path).set(str(key),str(value))

    def change(self,*path,**kwargs):

        """
        """

        explicit_path = (self.space[path[0]]+tuple(path[1:])) if path[0] in self.space else path
        self._mod(*explicit_path,**kwargs)

    def get(self,*path):

        """
        A wrapper for _traverse which checks the self.space (a type of "namespace")
        for a predefined path. If it finds a predfined path in self.space, it prepends that path to the
        incoming *path list.
        """

        if path[0] in self.space: return self._traverse(*(self.space[path[0]]+tuple(path[1:])))
        else: return self._traverse(*path)

    def write(self,outfn,override=False):

        """
        Write the current state of the COPASI simulation to disk.
        """

        if not override and os.path.isfile(outfn): raise Exception('[ERROR] file %s already exists'%outfn)
        else:
            ET.ElementTree(self.root).write(outfn,
                    encoding='utf-8')

    def catalog(self,root,*names,**kwargs):

        """
        Given the name of an element of comma-delimited 'cn' keys, extract all possible values from
        `Vector=NAME[VALUE]` substrings.
        """

        #---collect unique combinations of values for requested names
        debracket = lambda x : re.findall('\[([^\]]+)\]',x)
        if 'unique' in kwargs and kwargs['unique']:
            collected_values = []
            for elem in root:
                for j in elem.attrib['cn'].split(','):
                    extracts = dict([(name,uniq(debracket(uniq([j.split('=')[1]
                            for j in elem.attrib['cn'].split(',')
                            if re.match('^%s'%name,j.split('=')[1])]))))
                            for name in names])
                    collected_values.append(extracts)
        #---catalog values by name independently
        else:
            collected_values = dict([(i,[]) for i in names])
            for elem in root:
                for j in elem.attrib['cn'].split(','):
                    left,right = j.split('=')
                    inside = uniq(debracket(right))
                    for key in names:
                        if re.match('^%s'%key,right) and inside not in collected_values[key]:
                            collected_values[key].append(inside)
        return collected_values

    def catalog_DELETEME(self):

        """
        Traverse the tree and save the paths.
        UNDER DEVELOPMENTs
        """

        #---seed the incomplete list with the top level
        obj = self.root
        incs,paths = [],[]

        for i in obj.getchildren(): incs.append([self.noname(obj.tag),self.noname(i.tag)])
        while len(incs)>0:
            print incs
            inc = incs.pop()
            try: elem = self._traverse(*inc[1:])
            except: elem = []
            print elem
            for i in elem:
                if i.getchildren() == []:
                    new = inc+[self.noname(i.tag)]
                    if new not in paths: paths.append(new)
                else:
                    #---check if they are all the same
                    kids = i.getchildren()
                    if len(kids)!=1 and len(set([a.tag for a in kids]))==1:
                        #---? assumed at the end of the line but could do a is-one plunge?
                        for k in kids: paths.append(inc+[self.noname(i.tag),json.dumps(k.attrib)])
                    else:
                        #---if not all the same then we can just save their tags without redundancy
                        #---? where does this leave the possibilities for missing leaves?
                        for j in i.getchildren():
                            new = inc+[self.noname(i.tag)]
                            if new not in incs: incs.append(new)
        self.treepaths = paths

class CopasiSimulation():

    """
    Class which manages the states of a COPASI simulation.
    """

    def __init__(self):

        print "UNDER DEVELOPMENT"
