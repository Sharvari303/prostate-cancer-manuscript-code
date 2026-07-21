#!/usr/bin/env python

import re

#---first if only one otherwise None
uniq = lambda x : x[0] if len(x)==1 else None
#---first element or none
forn = lambda obj : ([obj][:1] or [None])[0]

def asciitree(obj,depth=0,wide=2,last=[],recursed=False):

    """
    Print a dictionary as a tree to the terminal.
    Includes some simuluxe-specific quirks.
    """

    corner = u'\u251C'
    horizo = u'\u2500'
    vertic = u'\u2502'

    spacer = {0:'\n',
        1:' '*(wide+1)*(depth-1)+corner+horizo*wide,
        2:' '*(wide+1)*(depth-1)
        }[depth] if depth <= 1 else (
        ''.join([(vertic if d not in last else ' ')+' '*wide for d in range(1,depth)])
        )+corner+horizo*wide
    if type(obj) in [str,float,int,bool]:
        if depth == 0: print spacer+str(obj)+'\n'+horizo*len(obj)
        else: print spacer+str(obj)
    elif type(obj) == dict and all([type(i) in [str,float,int,bool] for i in obj.values()]) and depth==0:
        asciitree({'HASH':obj},depth=1,recursed=True)
    elif type(obj) == list:
        for ind,item in enumerate(obj):
            if type(item) in [str,float,int,bool]: print spacer+str(item)
            elif item != {}:
                print spacer+'('+str(ind)+')'
                asciitree(item,depth=depth+1,
                    last=last+([depth] if ind==len(obj)-1 else []),
                    recursed=True)
            else: print 'unhandled tree object'
    elif type(obj) == dict and obj != {}:
        for ind,key in enumerate(obj.keys()):
            if type(obj[key]) in [str,float,int,bool]: print spacer+key+' = '+str(obj[key])
            #---special: print single-item lists of strings on the same line as the key
            elif type(obj[key])==list and len(obj[key])==1 and type(obj[key][0]) in [str,float,int,bool]:
                print spacer+key+' = '+str(obj[key])
            #---special: skip lists if blank dictionaries
            elif type(obj[key])==list and all([i=={} for i in obj[key]]):
                print spacer+key+' = (empty)'
            elif obj[key] != {}:
                #---fancy border for top level
                if depth == 0:
                    print '\n'+corner+horizo*(len(key)+0)+corner+spacer+vertic+str(key)+vertic+'\n'+\
                        corner+horizo*len(key)+corner+'\n'+vertic
                else: print spacer+key
                asciitree(obj[key],depth=depth+1,
                    last=last+([depth] if ind==len(obj)-1 else []),
                    recursed=True)
            elif type(obj[key])==list and obj[key]==[]:
                print spacer+'(empty)'
            else: print 'unhandled tree object'
    else: print 'unhandled tree object'
    if not recursed: print '\n'

#---remove tags
deltag = lambda r : re.findall('^{[^}]+}(.+)',r.tag).pop()

def xml2tree(obj):

    """
    Return tag tree for an XML element tree object.
    """

    kids = obj.getchildren()
    if kids == []: return deltag(obj)
    else:
        out = {}
        for kid in kids: out[deltag(kid)] = xml2tree(kid)
        return out

def text2array(text): return array(text.strip('\n').split()).astype(float)

def parse_modification(modtext):

    """
    Convert a human-readable modification script into a dictionary.
    """

    entry_types = ['change','remove']
    entry_processor = {
            'change':lambda x : {'path':tuple(x[1].split()),'key':x[2].split()[0],'value':x[2].split()[1]},
            }

    lines = modtext.strip('\n').split('\n')
    line_nos = [ii for ii,i in enumerate(lines) if re.match('^(%s)\s*$'%'|'.join(entry_types),i)]
    #---only valid record breaker is a single blank line
    break_nos = [ii for ii,i in enumerate(lines) if re.match('^$',i)]+[len(lines)]
    if len(line_nos)!=len(break_nos): raise Exception('parse error')
    entries = [lines[line_nos[a]:break_nos[a]] for a in range(len(line_nos))]
    mods = [entry_processor[entry[0]](entry) for entry in entries]
    return mods

def call(command,logfile=None,cwd=None,silent=False,inpipe=None,suppress_stdout=False):

    """
    Wrapper for system calls in a different directory with a dedicated log file.
    """

    if inpipe != None:
        if logfile == None: output=None
        else: output = open(('' if cwd == None else cwd)+logfile,'wb')
        if type(command) == list: command = ' '.join(command)
        p = subprocess.Popen(command,stdout=output,stdin=subprocess.PIPE,stderr=output,cwd=cwd,shell=True)
        catch = p.communicate(input=inpipe)[0]
    else:
        if type(command) == list: command = ' '.join(command)
        if logfile != None:
            output = open(('' if cwd == None else cwd)+logfile,'wb')
            if type(command) == list: command = ' '.join(command)
            if not silent: print 'executing command: "'+str(command)+'" logfile = '+logfile
            try:
                subprocess.check_call(command,
                        shell=True,
                        stdout=output,
                        stderr=output,
                        cwd=cwd)
            except:
                if logfile[-3:] == '-cg' and re.search('mdrun-em',logfile):
                    print 'warning: failed conjugate gradient descent but will procede'
                else: raise Exception('except: BASH execution error (probably GROMACS)!\nsee '+cwd+logfile)
            output.close()
        else:
            if not silent: print 'executing command: "'+str(command)+'"'
            if str(sys.stdout.__class__) == "<class 'amx.tools.tee'>": stderr = sys.stdout.files[0]
            else: stderr = sys.stdout
            try:
                if suppress_stdout:
                    devnull = open('/dev/null','w')
                    subprocess.check_call(command,shell=True,stderr=devnull,cwd=cwd,stdout=devnull)
                else: subprocess.check_call(command,shell=True,stderr=None,cwd=cwd)
            except:
                raise Exception('except: BASH execution error\ncommand: '+\
                        command+'\ncwd: '+cwd)
