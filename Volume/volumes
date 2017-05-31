import subprocess
import scipy.special
import numpy
#import glpk
from sympy import *
import sys



def ehrhart_from_dine(filename, t, method="lp", debug=True):

    # Read the dimension from the .ine
    d=0
    with open(filename, "r") as myfile:
        lines = myfile.readlines()
        for i in range(len(lines)):
            if "begin\n" in lines[i].split(" "):
                # The dimension is the (second number)-1
                #  in the line after begin
                ll=lines[i+1].split(" ")
                while "" in ll:
                    ll.remove("")  
                d=int(ll[1])-1
                break
    if debug: print "Dimension: ", d
    
   # Read the .ine file and write a .dine (dilated by t) 
    with open(filename, "r") as myfile:
        lines = myfile.readlines()
        j=0
        with open(filename+".dine", "w") as dine:
            jump = 0 
            for i in range(len(lines)):
                print "in line ", i
                
                if "begin\n" in lines[i].split(" "):
                    dine.write(lines[i])
                    dine.write(lines[i+1])
                    ll=lines[i+1].split(" ")
                    ll.remove("")
                    jump = int(ll[0])
                    for j in range(2,jump+2):
                        ll = lines[i+j].split(" ")
                        for st in ll:
                            st.replace(" ", "")
                        while "" in ll:
                            ll.remove("")    
                        if "\n" in ll:    
                            ll.remove("\n")            
                        dine.write( str(int(ll[0])*(t)) + " " + " ".join(ll[1:])+ "\n")
                    jump= i + jump + 1
                else:
                    if i>jump:
                        dine.write(lines[i])
    nexps=10
    output = subprocess.check_output(["./vol", "--exact", "-exp" , str(nexps), "-f1", filename+".dine"])
    if debug: print "Volume approximation output:\n" + output + "\n"
    #print output.split("[")[-1].split("]")[0].split(",")
    L = [float(i) for i in output.split("[")[-1].split("]")[0].split(",")]
    
    slack = 0
    #if debug: print "Bounds: ",L[0],L[1]
    #if debug: print "Dimension: ", d, t
    L = [int(round((t^(d))*i)) for i in L]
    L = [ L[0]-slack , L[1]+slack ]
 
    if debug: print "Bounds: ", L

    ie = [ [ L[0] ] + [ int(scipy.special.binom(t+d-i,d)) for i in range(d+1)], \
            [ -L[1] ] + [ -i for i in [ int(scipy.special.binom(t+d-i,d)) for i in range(d+1)]] ] \
        + [ [ 1 ] + [ 1 if j==0 else 0 for j in range(d+1)] ]\
        + [ [ -1 ] + [ -1 if j==0 else 0 for j in range(d+1)] ]\
        + [ [ 0] + [ 1 if i==j else 0 for j in range(d+1)] for i in range(1,d+1)]

    if debug: 
        print "System: "
        for i in ie:
            print i
    
    # To compute only one lattice point, optimizing a linear function
    if method=="lp":
        lpfilename = filename + ".lp"
        with open(lpfilename, "w") as f:
    #        f.write("maximize\n")
            f.write("minimize\n")
            f.write("+".join([ "x"+str(i) for i in range(2,len(ie[0]))])+"\n")
            f.write("subject to\n")
            for ieq in ie:
                c = ieq[0]
                ieq = ieq[1:]
                s = "+".join([ str(int(ieq[i]))+"x"+str(i+1)  for i in range(len(ieq)) if ieq[i]!=0])
                s = s.replace("+-", "-")+ "> " +str(int(c)) + "\n"
                f.write(s)			
            f.write("bounds\n")	
            f.write("\n".join([ "x"+str(i)+">0" for i in range(1,len(ie[0]))])+"\n")
            f.write("integer\n")	    
            f.write("\n".join([ "x"+str(i) for i in range(1,len(ie[0]))])+"\n")
            f.write("end\n")	    

        if debug: print "LP file created"
        subprocess.check_output(["glpsol", "--cpxlp", lpfilename, "-o", lpfilename +".out"])
        if debug: print "LP call completed"

        lpoutput=""
        with open(lpfilename +".out", "r") as myfile:
            lpoutput=myfile.read().replace('\n', '')
        L = lpoutput.split(' ')
        while "" in L:
            L.remove("")

        h=[0 for i in range(d+1)]
        while '*' in L:
            p = L.index('*')
            # zero index, glpk is 1-indexed
            v = int(L[p-1].replace("x",""))-1
            h[v] = int(float(L[p+1]))
            L.remove(L[p])
        if debug: print h
    # Compute all lattice points in the interpolation polytope
    else:
        with open(filename+"_new.ine", "w") as f:
            f.write("H-representation\n")
            f.write("begin\n")
            f.write( str(len(ie))+ " " + str(len(ie[0]))+" " + "real\n")
            
            for ineq in ie:
            	f.write(" ".join(str((-1)*ineq[i]) if i==0 else str(ineq[i]) for i in range(len(ineq)))+ "\n")
            f.write("end\ninput_incidence")
            f.close()
        output = subprocess.check_output(["/usr/local/bin/ConvertCDDineToLatte", filename+"_new.ine", filename+".new.latte"])
#            f.write( str(len(ie))+ " " + str(len(ie[0]))+ "\n")    
#            for ieq in ie:
#                f.write(" ".join([ str(i) for i in [ieq[0]]+ [(-1)*j for j in ieq[1:]]])+ "\n")			
        
        output = subprocess.check_output(["/usr/local/bin/count", "--multivariate-generating-function", filename+".new.latte"])
	numofpoints=0
        with open("numOfLatticePoints", "r") as f:
		tmp= "".join(f.readlines())
		numofpoints=int(tmp)
		print "------------------>", numofpoints
        output = [ s for s in output.split("\n") ]
        lattemultipoly=""
        with open(filename+".new.latte.rat", "r") as myfile:
            symbols("".join(["t"+str(i) for i in range(d)]))
            lattemultipoly =  "".join(myfile.readlines()).replace("\n","")
        for i in range(d+1):
            lattemultipoly =lattemultipoly.replace("x["+ str(i) +"]", "t"+str(i))

        print "ok"
        print "lattemultipoly: ",lattemultipoly

        stexp=str(expand(simplify(expand(simplify(expand(simplify(lattemultipoly)))))))
        #stexp = str(expand(lattemultipoly))
        #print "stepx: ",stexp
        #stexp = str(simplify(stexp))
        while "**" in stexp:
            stexp=stexp.replace("**","^")
        if debug: print "The multivariate gen fun", stexp 
        terms = [ term.split("*")  for term in stexp.split("+")]
        points=[]
        for cpoint in terms:
            point=[0 for i in range(d+1)]
            if "t" in cpoint[0]:
                for c in cpoint:
                    c = c.replace(" ","")
                    tmp = c.split("^")
                    if len(tmp)==1:    
                        point[int(tmp[0][1:])]=1
                    else:
                        point[int(tmp[0][1:])]=int(tmp[1])
            points.append(point)
        print points
        print "Ehrhart candidates"
        for i in range(len(points)):
            h= points[i]

            t = symbols("t")
            B = [ binomial(t+len(h)-1-i,len(h)-1) for i in range(len(h))]
            ehrhart = expand(combsimp(sum([ h[i]*B[i] for i in range(len(h))])))
            ehrhart = factor(ehrhart)
            print ehrhart

    t = symbols("t")
    B = [ binomial(t+len(h)-1-i,len(h)-1) for i in range(len(h))]
    ehrhart = expand(combsimp(sum([ h[i]*B[i] for i in range(len(h))])))
    ehrhart = factor(ehrhart)
    print ehrhart

    if debug: print "Ehrhart polynomial found:", ehrhart
    output = subprocess.check_output(["/usr/local/bin/ConvertCDDineToLatte", filename, filename+".latte"])

    output = subprocess.check_output(["/usr/local/bin/count", "--ehrhart-polynomial", filename+".latte"])
    output = [ s for s in output.split("\n") ]
    lattepoly =  expand(output[-3])
    print  "The polynomial of LattE:" , lattepoly
    print  "The polynomial of Volumes:" , expand(ehrhart)
    agree = (ehrhart == lattepoly)
    print "Correct: ", agree
    return agree

filename = sys.argv[1] 
t=int(sys.argv[2])
ehrhart_from_dine(filename,t,method="lp")
#ehrhart_from_dine("cube_3.ine",t)
#ehrhart_from_dine("cube10.ine",t)
