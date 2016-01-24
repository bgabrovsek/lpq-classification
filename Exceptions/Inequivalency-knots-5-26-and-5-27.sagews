︠56ad233b-1060-42a8-8722-0a5b7a3b42b9s︠
%python

from snappy import *

lens = [(1,2),(1,3),(1,4),(1,5),(2,5),(1,6),(1,7),(2,7),(1,8),(3,8),(1,9),(2,9),(1,10),(3,10),(1,11),(2,11),(3,11),(1,12),(5,12)]

def fill_component(m, pq, comp = 1):
    m_ = m.copy()
    m_.dehn_fill(pq,comp)
  
    return m_.filled_triangulation()
    
def non_inv(m, bp):
    if (m.verify_hyperbolicity(bits_prec=bp)[0] == True):
        t = m.canonical_retriangulation(verified = True)
        l = len(t.isomorphisms_to(t))
        #s = m.symmetry_group()
        if l == 1:
            print "knot is non-invertible (has trivial symmetry group)\n"
        else:
            print "symmetries found in canonical triangulation\n"
    else:
        print "Could not determine hyperbolicity.\n"
        
print "Checking invertibility for knot 5_26.\n"

L = Link("link-5-26.lnk")
M = L.exterior()

print "Solid Torus"
non_inv(M,)

for (i,pq) in enumerate(lens):
    print "Lens space", pq, ". Filling cusp with", pq, "Dehn filling."
    N = fill_component(M,pq)
    non_inv(N,100)

︡5ae8a2cf-9bc4-4198-b486-5c2fe5a20dbd︡︡{"stdout":"Checking invertibility for knot 5_26.\n\n","done":false}︡{"stdout":"Solid Torus\n","done":false}︡{"done":false,"stderr":"Error in lines 22-22\nTraceback (most recent call last):\n  File \"/projects/sage/sage-6.10/local/lib/python2.7/site-packages/smc_sagews/sage_server.py\", line 905, in execute\n    exec compile(block+'\\n', '', 'single') in namespace, locals\n  File \"\", line 1, in <module>\nTypeError: non_inv() takes exactly 2 arguments (1 given)\n"}︡{"done":true}
︠009f00cf-5e8b-491b-a030-c9256dbde86f︠









