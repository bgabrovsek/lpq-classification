︠1ff7e4cc-6fcc-48e1-9ad5-669a4d0bd1cdr︠

%python

from snappy import *

lens = [(1,2),(1,3),(1,4),(1,5),(2,5),(1,6),(1,7),(2,7),(1,8),(3,8),(1,9),(2,9),(1,10),(3,10),(1,11),(2,11),(3,11),(1,12),(5,12)]

def fill_triangulation(m, pq, component = 1):
    m_ = m.copy()
    m_.dehn_fill(pq,component)
    return m_.filled_triangulation()

def are_equivalent(m, n, bp = 500):
    hm = m.verify_hyperbolicity(bits_prec=bp)[0]
    hn = n.verify_hyperbolicity(bits_prec=bp)[0]
    if (hm == True) and (hn == True):
        tm = m.canonical_retriangulation(verified = True)
        tn = n.canonical_retriangulation(verified = True)
        l = len(tm.isomorphisms_to(tn))
        if l == 0:
            print "Knots are not equivalent (0 isomorphisms between their canonical triangulations).\n"
        else:
            print l, "isomorphisms between their canonical triangulations found.\n"
    else:
        print "Could not determine hyperbolicity.\n"

a = Link("link-5-76.lnk")
b = Link("link-5-76-mirror.lnk")

print "Checking knots 5-76 and mirror 5-76.\n"

ma = a.exterior()
mb = b.exterior()

ha = ma.verify_hyperbolicity()[0]
hb = mb.verify_hyperbolicity()[0]

if (ha == False) or (hb == False):
    print "Could not determine hyperbolicity."

ta = ma.canonical_retriangulation(verified = True)
tb = mb.canonical_retriangulation(verified = True)

print "Solid Torus"
are_equivalent(ta,tb)

for pq in lens:
    print "Lens space", pq, ". Filling cusp with", pq, "Dehn filling."
    fta = fill_triangulation(ma, pq)
    ftb = fill_triangulation(mb, pq)
    are_equivalent(fta,ftb)
    
    


︡6da42b94-2a40-4916-adfc-dca33e814562︡{"stdout":"Checking knots 4-23 and mirror 5-1.\n\n","done":false}︡{"stdout":"Solid Torus\n","done":false}︡{"stdout":"4","done":false}︡{"stdout":" isomorphisms between their canonical triangulations found.\n\n","done":false}︡{"stdout":"Lens space (1, 2) . Filling cusp with (1, 2) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 3) . Filling cusp with (1, 3) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 4) . Filling cusp with (1, 4) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 5) . Filling cusp with (1, 5) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (2, 5) . Filling cusp with (2, 5) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 6) . Filling cusp with (1, 6) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 7) . Filling cusp with (1, 7) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (2, 7) . Filling cusp with (2, 7) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 8) . Filling cusp with (1, 8) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (3, 8) . Filling cusp with (3, 8) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 9) . Filling cusp with (1, 9) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (2, 9) . Filling cusp with (2, 9) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 10) . Filling cusp with (1, 10) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (3, 10) . Filling cusp with (3, 10) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 11) . Filling cusp with (1, 11) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (2, 11) . Filling cusp with (2, 11) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (3, 11) . Filling cusp with (3, 11) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 12) . Filling cusp with (1, 12) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (2, 9) . Filling cusp with (2, 9) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 10) . Filling cusp with (1, 10) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (3, 10) . Filling cusp with (3, 10) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (5, 12) . Filling cusp with (5, 12) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\n","done":false}︡{"done":true}︡
︠︡{"stdout":"Checking knots 4-23 and mirror 5-1.\n\n","done":false}︡{"stdout":"Solid Torus\n","done":false}︡{"stdout":"4","done":false}︡{"stdout":" isomorphisms between their canonical triangulations found.\n\n","done":false}︡{"stdout":"Lens space (1, 2) . Filling cusp with (1, 2) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 3) . Filling cusp with (1, 3) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 7) . Filling cusp with (1, 7) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 4) . Filling cusp with (1, 4) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (2, 7) . Filling cusp with (2, 7) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 8) . Filling cusp with (1, 8) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (3, 8) . Filling cusp with (3, 8) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 9) . Filling cusp with (1, 9) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 5) . Filling cusp with (1, 5) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (2, 5) . Filling cusp with (2, 5) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}︡{"stdout":"\nLens space (1, 6) . Filling cusp with (1, 6) Dehn filling.\nKnots are not equivalent (0 isomorphisms between their canonical triangulations).\n","done":false}0df37f62-b769-42c0-a6e5-7be27d1fc815︠

︠99efee99-b349-45af-b3a2-146c12a14424︠

︠6a4482f5-f648-4722-b425-d55dd142e650︠

︡768ff35e-abd8-47d0-a1ab-46c85f43bd8f︡{"stdout":"q^(9/2) - q^(7/2) + q^(5/2) - 2*q^(3/2) + sqrt(q) - 2/sqrt(q) + 1/q^(3/2) - 1/q^(5/2)"}︡{"stdout":"\n"}︡{"stdout":"-1\n"}︡{"stdout":"(t1^4*t2 - t1^4 + t1^3 - t1^2*t2 + t1*t2^2 - t2^2 + t2)/(t1^2*t2)"}︡{"stdout":"\n"}︡{"stdout":"3\n"}︡
︠1ca3e601-8b0c-4241-9918-25674f2ba758︠









