import epistasiscomponents
from epistasiscomponents.epistasis_objects import sire 
from epistasiscomponents.epistasis_functions import gene_functions as gf

parent = sire.Sire("AAAAAAAAA", 'S100 test')

mutation = (3, "G")

mutant = gf.mutate(parent.sequence, mutation)

print(mutant)

print(parent.create_mutants(3))

