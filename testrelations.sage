import os
import sys
import nprofile
import units
import goodprime
import goodprimecheat
import mult
import div
import subsetprod
import powerprod
import sqrt
import field
import relations
import timeit
import fb
import trees
from sage.combinat.subset import SubsetsSorted
import argparse

def parse_files(d, dick, food = None, gens = None):
  if food == None:
    food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "f": 41, "g": 53, "h": 61}
  if gens == None:
    gens = list(food.keys())
  n = len(d)
  file = gens[0] + '_'
  aDone = False
  for i in SubsetsSorted(gens[1:n])[1:]:
    #print(i, file)
    count = 0
    counting = False
    counter = 0
    currentfile = file + ''.join(i)
    f = open("trees/" + currentfile + ".christine","r")
    #print('openning ', "trees/" + currentfile + ".christine")
    for line in f:
      if line[0] == 'F':
        line = line[:-1]
        spl = line.split(' ')
        lst = spl[2]
        field = []
        for elt in lst.split('_'):
          field += [prod(int(food[k]^elt.count(k)) for k,_ in food.items())]
        if len(field) == 1:
          #if field[0] == 5 and aDone == True:
          #if field[0] == food.values()[0] and aDone == True: 
          if field[0] == d[0] and aDone == True: 
            counting = False
            continue 
          #if field[0] == 5 and aDone == False: aDone = True 
          #if field[0] == food.values()[0] and aDone == False: aDone = True
          if field[0] == d[0] and aDone == False: aDone = True
          counting = True
          counter = field[0]
        else:
          counting = False
      else:
        if counting == True: 
          dick[counter] += 1/3
    #print 'dick:', dick
    f.close()
  #for k, v in dick.iteritems():
  for k, v in dick.items():
    dick[k] = ZZ(dick[k] + 1)
    #print k, v, dick
  return dick

def make_dictionary(d):
   dick = {}
   for i in Subsets(d)[1:]:
      dick[prod(i)] = 0
   return dick

#returns (# quad relations, dictionary of dick[gennorm] = (#rels Q(gen,), where in vec start), list of pointers, list of gens, dictionary of conjugation vectors)
def get_info(d, pars):
   I = [1]
   for di in d:
     I += [di*j for j in I]
   aism = 0
   point = [] #S: offsets of S-units of Q(sqrt(di)) in ais
   for di in I[1:]:
      point += [aism]
      #aism += pars[di]
      if di > 0:
        #S: quadratic real case, the number of S-unit group generators is |S| + 1
        pars[di] = (point[-1], pars[di])
      else:
        #S: quadratic imaginary case, the number of S-unit group generators is |S|
        K = NumberField(x^2 - di, names="z") 
        if K.unit_group(proof=False).order() == 2: #S: TODO, replace this by explicit conditions
          pars[di] = (point[-1], pars[di] - 1)
        else:
          print("[info] Quadratic subfield with di = {} has non-trivial units".format(di))
          pars[di] = (point[-1], pars[di]) #S: there are non-trivial units
      aism += pars[di][1]
			#print('pars interm:', pars)
   sigmas = {}
   for di in I[1:]:
     sigmas[di] = get_sigma(I, di, aism, pars)
   return aism, pars, point, I[1:], sigmas

def get_sigma(I, di, aism, pars):
  vec = aism*[1]
  for i in I:
    if i % di == 0:
      vec[pars[i][0]:pars[i][0] + pars[i][1]] = pars[i][1]*[-1]
  return vec

food = trees.get_food() # load food from trees

seed = 1
def testrelations(food, seed):
	d = tuple(sorted(food.values(), key=lambda di: abs(di)))
	gens = list(sorted(food.keys()))
	#print('d:', d, 'gens:', gens)

	K = field.field(d)
	print("Regulator: ", units.approxregulator(d))
	dick = make_dictionary(d)
	print("dick:", dick)
	pars = parse_files(tuple(d), dick, food = food)
	print('pars:', pars)
	aism, pars, points, I, sigmas = get_info(d, pars)
	print('pars:', pars)
	relations.set_vars(aism, sigmas, pars, d)
	t = walltime()
	file = "_".join(gens)
	res = relations.relations_internal(tuple(d), file, food = food, seed = seed)
	print('len(relations)', len(res))
	#cProfile.run('relations.relations_internal(tuple(d), file, food = food)')
	print("total computation: ", walltime(t))
	# f = open("relations/" + str(d)+'_relations','w')
	# for i in range(len(res)):
	# 	f.write(str(list(res[i].factors))+ "\n" )
	# f.close()

	return 1

parser = argparse.ArgumentParser(description='Computing class group for a multiquadratic field.')
parser.add_argument("food", type=str, nargs="*", help="dictionary of the form {'a': d_1, ..., 'e': d_n} containing information of multiquadratic field and corresponding labels for generators")
parser.add_argument("--store", action=argparse.BooleanOptionalAction, default=True, help="store (or not) relations matrices")

parser.add_argument('--shorten', '-s', choices=['hnf', 'lll', 'bkz'], default="lll", help="algorithm for shortening of intermediate relation matrices.")
parser.add_argument('--shorten-final', '-sf', choices=['hnf', 'lll', 'bkz'], default="lll", help="algorithm for shortening of final relation matrix.")
parser.add_argument('--bkz-bsize', type=Integer, help="Block size for BKZ (default = log(d_1 * ... * d_n) / 2)")

parser.add_argument("--schirokauer", action=argparse.BooleanOptionalAction, default=False, help="use Schirokauer maps")

parser.add_argument("--quadchars", action=argparse.BooleanOptionalAction, default=False, help="use additional quadratic characters to avoid obstructions")

args = parser.parse_args()

relations.set_shortening_alg(args.shorten)
relations.set_final_shortening_alg(args.shorten_final)

relations.USE_SCHIROKAUER_MAPS = args.schirokauer
relations.USE_ADDITIONAL_QUAD_CHARACTERS = args.quadchars

if args.food != []:
  food = eval(preparse(" ".join(args.food)))

trees_food = trees.get_food()

if trees_food != None:
  assert food == trees_food, f"Run trees generation first! {food} != {trees_food}."

if args.shorten == "bkz" or args.shorten_final == 'bkz':
  if args.bkz_bsize != None:
    relations.set_bkz_block_size(args.bkz_bsize)
  else:
    dd = prod(food.values()).abs()
    #bs = floor(log(dd)/2)
    bs = floor((log(dd))^(1/2) * log(log(dd)) * len(food.values()))
    # for small dd we use default value of block size
    if bs < 40:
      bs = 40
    relations.set_bkz_block_size(bs)

trees.clear_dir("relations")
print(pari.allocatemem(600*1024^3))
testrelations(food, seed)
#nprofile.output([goodprime,units,relations])
