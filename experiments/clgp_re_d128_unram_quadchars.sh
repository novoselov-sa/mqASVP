cd ..
make


for gm_add in {100,200,250,300}; do 
#echo $gm_add
fn="experiments/d128_re_gm${gm_add}_no_hard_unramified_quadchars"
sage trees_generation.sage 5 13 17 29 37 41 53 --gm-add=$gm_add --no-ramified-primes --no-hard > ${fn}_trees.log 2>&1
sage testrelations.sage --quadchars > ${fn}_rels.log 2>&1
sage relations_norms.sage > ${fn}_norms.log 2>&1
sage snf_precompute.sage > ${fn}_snf.log 2>&1
zip -9 -r ${fn}_trees_rels.zip trees relations ${fn}_rels.log ${fn}_norms.log ${fn}_trees.log ${fn}_snf.log
done