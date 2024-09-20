cd ..
make

for gm_add in {100,200,300}; do
    # d_parent = (-3, -7, -11, -19, -23, -31, -43)
    #gm_add=100
    fn="experiments/d64_re_max_gm_add_${gm_add}_no_hard"
    sage trees_generation.sage 21 77 209 437 713 1333 --gm-add=$gm_add --no-hard > ${fn}_trees.log 2>&1
    sage testrelations.sage > ${fn}_rels.log 2>&1
    sage relations_norms.sage > ${fn}_norms.log 2>&1
    #sage snf_precompute.sage > ${fn}_snf.log 2>&1
    zip -9 -r ${fn}_trees_rels.zip trees relations ${fn}_rels.log ${fn}_norms.log ${fn}_trees.log ${fn}_snf.log


    # d_parent = (-3, -7, -11, -19, -23, -31, -43, -47)
    #gm_add=100
    fn="experiments/d128_re_max_gm_add_${gm_add}_no_hard"
    sage trees_generation.sage 21 77 209 437 713 1333 2021 --gm-add=$gm_add --no-hard > ${fn}_trees.log 2>&1
    sage testrelations.sage > ${fn}_rels.log 2>&1
    sage relations_norms.sage > ${fn}_norms.log 2>&1
    #sage snf_precompute.sage > ${fn}_snf.log 2>&1
    zip -9 -r ${fn}_trees_rels.zip trees relations ${fn}_rels.log ${fn}_norms.log ${fn}_trees.log ${fn}_snf.log
done