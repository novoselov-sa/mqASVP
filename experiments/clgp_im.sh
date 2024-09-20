cd ..
make

#sage trees_generation.sage -3 -7 -11 -19 -23 --gm-add=100 > experiments/d32_im_gm_add_100_trees.log 2>&1
#sage testrelations.sage > experiments/d32_im_gm_add_100_rels.log 2>&1
#tar -cf experiments/d32_im_gm_add_100_trees_rels.tar trees relations

#sage trees_generation.sage -3 -7 -11 -19 -23 -31 --gm-add=100 > experiments/d64_im_gm_add_100_trees.log 2>&1
#sage testrelations.sage > experiments/d64_im_gm_add_100_rels.log 2>&1
#tar -cf experiments/d64_im_gm_add_100_trees_rels.tar trees relations

# sage trees_generation.sage -3 -7 -11 -19 -23 -31 -43 --gm-add=300 --no-hard > experiments/d128_im_gm_add_300_no_hard_trees.log 2>&1
# sage testrelations.sage > experiments/d128_im_gm_add_300_no_hard_rels.log 2>&1
# zip -9 -r experiments/d128_im_gm_add_300_no_hard_trees_rels.zip trees relations experiments/d128_im_gm_add_300_no_hard_rels.log experiments/d128_im_gm_add_300_no_hard_trees.log

# sage trees_generation.sage -3 -7 -11 -19 -23 -31 -43 -47 --gm-add=450 --no-hard > experiments/d256_im_gm_add_450_no_hard_trees.log 2>&1
# sage testrelations.sage > experiments/d256_im_gm_add_450_no_hard_rels.log 2>&1
# zip -9 -r experiments/d256_im_gm_add_450_no_hard_trees_rels.zip trees relations experiments/d256_im_gm_add_450_no_hard_rels.log experiments/d256_im_gm_add_450_no_hard_trees.log

# sage trees_generation.sage -3 -7 -11 -19 -23 -31 -43 -47 --gm-add=550 --no-hard > experiments/d256_im_gm_add_550_no_hard_trees.log 2>&1
# sage testrelations.sage > experiments/d256_im_gm_add_550_no_hard_rels.log 2>&1
# zip -9 -r experiments/d256_im_gm_add_550_no_hard_trees_rels.zip trees relations experiments/d256_im_gm_add_550_no_hard_rels.log experiments/d256_im_gm_add_550_no_hard_trees.log


sage trees_generation.sage -3 -7 -11 -19 -23 -31 -43 -47 --gm-add=520 --no-hard > experiments/d256_im_gm_add_520_no_hard_trees.log 2>&1
sage testrelations.sage > experiments/d256_im_gm_add_520_no_hard_rels.log 2>&1
zip -9 -r experiments/d256_im_gm_add_520_no_hard_trees_rels.zip trees relations experiments/d256_im_gm_add_520_no_hard_rels.log experiments/d256_im_gm_add_520_no_hard_trees.log