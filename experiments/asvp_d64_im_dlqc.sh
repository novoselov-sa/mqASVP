cd ..
make

gm_add=170
loops=10
fprecision=14000
finc=2000

m=64
n=6
# clgp
fn="experiments/clgp_d${m}_im_gm${gm_add}_no_hard_unramified_quadchars"
sage trees_generation.sage -3 -7 -11 -19 -23 -31 --gm-add=$gm_add --no-ramified-primes --no-hard > ${fn}_trees.log 2>&1
sage testrelations.sage --quadchars > ${fn}_rels.log 2>&1
sage relations_norms.sage > ${fn}_norms.log 2>&1
sage snf_precompute.sage > ${fn}_snf.log 2>&1
zip -9 -r ${fn}_trees_rels.zip trees relations ${fn}_rels.log ${fn}_norms.log ${fn}_trees.log ${fn}_snf.log

for ((i=1;i<=loops;i++)); do
    # dlog
    fn="experiments/dlog_d${m}_im_gm${gm_add}_no_hard_unramified_quadchars_dlqc_${i}"
    sage dlog_cyc_rat.sage --quadchars > ${fn}.log 2>&1
    zip -9 -r ${fn}.zip ${fn}.log experiments/dlogs
    
    # asvp
    fn="experiments/asvp_d${m}_im_gm${gm_add}_no_hard_unramified_quadchars_dlqc_${i}"
    sage ideal_asvp.sage ${n} -f ${fprecision} -fin ${finc} > ${fn}.log 2>&1
    zip -9 -r ${fn}.zip ${fn}.log
done
