cd ..
make

gm_add=300
fn="experiments/dlog_d128_re_gm_add_${gm_add}_no_hard"
sage dlog_cyc_rat.sage > ${fn}.log 2>&1
zip -9 -r ${fn}.zip ${fn}.log experiments/dlogs
