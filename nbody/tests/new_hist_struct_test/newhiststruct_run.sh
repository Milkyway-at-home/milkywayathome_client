# shell script for running new hist struct tests
# change to your own pathways

cd /mnt/c/users/"Emily Crook"/desktop/new_hist_struct/milkywayathome_client/build/bin

./milkyway_nbody \
-f /mnt/c/users/"Emily Crook"/desktop/new_hist_struct/milkywayathome_client/nbody/tests/new_hist_struct_test/newhiststruct_for_developers.lua \
-z /mnt/c/users/"Emily Crook"/desktop/new_hist_struct/milkywayathome_client/nbody/tests/new_hist_struct_test/input_hist.hist \
-o /mnt/c/users/"Emily Crook"/desktop/new_hist_struct/milkywayathome_client/nbody/tests/new_hist_struct_test/input_out.out \
-n 8 -b -e $13\
-i $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11\