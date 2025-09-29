study="./data/pvals/exampleOR/"
study_random="./data/pvals/exampleRR/"
modules="./data/modules/ker_based/"
num_permutations=200

./fishnet_multi.sh \
    --study $study \
    --study-random $study_random \
    --modules $modules \
    --docker \
    --conda \
    --conda_env fishnet \
    --num-permutations $num_permutations

