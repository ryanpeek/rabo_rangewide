# Tasks

## Completed:

### Admixture:

> `sbatch -t 3600 --mem=60G -p high 04_get_admix.sh all_rabo_100k 7`
 - Run for both 75k and 100k, this is not the filtered version so outliers are included (only ~6 individuals)

### PCA

 > `sbatch -p high -t 3600 --mem=32G --mail-type ALL --mail-user rapeek@ucdavis.edu 03_pca_ibs.sh all_rabo_filt_75k`

 - **all_rabo_filt_75k/100k** both run, same as above they've been cleaned filtered to remove outliers and contain >2 samples per site. No upper max.
 - **all_rabo_50k** run with >2 samples, does not have outliers removed.

### FST with realSFS SAFs (unfolded)

Using this code:

 > `sbatch -t 100 -p high 09_get_fst_safs.sh subpops_list_100k_filtered 100k`

This generates FST based on safs (instead of the 2DSFS).
 - run for 100k_filtered only, approx 63 sites, so !63 pairwise comparisons.

### FST with realSFS SAFs (Folded)

 > `sbatch -t 2700 -p high 09_get_fst_folded_safs.sh subpops_list_100k_filtered 100k`

Used this to generate Fst *without* 2DSFS.
This method generates a max likelihood file first from the saf.idx files (`realSFS file1 file2 > file3.ml`), then uses that to index and finally calc fsts (realSFS fst index file1 file2 file3.ml), then (realSFS fst stats *.idx > finalFSTout)

### Thetas (FOLDED)

 - Generate theta stats from folded SFS, **100k filtered** only.

 > `sbatch -t 200 05b_get_folded_sfs.sh subpops_list_100k_filtered 100k`

 > `sbatch -t 300 08a_get_diversity.sh subpops_list_100k_filtered 100k`

 > `sbatch -t 300 -p med 08b_get_theta_counts  all_rabo_100k_gw_theta`

To take mean of a col:

`awk '{ total += $2 } END { print total/NR }' yourFile.whatever`

Take SD of a col:
`awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}'`

results_thetas/ala-arroyomocho_100k_gw_thetasWin.gz.pestPG 0.00677853 0.00683759 -0.0734949 2857858


### Make 2DSFS

  > `sbatch -t 200 -p high 07b_get_folded_2dsfs.sh subpops_list_100k_filtered 100k`

 To generate 2DSFS for calculating FST, used the above.

### FST with 2DSFS

 > `realSFS fst index results_folded_sfs/${pop1}_${thresh}.folded.saf.idx results_folded_sfs/${pop2}_${thresh}.folded.saf.idx -sfs results_folded_sfs/${pop1}.${pop2}_${thresh}.folded.2dsfs -fstout results_fst/${pop1}.${pop2}_${thresh}.folded`

 > `realSFS fst stats results_fst/${pop1}.${pop2}_${thresh}.folded.fst.idx 2> results_fst/${pop1}.${pop2}_${thresh}_folded_global.fst`

 Need to generate 2DSFS first before using this to calc Fst.
