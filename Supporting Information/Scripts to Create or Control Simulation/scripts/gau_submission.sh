#!/bin/bash --login
#$-cwd
#$-o grp_gpu.$JOB_ID.out
#$-N mptest
#$-j y
#$-l h_rt=172800
#$-pe mp 4

module load g09

g09 < opt.gau > opt.out
