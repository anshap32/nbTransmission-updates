#!/bin/bash

#Arguments: infection (1) or sample (2) date and HH prob

qsub -N corrt2 CorrSimEstQsub.qsub 300 2 2

qsub -N corrt3 CorrSimEstQsub.qsub 300 2 3

qsub -N corrt4 CorrSimEstQsub.qsub 300 2 4

qsub -N corrt5 CorrSimEstQsub.qsub 300 2 5
