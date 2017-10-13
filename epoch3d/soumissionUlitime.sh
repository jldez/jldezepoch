#!/bin/bash

material="silica"
lambda=800

for F in 2.51
do

        for tau in 10
        do

                dir=$material"_"$lambda"_"$F"_"$tau
                F_line="fluence="$F
                tau_line="tau="$tau

                perl -i -pe "s/.*/$F_line/ if $. == 6" multifeedback.sh
                perl -i -pe "s/.*/$tau_line/ if $. == 7" multifeedback.sh

                bash multifeedback.sh
                qsub $dir/job_pulse1.sh

        done

done

