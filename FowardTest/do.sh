sta=test
refE=0.035

dir1=target
dir2=tar
Ngroup=3
a=2.5 #gaussian parameter for rf
rayp=0.06 #ray parameter for rf

fcontrol=${sta}_forward.control

fmod=${sta}.mod
fa=template_a.dat
fp=template_p.dat
fg=template_g.dat
fe=template_e.dat
frf=template_rf.dat
echo "model ${Ngroup} ${fmod}" >${fcontrol}
echo "para in.para" >> ${fcontrol}
echo "disp R 4 p $fp g $fg e $fe a $fa"     >> ${fcontrol}
echo "rf $a $rayp $frf" >> ${fcontrol}
echo "rfweight 0.4" >> ${fcontrol}
echo "hk tar_hk.lst 1 1" >> ${fcontrol}
echo "hkweight 0.3 0.4 0.3" >> ${fcontrol}
echo "monol" >> ${fcontrol}
echo "Eweight 1. $refE" >> ${fcontrol}
echo "#model -1" >> ${fcontrol}
echo "#search 30" >> ${fcontrol}
echo "outdir $dir1 $dir2" >> ${fcontrol}
echo "end" >> ${fcontrol}

dir=/Volumes/X9Pro/GeoInverse/src
$dir/MC_main<<EOF 
${fcontrol} 1
EOF
