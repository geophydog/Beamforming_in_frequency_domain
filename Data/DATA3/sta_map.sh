R=-5/60/-20/20
J=M10i
PS=sta.ps
PDF=sta.pdf

gmt gmtset FONT_ANNOT_PRIMARY 20,5
gmt gmtset MAP_FRAME_TYPE plain

gmt psxy -R$R -J$J -K -T > $PS
gmt pscoast -R -J -K -O -B5f0.5 -By5f0.5 -BWSen -Swhite -G220/220/220 >> $PS
saclst stlo stla f DATA2/*.Z.SAC | awk '{print $2,$3}' | \
    gmt psxy -R -J -K -O -W0.3p,100/100/100 -Fr0/0 >> $PS
saclst stlo stla f DATA2/*.Z.SAC | awk '{print $2,$3}' | \
    gmt psxy -R -J -K -O -Si0.75c -G130/130/130 -W1.5p >> $PS
echo "0 0" | gmt psxy -R -J -K -O -Sa1c -G50/50/50 -W1.5p >> $PS
gmt psxy -R -J -O -T >> $PS

ps2pdf $PS $PDF
evince $PDF
rm gmt.*
