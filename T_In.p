set terminal pngcairo size 1600,800 enhanced font "Helvetica, 20"
set output "./T_In.png"


grid_color = "#d5e0c9"
text_color = "#6a6a6a"
my_font_title = "SVBasic Manual,30"

#unset key
set ylabel offset 2.8,0
set xlabel offset 0,0.7

my_font_label="SVBasic Manual,13"
my_font_tics="SVBasic Manual,12"

set xlabel textcolor rgb text_color font my_font_label
set ylabel textcolor rgb text_color font my_font_label
set grid lc rgb grid_color
set xr [0:1000]
set yr [90:180]

set xtics textcolor rgb text_color font my_font_tics
set ytics textcolor rgb text_color font my_font_tics
set format y "%4g"

set style fill transparent solid 0.3 bo
set grid


## color scheme ##
blue_1="#009999"
blue_2="#079191"
yellow_1="#EDBF00"
yellow_2="#F2C60C"
yellow_3="#FFCE00"
purple_1="#C40060"
purple_2="#C90A68"
purple_3="#D30068"
blue_3="#0FA3A3"
blue_4="#08A3A3"
blue_5="#008E8E"

set style line 1 	lt rgb blue_1		lw 2
set style line 2 	lt rgb blue_2 		lw 2	
set style line 3 	lt rgb yellow_1  	lw 2	
set style line 4 	lt rgb yellow_2	 	lw 2	
set style line 5 	lt rgb yellow_3		lw 2	
set style line 6 	lt rgb purple_1 	lw 2	
set style line 7 	lt rgb purple_2 	lw 2	
set style line 8 	lt rgb purple_3 	lw 2	
set style line 9 	lt rgb blue_3		lw 2	
set style line 10 	lt rgb blue_4		lw 2	
set style line 11	lt rgb blue_5		lw 2	

set ylabel "Temperature (K)"
set xlabel "IN"

plot 'out.txt' using 1:2 w lines title 'Surface temperature' ls 4,\
	 'out.txt' using 1:3 w lines title 'External temperature' ls 1