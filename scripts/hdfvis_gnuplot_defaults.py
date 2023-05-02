# Defaults for gnuplot
hdfvis_gnuplot_anima = [
    "set style data lines",
    "set grid",
    'set xlabel "r"',
    'set ylabel "value of function"',
]
hdfvis_gnuplot_plot = hdfvis_gnuplot_anima
hdfvis_gnuplot_error = [
    'set xlabel "r"',
    'set ylabel "Log_2(|difference|)"',
    "set y2tics nomirror",
    "set y2tics",
    "set style data lines",
    "set grid",
    "set nokey",  #'set key left'\
]
hdfvis_gnuplot_error2D = [
    'set xlabel "r"',
    'set ylabel "Log_2(|difference|)"',
    "set y2tics nomirror",
    "set y2tics",
    "set style data lines",
    "set grid",
    "set key left",
]
hdfvis_gnuplot_plot2D = ["set zrange [-2:2]", "set style data lines"]

# Defaults for gnuplot output
# Please ensure that the file extension is appropriate for the terminal
hdfvis_gnuplot_anima_nodis = ["set terminal gif animate optimize"]
hdfvis_gnuplot_anima_nodis_ext = "gif"
hdfvis_gnuplot_plot_nodis = ["set terminal postscript eps color"]
hdfvis_gnuplot_plot_nodis_ext = "eps"
hdfvis_gnuplot_plot2D_nodis = hdfvis_gnuplot_plot_nodis
hdfvis_gnuplot_plot2D_nodis_ext = hdfvis_gnuplot_plot_nodis_ext
hdfvis_gnuplot_error_nodis = ["set terminal postscript eps color"]
hdfvis_gnuplot_error_nodis_ext = "eps"
hdfvis_gnuplot_error2D_nodis = hdfvis_gnuplot_error_nodis
hdfvis_gnuplot_error2D_nodis_ext = hdfvis_gnuplot_error_nodis_ext
