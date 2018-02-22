//
//  csim_plot.cpp
//  csim
//
//  Created by steve on 4/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <assert.h>

#include "csim_plot.hpp"

void draw_mat(arma::mat dMat, const char *title, const char *colormap) {
    
    draw_mat(dMat, 0, dMat.n_cols-1, 0, dMat.n_rows-1, title, colormap);
}

void draw_mat(arma::mat dMat, double x0, double x1, double y0, double y1, const char *title, const char *colormap) {
    
    FILE *gnuDataFile = fopen("gnuData.txt", "w");
    for (int i = 0; i < dMat.n_rows; i++) {
        for (int j = 0; j < dMat.n_cols; j++) {
            fprintf(gnuDataFile, "%g ", dMat(i, j));
        }
        fprintf(gnuDataFile, "\n");
    }
    fclose(gnuDataFile);
    
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    // FILE * gnuplotPipe = stdout;
    fprintf(gnuplotPipe, "unset key \n");
    fprintf(gnuplotPipe, "set size ratio -1 \n");
    fprintf(gnuplotPipe, "set title \"%s\" \n", title);
    fprintf(gnuplotPipe, "set xrange [%f:%f] \n", x0, x1);
    fprintf(gnuplotPipe, "set yrange [%f:%f] \n", y0, y1);
    if (!strcmp(colormap, "default"))
        fprintf(gnuplotPipe, "set palette rgb 7,5,15 \n"); // default
    else if (!strcmp(colormap, "hot"))
        fprintf(gnuplotPipe, "set palette rgb 21,22,23 \n"); // hot
    else if (!strcmp(colormap, "rainbow"))
        fprintf(gnuplotPipe, "set palette rgb 33,13,10 \n"); // rainbow
    else if (!strcmp(colormap, "matlab"))
        fprintf(gnuplotPipe, "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e') \n"); // parula, matlab default after 2014
    else if (!strcmp(colormap, "gray"))
        fprintf(gnuplotPipe, "set palette gray \n"); // gray
    else {
        std::cout << "unknown color map, setting to default" << std::endl;
        fprintf(gnuplotPipe, "set palette rgb 7,5,15 \n"); // default
    }
    double xb = (x1-x0)/(dMat.n_cols-1);
    double xa = x0/xb;
    double yb = (y1-y0)/(dMat.n_rows-1);
    double ya = y0/yb;
    fprintf(gnuplotPipe, "plot 'gnuData.txt' u (($1+%f)*%f):(($2+%f)*%f):3 matrix with image\n", xa, xb, ya, yb);
    pclose(gnuplotPipe);
}

void plot_vec(arma::vec x, arma::vec y, const char *title, const char *type) {
    plot_vec(x, y, -1, -1, -1, -1, title, type);
}

void plot_vec(arma::vec x, arma::vec y, double x0, double x1, double y0, double y1, const char *title, const char *type) {
    
    FILE *gnuDataFile = fopen("gnuData.txt", "w");
    for (int i = 0; i < x.n_elem; i++) {
        fprintf(gnuDataFile, "%g %g\n", x[i], y[i]);
    }
    fclose(gnuDataFile);
    
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
//    FILE * gnuplotPipe = fopen ("gnuplot.txt", "w");
    fprintf(gnuplotPipe, "unset key \n");
    fprintf(gnuplotPipe, "set title \"%s\" \n", title);
    
    if (x1 != -1 || x0 != x1)
        fprintf(gnuplotPipe, "set xrange [%f:%f] \n", x0, x1);
    if (y0 != -1 || y0 != y1)
        fprintf(gnuplotPipe, "set yrange [%f:%f] \n", y0, y1);
    
    if (!strcmp(type, "semilogy")) {
        fprintf(gnuplotPipe, "set logscale y \n");
        fprintf(gnuplotPipe, "set format y \"10^{%%+02T}\" \n");
    }
    fprintf(gnuplotPipe, "set grid xtics lt 0 lw 1 lc rgb \"#bbbbbb\" \n");
    fprintf(gnuplotPipe, "set grid ytics lt 0 lw 1 lc rgb \"#bbbbbb\" \n");

    fprintf(gnuplotPipe, "plot 'gnuData.txt' with lines \n");
    pclose(gnuplotPipe);
    
}
