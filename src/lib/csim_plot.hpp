//
//  csim_plot.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef csim_plot_hpp
#define csim_plot_hpp

#include "armadillo"

void draw_mat(arma::mat dMat, const char *title = "", const char *colormap = "rainbow");
void draw_mat(arma::mat dMat, double x0, double x1, double y0, double y1, const char *title = "", const char *colormap = "rainbow");
void plot_vec(arma::vec x, arma::vec y, const char *title, const char *type = "default");
void plot_vec(arma::vec x, arma::vec y, double x0, double x1, double y0, double y1, const char *title, const char *type = "default");


#endif /* csim_plot_hpp */
