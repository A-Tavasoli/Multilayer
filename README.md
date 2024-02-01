# Multilayer Network Analysis

This repository contains MATLAB code for generating and analyzing multilayer networks. 
The codes are used to reproduce the results presented in the paper titled "The Art of Interconnections: Achieving Maximum Algebraic Connectivity in Multilayer Networks."

## Guidelines for Running Codes

The cvx package is used to solve the convex primal and dual programmings: http://cvxr.com/cvx/ 


### Primal Programming

The file `main_cvx_Spectrum_of_Laplacian.m` is used to generate the results of primal programming in sections 2-5 of the paper. 
The input to the code is the adjacency matrix for each layer. 
Random graph generation for BA, ER, Geo, and WS types is also supported. 
To reproduce specific figures from the paper, follow these instructions:

- **Case 1:** Load `A1GeoCase1.mat` and `A2GeoCase1.mat` for figures 2 and 3.
- **Case 2:** Load `A1GeoCase2.mat` and `A2GeoCase2.mat` for figures D.13 and D.14.
- **Case 3:** Load `A1GeoCase3.mat` and `A2GeoCase3.mat` for figures D.15 and D.16.

### Dual Programming and Embedding

The file `main_Embeding.m` is used to generate the results of dual programming and embedding in sections 2-5 of the paper. 
Similar to primal programming, the input is the adjacency matrix for each layer, and random graph generation is supported. 
To reproduce specific figures from the paper, follow these instructions:

- **Case 1:** Load `A1GeoCase1.mat` and `A2GeoCase1.mat` for figure 4. Load `A1WSCase1.mat` and `A2WSCase1.mat` for figure D.17.
- **Case 2:** Load `A1GeoCase2.mat` and `A2GeoCase2.mat` for figure D.18.
- **Case 3:** Load `A1GeoCase3.mat` and `A2GeoCase3.mat` for figure D.19.

### Well-Interconnected Graphs

The file `main_WellConnectedGeneral.m` is used to generate the results in section 6 of the paper. 
The input and random graph generation are similar to the previous sections. 
For instance, to regenerate Figure 9, load `A_Small1.mat` and `A_Small2.mat`.



Feel free to explore and adapt the code for your own research or projects.
