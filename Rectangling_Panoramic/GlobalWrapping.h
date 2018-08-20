#pragma once
#pragma once

#ifndef GlobalWrapping_h
#define GlobalWrapping_h

#include"Common.h"

#define clamp(x,a,b)    (  ((a)<(b))				\
? ((x)<(a))?(a):(((x)>(b))?(b):(x))	\
: ((x)<(b))?(b):(((x)>(a))?(a):(x))	\
)

struct Line_rotate {
	Vector2d pstart = Vector2d::Zero();
	Vector2d pend = Vector2d::Zero();
	double angle = 0;
	Line_rotate(Vector2d pstart, Vector2d pend, double angle) {
		this->pstart = pstart;
		this->pend = pend;
		this->angle = angle;
	}
};

struct BilinearWeights {
	double s;
	double t;
};

SpareseMatrixD_Row get_shape_mat(vector<vector<CoordinateDouble>> mesh, Config config);
SpareseMatrixD_Row get_vertex_to_shape_mat(vector<vector<CoordinateDouble>> mesh, Config config);
pair<SpareseMatrixD_Row, VectorXd> get_boundary_mat(CVMat src, vector<vector<CoordinateDouble>> mesh, Config config);
VectorXd get_vertice(int row, int col, vector<vector<CoordinateDouble>> mesh);
vector<vector<vector<LineD>>> init_line_seg(CVMat src, CVMat mask, Config config, vector < LineD > &lineSeg_flatten,vector<vector<CoordinateDouble>> mesh, vector<pair<int, double>>&id_theta, vector<double> &rotate_theta);
SpareseMatrixD_Row get_line_mat(CVMat src, CVMat mask, vector<vector<CoordinateDouble>> mesh,
	vector<double>rotate_theta, vector<vector<vector<LineD>>> lineSeg, vector<pair<MatrixXd, MatrixXd>>& BilinearVec,
	Config config, int &linenum, vector<bool>& bad);
#endif
