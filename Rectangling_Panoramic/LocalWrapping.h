#pragma once

#ifndef LocalWrapping_h
#define LocalWrapping_h


#include"Common.h"

enum Border {
	BORDER_TOP = 0,
	BORDER_BOTTOM = 1,
	BORDER_LEFT = 2,
	BORDER_RIGHT = 3
};

enum SeamDirection {
	SEAM_VERTICAL = 0,
	SEAM_HORIZONTAL = 1
};

vector<vector<Coordinate>> Local_wrap(CVMat src, CVMat& wrap_img, CVMat mask);
CVMat Insert_local_seam(CVMat src, CVMat& mask, int* seam, SeamDirection seamdirection, pair<int, int> begin_end, bool shiftToend);
int* Get_local_seam(CVMat src, CVMat mask, SeamDirection seamdirection, pair<int, int> begin_end);
vector<vector<Coordinate>> Get_Local_wrap_displacement(CVMat src, CVMat mask);
pair<int, int> Choose_longest_border(CVMat src, CVMat mask, Border& direction);
void wrap_mesh_back(vector<vector<CoordinateDouble>>& mesh, vector<vector<Coordinate>> displacementMap, Config config);
vector<vector<CoordinateDouble>> get_rectangle_mesh(CVMat src, Config config);

#endif
