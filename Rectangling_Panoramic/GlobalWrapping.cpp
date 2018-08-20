#include"GlobalWrapping.h"

MatrixXd BilinearWeightsToMatrix(BilinearWeights w) {
	MatrixXd mat(2,8);
	double v1w= 1 - w.s - w.t + w.s*w.t;
	double v2w = w.s - w.s*w.t;
	double v3w = w.t - w.s*w.t;
	double v4w = w.s*w.t;
	mat << v1w, 0, v2w, 0, v3w, 0, v4w, 0,
		0, v1w, 0, v2w, 0, v3w, 0, v4w;
	return mat;
}

BilinearWeights get_bilinear_weights(CoordinateDouble point, Coordinate upperLeftIndices, vector<vector<CoordinateDouble>> mesh){
	CoordinateDouble p1 = mesh[upperLeftIndices.row][upperLeftIndices.col]; // topLeft
	CoordinateDouble p2 = mesh[upperLeftIndices.row][upperLeftIndices.col + 1]; // topRight
	CoordinateDouble p3 = mesh[upperLeftIndices.row + 1][upperLeftIndices.col]; // bottomLeft
	CoordinateDouble p4 = mesh[upperLeftIndices.row + 1][upperLeftIndices.col + 1]; // bottomRight

	//test
	//cout << p1 << " " << p2 << endl << p3 << " " << p4 << endl<<point;
	//

	double slopeTop = (p2.row - p1.row) / (p2.col - p1.col);
	double slopeBottom = (p4.row - p3.row) / (p4.col - p3.col);
	double slopeLeft = (p1.row - p3.row) / (p1.col - p3.col);
	double slopeRight = (p2.row - p4.row) / (p2.col - p4.col);

	double quadraticEpsilon = 0.01;

	if (slopeTop == slopeBottom && slopeLeft == slopeRight) {

		// method 3
		Matrix2d mat1;
		mat1 << p2.col - p1.col, p3.col - p1.col,
			p2.row - p1.row, p3.row - p1.row;
		
		MatrixXd mat2(2,1);
		mat2 << point.col - p1.col, point.row - p1.row;

		MatrixXd matsolution = mat1.inverse()*mat2;

		BilinearWeights weights;
		weights.s = matsolution(0,0);
		weights.t = matsolution(1,0);
		return weights;
	}
	else if (slopeLeft == slopeRight) {

		// method 2
		double a = (p2.col - p1.col)*(p4.row - p3.row) - (p2.row - p1.row)*(p4.col - p3.col);
		double b = point.row*((p4.col - p3.col) - (p2.col - p1.col)) - point.col*((p4.row - p3.row) - (p2.row - p1.row)) + p1.col*(p4.row - p3.row) - p1.row*(p4.col - p3.col) + (p2.col - p1.col)*(p3.row) - (p2.row - p1.row)*(p3.col);
		double c = point.row*(p3.col - p1.col) - point.col*(p3.row - p1.row) + p1.col*p3.row - p3.col*p1.row;

		double s1 = (-1 * b + sqrt(b*b - 4 * a*c)) / (2 * a);
		double s2 = (-1 * b - sqrt(b*b - 4 * a*c)) / (2 * a);
		double s;
		if (s1 >= 0 && s1 <= 1) {
			s = s1;
		}
		else if (s2 >= 0 && s2 <= 1) {
			s = s2;
		}
		else {

			if ((s1 > 1 && s1 - quadraticEpsilon < 1) ||
				(s2 > 1 && s2 - quadraticEpsilon < 1)) {
				s = 1;
			}
			else if ((s1 < 0 && s1 + quadraticEpsilon > 0) ||
				(s2 < 0 && s2 + quadraticEpsilon > 0)) {
				s = 0;
			}
			else {
				// this case should not happen
				cerr << "   Could not interpolate s weight for coordinate (" << point.col << "," << point.row << ")." << endl;
				s = 0;
			}
		}

		double val = (p3.row + (p4.row - p3.row)*s - p1.row - (p2.row - p1.row)*s);
		double t = (point.row - p1.row - (p2.row - p1.row)*s) / val;
		double valEpsilon = 0.1; // 0.1 and 0.01 appear identical
		if (fabs(val) < valEpsilon) {
			// Py ~= Cy because Dy - Cy ~= 0. So, instead of interpolating with y, we use x.
			t = (point.col - p1.col - (p2.col - p1.col)*s) / (p3.col + (p4.col - p3.col)*s - p1.col - (p2.col - p1.col)*s);
		}

		BilinearWeights weights;
		weights.s = s;
		weights.t = t;		
		return weights;
	}
	else {

		// method 1
		double a = (p3.col - p1.col)*(p4.row - p2.row) - (p3.row - p1.row)*(p4.col - p2.col);
		double b = point.row*((p4.col - p2.col) - (p3.col - p1.col)) - point.col*((p4.row - p2.row) - (p3.row - p1.row)) + (p3.col - p1.col)*(p2.row) - (p3.row - p1.row)*(p2.col) + (p1.col)*(p4.row - p2.row) - (p1.row)*(p4.col - p2.col);
		double c = point.row*(p2.col - p1.col) - (point.col)*(p2.row - p1.row) + p1.col*p2.row - p2.col*p1.row;

		double t1 = (-1 * b + sqrt(b*b - 4 * a*c)) / (2 * a);
		double t2 = (-1 * b - sqrt(b*b - 4 * a*c)) / (2 * a);
		double t;
		if (t1 >= 0 && t1 <= 1) {
			t = t1;
		}
		else if (t2 >= 0 && t2 <= 1) {
			t = t2;
		}
		else {
			if ((t1 > 1 && t1 - quadraticEpsilon < 1) ||
				(t2 > 1 && t2 - quadraticEpsilon < 1)) {
				t = 1;
			}
			else if ((t1 < 0 && t1 + quadraticEpsilon > 0) ||
				(t2 < 0 && t2 + quadraticEpsilon > 0)) {
				t = 0;
			}
			else {
				// this case should not happen
				cerr << "   Could not interpolate t weight for coordinate (" << point.col << "," << point.row << ")." << endl;
				t = 0;
			}
		}

		double val = (p2.row + (p4.row - p2.row)*t - p1.row - (p3.row - p1.row)*t);
		double s = (point.row- p1.row - (p3.row - p1.row)*t) / val;
		double valEpsilon = 0.1; // 0.1 and 0.01 appear identical
		if (fabs(val) < valEpsilon) {
			// Py ~= Ay because By - Ay ~= 0. So, instead of interpolating with y, we use x.
			s = (point.col - p1.col - (p3.col - p1.col)*t) / (p2.col + (p4.col - p2.col)*t - p1.col - (p3.col - p1.col)*t);
		}

		BilinearWeights weights;
		weights.s = clamp(s, 0, 1);
		weights.t = clamp(t, 0, 1);
		return weights;
	}
}

pair<SpareseMatrixD_Row, VectorXd> get_boundary_mat(CVMat src, vector<vector<CoordinateDouble>> mesh, Config config) {
	//Vq=[x0 y0,x1,y1...]
	int rows = config.rows;
	int cols = config.cols;
	int numMeshRow = config.meshNumRow;
	int numMeshCol = config.meshNumCol;
	int vertexnum = numMeshRow * numMeshCol;

	VectorXd dvec = VectorXd::Zero(vertexnum * 2);
	VectorXd B = VectorXd::Zero(vertexnum * 2);
	for (int i = 0; i < vertexnum * 2; i += numMeshCol * 2) {//left
		dvec(i) = 1;
		B(i) = 0;
	}//x
	for (int i = numMeshCol * 2 - 2; i < vertexnum * 2; i += numMeshCol * 2) {//right
		dvec(i) = 1;
		B(i) = cols - 1;
	}//y

	for (int i = 1; i < 2 * numMeshCol; i += 2) {//top
		dvec(i) = 1;
		B(i) = 0;
	}

	for (int i = 2 * vertexnum - 2 * numMeshCol + 1; i < vertexnum * 2; i += 2) {//bottom
		dvec(i) = 1;
		B(i) = rows - 1;
	}

	//diag sparse;
	SpareseMatrixD_Row diag(dvec.size(), dvec.size());
	for (int i = 0; i < dvec.size(); i++) {
		diag.insert(i, i) = dvec(i);
	}
	diag.makeCompressed();
	return make_pair(diag, B);
};

VectorXd get_vertice(int row, int col, vector<vector<CoordinateDouble>> mesh) {//y0,x0,y1,x1...
	VectorXd Vq = VectorXd::Zero(8);
	CoordinateDouble p0 = mesh[row][col];//左上

	CoordinateDouble p1 = mesh[row][col + 1];//右上
	CoordinateDouble p2 = mesh[row + 1][col];//左下
	CoordinateDouble p3 = mesh[row + 1][col + 1];//右下
	Vq << p0.col, p0.row, p1.col, p1.row, p2.col, p2.row, p3.col, p3.row;
	return Vq;
}

SpareseMatrixD_Row get_shape_mat(vector<vector<CoordinateDouble>> mesh, Config config) {
	int numMeshRow = config.meshNumRow;
	int numMeshCol = config.meshNumCol;
	int numQuadRow = config.meshQuadRow;
	int numQuadCol = config.meshQuadRow;
	SpareseMatrixD_Row Shape_energy(8 * numQuadRow*numQuadCol, 8 * numQuadRow*numQuadCol);
	for (int row = 0; row < numQuadRow; row++) {
		for (int col = 0; col < numQuadCol; col++) {
			CoordinateDouble p0 = mesh[row][col];//左上
			CoordinateDouble p1 = mesh[row][col + 1];//右上
			CoordinateDouble p2 = mesh[row + 1][col];//左下
			CoordinateDouble p3 = mesh[row + 1][col + 1];//右下
			MatrixXd Aq(8, 4);//见论文
			Aq << p0.col, -p0.row, 1, 0,
				p0.row, p0.col, 0, 1,
				p1.col, -p1.row, 1, 0,
				p1.row, p1.col, 0, 1,
				p2.col, -p2.row, 1, 0,
				p2.row, p2.col, 0, 1,
				p3.col, -p3.row, 1, 0,
				p3.row, p3.col, 0, 1;
			//计算系数

			MatrixXd Aq_trans = Aq.transpose();
			MatrixXd Aq_trans_mul_Aq_reverse = (Aq_trans * Aq).inverse();
			MatrixXd I = MatrixXd::Identity(8, 8);
			MatrixXd coeff = (Aq*(Aq_trans_mul_Aq_reverse)*Aq_trans - I);

			int left_top_x = (row*numQuadCol + col) * 8;
			for (int i = 0; i < 8; i++) {
				for (int j = 0; j < 8; j++) {
					Shape_energy.insert(left_top_x + i, left_top_x + j) = coeff(i, j);
				}
			}
		}
	}
	Shape_energy.makeCompressed();
	return Shape_energy;
}

SpareseMatrixD_Row get_vertex_to_shape_mat(vector<vector<CoordinateDouble>> mesh, Config config) {
	int numMeshRow = config.meshNumRow;
	int numMeshCol = config.meshNumCol;
	int numQuadRow = config.meshQuadRow;
	int numQuadCol = config.meshQuadCol;
	SpareseMatrixD_Row Q(8 * numQuadRow*numQuadCol, 2 * numMeshRow*numMeshCol);
	for (int row = 0; row < numQuadRow; row++) {
		for (int col = 0; col < numQuadCol; col++) {
			int quadid = 8 * (row*numQuadCol + col);
			int topleftvertexId = 2 * (row*numMeshCol + col);
			Q.insert(quadid, topleftvertexId) = 1;
			Q.insert(quadid+1, topleftvertexId+1) = 1;
			Q.insert(quadid+2, topleftvertexId+2) = 1;
			Q.insert(quadid+3, topleftvertexId+3) = 1;
			Q.insert(quadid + 4, topleftvertexId + 2 * numMeshCol) = 1;
			Q.insert(quadid + 5, topleftvertexId + 2 * numMeshCol+1) = 1;
			Q.insert(quadid + 6, topleftvertexId + 2 * numMeshCol+2) = 1;
			Q.insert(quadid + 7, topleftvertexId + 2 * numMeshCol+3) = 1;
		}
	}
	Q.makeCompressed();
	/*
	for (int k = 0; k < Q.outerSize(); ++k) {
	for (SparseMatrix<double>::InnerIterator it(Q, k); it; ++it)
	{
	std::cout << it.row() << " " << it.col() << " : " << it.value() << std::endl;
	}
	system("pause");
	std::cout << std::endl;
	}*/
	return Q;
}

bool between(double a, double X0, double X1){
	double temp1 = a - X0;
	double temp2 = a - X1;
	if ((temp1<1e-8 && temp2>-1e-8) || (temp2<1e-6 && temp1>-1e-8)){
		return true;
	}
	else{
		return false;
	}
}

Vector2d detectIntersect(Matrix2d line1, Matrix2d line2, bool& isintersection) {
	double line_x = 0, line_y = 0; //交点  
	double p1_x = line1(0, 1), p1_y = line1(0, 0), p2_x = line1(1, 1), p2_y = line1(1, 0);
	double p3_x = line2(0, 1), p3_y = line2(0, 0), p4_x = line2(1, 1), p4_y = line2(1, 0);
	if ((fabs(p1_x - p2_x)<1e-6) && (fabs(p3_x - p4_x)<1e-6)) {
		isintersection = false;
	}
	else if ((fabs(p1_x - p2_x)<1e-6)) { //如果直线段p1p2垂直与y轴  
		if (between(p1_x, p3_x, p4_x)) {
			double k = (p4_y - p3_y) / (p4_x - p3_x);
			line_x = p1_x;
			line_y = k * (line_x - p3_x) + p3_y;

			if (between(line_y, p1_y, p2_y)) {
				isintersection = true;
			}
			else {
				isintersection = false;
			}
		}
		else {
			isintersection = false;
		}
	}
	else if ((fabs(p3_x - p4_x)<1e-6)) { //如果直线段p3p4垂直与y轴  
		if (between(p3_x, p1_x, p2_x)) {
			double k = (p2_y - p1_y) / (p2_x - p1_x);
			line_x = p3_x;
			line_y = k * (line_x - p2_x) + p2_y;
			if (between(line_y, p3_y, p4_y)) {
				isintersection = true;
			}
			else {
				isintersection = false;
			}
		}
		else {
			isintersection = false;
		}
	}
	else {
		double k1 = (p2_y - p1_y) / (p2_x - p1_x);
		double k2 = (p4_y - p3_y) / (p4_x - p3_x);

		if (fabs(k1 - k2)<1e-6) {
			isintersection = false;
		}
		else {
			line_x = ((p3_y - p1_y) - (k2*p3_x - k1 * p1_x)) / (k1 - k2);
			line_y = k1 * (line_x - p1_x) + p1_y;
		}

		if (between(line_x, p1_x, p2_x) && between(line_x, p3_x, p4_x)) {
			isintersection = true;
		}
		else {
			isintersection = false;
		}
	}
	Vector2d p;
	p << line_y, line_x;
	return p;
}

void revise_mask_for_lines(CVMat &mask) {
	//边缘的检测不予关注
	//对mask腐蚀
	int rows = mask.rows;
	int cols = mask.cols;
	for (int row = 0; row < rows; row++) {
		mask.at<uchar>(row, 0) = 255;
		mask.at<uchar>(row, cols - 1) = 255;
	}
	for (int col = 0; col < cols; col++) {
		mask.at<uchar>(0, col) = 255;
		mask.at<uchar>(rows - 1, col) = 255;
	}
	CVMat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(8, 8));
	cv::dilate(mask, mask, element);
	cv::dilate(mask, mask, element);
}

bool is_in_quad(CoordinateDouble point, CoordinateDouble topLeft, CoordinateDouble topRight,
	CoordinateDouble bottomLeft, CoordinateDouble bottomRight){
	// the point must be to the right of the left line, below the top line, above the bottom line,
	// and to the left of the right line

	// must be right of left line
	if (topLeft.col == bottomLeft.col) {
		if (point.col < topLeft.col) {
			return false;
		}
	}
	else {
		double leftSlope = (topLeft.row - bottomLeft.row) / (topLeft.col - bottomLeft.col);
		double leftIntersect = topLeft.row - leftSlope * topLeft.col;
		double yOnLineX = (point.row - leftIntersect) / leftSlope;
		if (point.col < yOnLineX) {
			return false;
		}
	}
	// must be left of right line
	if (topRight.col == bottomRight.col) {
		if (point.col > topRight.col) {
			return false;
		}
	}
	else {
		double rightSlope = (topRight.row - bottomRight.row) / (topRight.col - bottomRight.col);
		double rightIntersect = topRight.row - rightSlope * topRight.col;
		double yOnLineX = (point.row - rightIntersect) / rightSlope;
		if (point.col > yOnLineX) {
			return false;
		}
	}
	// must be below top line
	if (topLeft.row == topRight.row) {
		if (point.row < topLeft.row) {
			return false;
		}
	}
	else {
		double topSlope = (topRight.row - topLeft.row) / (topRight.col - topLeft.col);
		double topIntersect = topLeft.row - topSlope * topLeft.col;
		double xOnLineY = topSlope * point.col+topIntersect;
		if (point.row < xOnLineY) {
			return false;
		}
	}
	// must be above bottom line
	if (bottomLeft.row == bottomRight.row) {
		if (point.row > bottomLeft.row) {
			return false;
		}
	}
	else {
		double bottomSlope = (bottomRight.row - bottomLeft.row) / (bottomRight.col - bottomLeft.col);
		double bottomIntersect = bottomLeft.row - bottomSlope * bottomLeft.col;
		double xOnLineY = bottomSlope * point.col + bottomIntersect;
		if (point.row > xOnLineY) {
			return false;
		}
	}
	// if all four constraints are satisfied, the point must be in the quad
	return true;
}

bool line_in_mask(CVMat mask, LineD line) {
	if (mask.at<uchar>(line.row1, line.col1) == 0 && mask.at<uchar>(line.row2, line.col2) == 0) {
		return true;
	}
	return false;
}

vector<LineD> lsd_detect(CVMat src,CVMat mask) {
	//CVMat src2;
	//src.copyTo(src2);
	int rows = src.rows;
	int cols = src.cols;
	CVMat gray_img;
	cv::cvtColor(src, gray_img, CV_BGR2GRAY);
	double *image = new double[gray_img.rows*gray_img.cols];
	for (int row = 0; row < gray_img.rows; row++) {
		for (int col = 0; col < gray_img.cols; col++) {
			image[row*gray_img.cols + col] = gray_img.at<uchar>(row, col);
		}
	}
	vector<LineD> lines;
	double * out;
	int num_lines;
	out = lsd(&num_lines, image, gray_img.cols, gray_img.rows);
	for (int i = 0; i < num_lines; i++) {
		LineD line(out[i * 7 + 1], out[i * 7 + 0], out[i * 7 + 3], out[i * 7 + 2]);
		//CoordinateDouble start(out[i * 7 + 1], out[i * 7 + 0]);
		//CoordinateDouble end(out[i * 7 + 3], out[i * 7 + 2]);
		if (line_in_mask(mask, line)) {
			lines.push_back(line);
		}
		
		//DrawLine(src, start, end);
		/*DrawLine(src, start, end);
		cv::namedWindow("Border", CV_WINDOW_AUTOSIZE);
		cv::imshow("Border", src);
		cv::waitKey(0);*/
	}
	
	/*cv::namedWindow("Border", CV_WINDOW_AUTOSIZE);
	cv::imshow("Border", src);
	cv::waitKey(0);
	*/
	return lines;
}

bool does_segment_intersect_line(LineD lineSegment, double slope, double intersect,
	bool vertical, CoordinateDouble& intersectPoint){
	// calculate line segment m and b
	double lineSegmentSlope = INF;
	if (lineSegment.col1 != lineSegment.col2) {
		lineSegmentSlope = (lineSegment.row2 - lineSegment.row1) / (lineSegment.col2 - lineSegment.col1);
	}
	double lineSegmentIntersect = lineSegment.row1 - lineSegmentSlope * lineSegment.col1;

	// calculate intersection
	if (lineSegmentSlope == slope) {
		if (lineSegmentIntersect == intersect) {
			// same line
			intersectPoint.col = lineSegment.col1;
			intersectPoint.row = lineSegment.row1;
			return true;
		}
		else {
			return false;
		}
	}
	double intersectX = (intersect - lineSegmentIntersect) / (lineSegmentSlope - slope);
	double intersectY = lineSegmentSlope * intersectX + lineSegmentIntersect;
	// check if intersection is in the bounds of the line segment
	if (vertical) {
		if ((intersectY <= lineSegment.row1 && intersectY >= lineSegment.row2) ||
			(intersectY <= lineSegment.row2 && intersectY >= lineSegment.row1)) {
			intersectPoint.col = intersectX;
			intersectPoint.row = intersectY;
			return true;
		}
		else {
			return false;
		}
	}
	else {
		if ((intersectX <= lineSegment.col1 && intersectX >= lineSegment.col2) ||
			(intersectX <= lineSegment.col2 && intersectX >= lineSegment.col1)) {
			intersectPoint.col = intersectX;
			intersectPoint.row = intersectY;
			return true;
		}
		else {
			return false;
		}
	}
}

vector<CoordinateDouble> intersections_with_quad(LineD lineSegment, CoordinateDouble topLeft,
	CoordinateDouble topRight, CoordinateDouble bottomLeft,
	CoordinateDouble bottomRight){
	vector<CoordinateDouble> intersections;

	// left
	double leftSlope = INF;
	if (topLeft.col != bottomLeft.col) {
		leftSlope = (topLeft.row - bottomLeft.row) / (topLeft.col - bottomLeft.col);
	}
	double leftIntersect = topLeft.row - leftSlope * topLeft.col;
	// check
	CoordinateDouble leftIntersectPoint;
	if (does_segment_intersect_line(lineSegment, leftSlope, leftIntersect, true, leftIntersectPoint)) {
		if (leftIntersectPoint.row >= topLeft.row && leftIntersectPoint.row <= bottomLeft.row) {
			intersections.push_back(leftIntersectPoint);
		}
	}

	// right
	double rightSlope = INF;
	if (topRight.col != bottomRight.col) {
		rightSlope = (topRight.row - bottomRight.row) / (topRight.col - bottomRight.col);
	}
	double rightIntersect = topRight.row - rightSlope * topRight.col;
	// check
	CoordinateDouble rightIntersectPoint;
	if (does_segment_intersect_line(lineSegment, rightSlope, rightIntersect, true, rightIntersectPoint)) {
		if (rightIntersectPoint.row >= topRight.row && rightIntersectPoint.row <= bottomRight.row) {
			intersections.push_back(rightIntersectPoint);
		}
	}

	// top
	double topSlope = INF;
	if (topLeft.col != topRight.col) {
		topSlope = (topRight.row - topLeft.row) / (topRight.col - topLeft.col);
	}
	double topIntersect = topLeft.row - topSlope * topLeft.col;
	// check
	CoordinateDouble topIntersectPoint;
	if (does_segment_intersect_line(lineSegment, topSlope, topIntersect, false, topIntersectPoint)) {
		if (topIntersectPoint.col >= topLeft.col && topIntersectPoint.col <= topRight.col) {
			intersections.push_back(topIntersectPoint);
		}
	}

	// bottom
	double bottomSlope = INF;
	if (bottomLeft.col != bottomRight.col) {
		bottomSlope = (bottomRight.row - bottomLeft.row) / (bottomRight.col - bottomLeft.col);
	}
	double bottomIntersect = bottomLeft.row - bottomSlope * bottomLeft.col;
	// check
	CoordinateDouble bottomIntersectPoint;
	if (does_segment_intersect_line(lineSegment, bottomSlope, bottomIntersect, false, bottomIntersectPoint)) {
		if (bottomIntersectPoint.col >= bottomLeft.col && bottomIntersectPoint.col <= bottomRight.col) {
			intersections.push_back(bottomIntersectPoint);
		}
	}

	return intersections;
}

vector<vector<vector<LineD>>> segment_line_in_quad(CVMat src,vector<LineD> lines, vector<vector<CoordinateDouble>> mesh,Config config) {
	int QuadnumRow = config.meshQuadRow;
	int QuadnumCol = config.meshQuadCol;
	vector<vector<vector<LineD>>> quad_line_seg;
	CVMat src2;
	src.copyTo(src2);
	
	for (int row = 0; row < QuadnumRow; row++) {
		vector<vector<LineD>> vec_row;
		for (int col = 0; col < QuadnumCol; col++) {
			CoordinateDouble lefttop = mesh[row][col];
			CoordinateDouble righttop = mesh[row][col+1];
			CoordinateDouble leftbottom = mesh[row+1][col];
			CoordinateDouble rightbottom = mesh[row+1][col+1];
			
			vector<LineD> lineInQuad;
			for (int i = 0; i < lines.size(); i++) {
				LineD line = lines[i];
				CoordinateDouble point1(line.row1, line.col1);
				CoordinateDouble point2(line.row2, line.col2);
				bool p1InQuad = is_in_quad(point1, lefttop, righttop, leftbottom, rightbottom);
				bool p2InQuad = is_in_quad(point2, lefttop, righttop, leftbottom, rightbottom);
				if (p1InQuad && p2InQuad) {
					lineInQuad.push_back(line);
				}
				else if (p1InQuad) {
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					if (intersections.size() != 0) {
						LineD cutLine(point1,intersections[0]);
						lineInQuad.push_back(cutLine);
					}
				}
				else if (p2InQuad) {
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					if (intersections.size() != 0) {
						LineD cutLine(point2, intersections[0]);
						lineInQuad.push_back(cutLine);
					}
				}
				else {
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					if (intersections.size() ==2) {
						LineD cutLine(intersections[0], intersections[1]);
						lineInQuad.push_back(cutLine);
					}
				}
			}
			vec_row.push_back(lineInQuad);
			//TEST
			//DRAW BORDER
			//DrawLine(src, lefttop, righttop);
			//DrawLine(src, lefttop, leftbottom);
			//DrawLine(src, righttop, rightbottom);
			//DrawLine(src, leftbottom, rightbottom);
			/*for (int i = 0; i < lineInQuad.size(); i++) {
				LineD line = lineInQuad[i];
				DrawLine(src2, line);
			}
			cv::namedWindow("quad", CV_WINDOW_AUTOSIZE);
			cv::imshow("quad", src);
			cv::namedWindow("line", CV_WINDOW_AUTOSIZE);
			cv::imshow("line", src2);

			cv::waitKey(0);
			*/
		}
		quad_line_seg.push_back(vec_row);
	}
	return quad_line_seg;
}

void flatten(vector<vector<vector<LineD>>> lineSeg, vector<LineD>& line_vec, Config config) {
	int numQuadRow = config.meshQuadRow;
	int numQuadCol = config.meshQuadCol;
	for (int row = 0; row < numQuadRow; row++) {
		for (int col = 0; col < numQuadCol; col++) {
			for (int k = 0; k < lineSeg[row][col].size(); k++) {
				line_vec.push_back(lineSeg[row][col][k]);
			}
		}
	}
}

SpareseMatrixD_Row block_diag(SpareseMatrixD_Row origin, MatrixXd addin, int QuadID, Config config) {
	int cols_total = 8 * config.meshQuadRow*config.meshQuadCol;
	SpareseMatrixD_Row res(origin.rows() + addin.rows(), cols_total);
	res.topRows(origin.rows()) = origin;
	
	int lefttop_row = origin.rows();
	int lefttop_col = 8 * QuadID;
	for (int row = 0; row < addin.rows(); row++) {
		for (int col = 0; col < addin.cols(); col++) {
			res.insert(lefttop_row+row, lefttop_col+col) = addin(row,col);
		}
	}
	res.makeCompressed();
	return res;
	
}
vector<vector<vector<LineD>>> init_line_seg(CVMat src, CVMat mask,Config config, vector < LineD > &lineSeg_flatten,
	vector<vector<CoordinateDouble>> mesh, vector<pair<int, double>>&id_theta,vector<double> &rotate_theta ) {
	double thetaPerbin = PI / 49;
	revise_mask_for_lines(mask);
	vector<LineD> lines = lsd_detect(src, mask);//first step:detect line except border
												//step2: segment line in each quad
	vector<vector<vector<LineD>>> lineSeg = segment_line_in_quad(src,lines, mesh, config);
	//step3: construct sparsematrix
	flatten(lineSeg,lineSeg_flatten, config);
	
	for (int i = 0; i < lineSeg_flatten.size(); i++) {
		LineD line = lineSeg_flatten[i];
		double theta = atan((line.row1 - line.row2) / (line.col1 - line.col2));
		int lineSegmentBucket =(int) round((theta + PI / 2) / thetaPerbin);
		assert(lineSegmentBucket < 50);
		id_theta.push_back(make_pair(lineSegmentBucket, theta));
		rotate_theta.push_back(0);
	}
	return lineSeg;
}

SpareseMatrixD_Row get_line_mat(CVMat src, CVMat mask,vector<vector<CoordinateDouble>> mesh,
	vector<double>rotate_theta, vector<vector<vector<LineD>>> lineSeg,vector<pair<MatrixXd,MatrixXd>>& BilinearVec,
	Config config,int &linenum,vector<bool>& bad) {
	
	int linetmpnum = -1;
	int rows = config.rows;
	int cols = config.cols;
	int QuadnumRow = config.meshQuadRow;
	int QuadnumCol = config.meshQuadCol;
	double gridcols = config.colPermesh;
	double gridrows = config.rowPermesh;

	/*
	for (int i = 0; i < lines.size(); i++) {
		LineD line = lines[i];
		
		CoordinateDouble start(line.row1, line.col1);
		CoordinateDouble end(line.row2,line.col2 );
		if (line_in_mask(mask, line)) {
			DrawLine(src, start, end);
			cv::namedWindow("Border", CV_WINDOW_AUTOSIZE);
			cv::imshow("Border", src);
			cv::waitKey(0);
		}
	}*/

	
	SpareseMatrixD_Row energy_line;
	for (int row = 0; row < QuadnumRow; row++) {
		for (int col = 0; col < QuadnumCol; col++) {
			vector<LineD> linesegInquad = lineSeg[row][col];
			int QuadID = row * QuadnumCol + col;
			if (linesegInquad.size() == 0) {
				continue;
			}
			else {
				Coordinate topleft(row, col);
				MatrixXd C_row_stack(0,8);
				/*if (linesegInquad.size() > 2) {
					cout << endl<<QuadID<<" "<< linesegInquad.size();
					system("pause");
				}*/
				for (int k = 0; k < linesegInquad.size(); k++) {
					linetmpnum++;
					LineD line = linesegInquad[k];
					CoordinateDouble linestart(line.row1, line.col1);
					CoordinateDouble lineend(line.row2, line.col2);
					
					//test
					BilinearWeights startWeight = get_bilinear_weights(linestart, topleft, mesh);//s t2n t t1n
					MatrixXd start_W_mat = BilinearWeightsToMatrix(startWeight);
					BilinearWeights endWeight = get_bilinear_weights(lineend, topleft, mesh);
					MatrixXd end_W_mat = BilinearWeightsToMatrix(endWeight);
					//cout << startWeight.s << " " << startWeight.t << endl;//test
					//test
					VectorXd S= get_vertice(row, col, mesh);
					Vector2d ans = start_W_mat * S-Vector2d(linestart.col,linestart.row);
					Vector2d ans2 = end_W_mat * S - Vector2d(lineend.col, lineend.row);
					
					if (ans2.norm() >= 0.0001||ans.norm()>=0.0001) {//error case
						bad.push_back(true);
						BilinearVec.push_back(make_pair(MatrixXd::Zero(2,8), MatrixXd::Zero(2, 8)));
						continue;
					}
					assert(ans.norm() < 0.0001);
					assert(ans2.norm() < 0.0001);
					bad.push_back(false);
					//end test
					//system("pause");
					double theta = rotate_theta[linetmpnum];	
					BilinearVec.push_back(make_pair(start_W_mat, end_W_mat));
					Matrix2d R;
					R << cos(theta), -sin(theta), 
						sin(theta), cos(theta);
					MatrixXd ehat(2,1);
					ehat << line.col1 - line.col2, line.row1-line.row2;
					MatrixXd tmp = (ehat.transpose()*ehat).inverse();
					Matrix2d I = Matrix2d::Identity();
					MatrixXd C = R * ehat*tmp*(ehat.transpose())*(R.transpose()) - I;
					MatrixXd CT = C * (start_W_mat - end_W_mat);
					C_row_stack = row_stack(C_row_stack, CT);
				}
				energy_line = block_diag(energy_line, C_row_stack, QuadID,config);
			}
		}
	}
	linenum = linetmpnum;
	return energy_line;
}

vector < vector<vector<pair<int, double>>>> cal_theta(vector<vector<vector<Line_rotate>>> lineSeg, Config config) {
	vector < vector<vector<pair<int, double>>>> lineGroup;
	for (int row = 0; row < config.meshQuadRow; row++) {
		vector<vector<pair<int, double>>> row_vec;
		for (int col = 0; col < config.meshQuadCol; col++) {
			vector<pair<int, double>> vec;
			row_vec.push_back(vec);
		}
		lineGroup.push_back(row_vec);
	}
	double qstep = PI / 49;
	for (int row = 0; row < config.meshQuadRow; row++) {
		for (int col = 0; col < config.meshQuadCol; col++) {
			vector<Line_rotate> linevec = lineSeg[row][col];
			int linenum = linevec.size();
			for (int i = 0; i < linenum; i++) {
				Line_rotate line = linevec[i];
				Vector2d pstart = line.pstart;//y x
				Vector2d pend = line.pend;
				double theta = atan((pstart(0) - pend(0)) / (pstart(1) - pend(1)));
				int groupid = (int)(round((theta + PI / 2) / qstep) + 1);
				lineGroup[row][col].push_back(make_pair(groupid, theta));
			}
		}
	}
	return lineGroup;
}