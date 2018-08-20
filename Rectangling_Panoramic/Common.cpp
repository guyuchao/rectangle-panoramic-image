#include"Common.h"

void fillHole(const CVMat srcBw, CVMat &dstBw) {
	cv::Size m_Size = srcBw.size();
	CVMat Temp = CVMat::zeros(m_Size.height + 2, m_Size.width + 2, srcBw.type());//—”’πÕºœÒ  
	srcBw.copyTo(Temp(cv::Range::Range(1, m_Size.height + 1), cv::Range::Range(1, m_Size.width + 1)));

	cv::floodFill(Temp, cv::Point(0, 0), cv::Scalar(255));

	CVMat cutImg;//≤√ºÙ—”’πµƒÕºœÒ  
	Temp(cv::Range::Range(1, m_Size.height + 1), cv::Range::Range(1, m_Size.width + 1)).copyTo(cutImg);

	dstBw = srcBw | (~cutImg);
}

CVMat Mask_contour(CVMat src) {
	CVMat bw;
	cv::cvtColor(src, src, CV_BGR2GRAY);
	

	uchar thr = 252;
	CVMat mask = CVMat::zeros(src.size(), CV_8UC1);
	for (int row = 0; row < src.rows; row++) {
		for (int col = 0; col < src.cols; col++) {
			if (src.at<uchar>(row, col)<thr) {
				mask.at<uchar>(row, col) = 255;
			}
		}
	}

	fillHole(mask, bw);
	bw = ~bw;
	
	CVMat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5, 5));
	CVMat dilate_out;//≈Ú’Õ
	cv::dilate(bw, dilate_out, element);
	cv::dilate(dilate_out, dilate_out, element);
	cv::dilate(dilate_out, dilate_out, element);

	CVMat erode_out;//∏Ø ¥
	erode(dilate_out, erode_out, element);
	return erode_out;
}

VectorXd mesh_to_vector(vector<vector<CoordinateDouble>> mesh, Config config) {
	int numMeshRow = config.meshNumRow;
	int numMeshCol = config.meshNumCol;
	VectorXd vec = VectorXd::Zero(numMeshRow*numMeshCol * 2);
	for (int row = 0; row < numMeshRow; row++) {
		for (int col = 0; col < numMeshCol; col++) {
			CoordinateDouble coord = mesh[row][col];
			vec((row*numMeshCol + col) * 2) = coord.col;
			vec((row*numMeshCol + col) * 2 + 1) = coord.row;
		}
	}
	return vec;
}

vector<vector<CoordinateDouble>> vector_to_mesh(VectorXd x, Config config) {
	int numMeshRow = config.meshNumRow;
	int numMeshCol = config.meshNumCol;
	vector<vector<CoordinateDouble>> mesh;
	for (int row = 0; row < numMeshRow; row++) {
		vector<CoordinateDouble> meshRow;
		for (int col = 0; col < numMeshCol; col++) {
			int xid = (row * numMeshCol + col) * 2;
			CoordinateDouble coord;
			coord.row = x(xid + 1);
			coord.col = x(xid);
			meshRow.push_back(coord);
		}
		mesh.push_back(meshRow);
	}
	return mesh;
}

void print_sparse_mat(SparseMatrixD Q) {
	for (int k = 0; k < Q.outerSize(); ++k) {
		for (SparseMatrixD::InnerIterator it(Q, k); it; ++it){
			std::cout << it.row() << " " << it.col() << " : " << it.value() << std::endl;
		}
		std::cout << std::endl;
	}
}

SpareseMatrixD_Row row_stack(SparseMatrixD origin, SpareseMatrixD_Row diag) {
	SpareseMatrixD_Row res(origin.rows() + diag.rows(), origin.cols());
	res.topRows(origin.rows()) = origin;
	res.bottomRows(diag.rows()) = diag;
	return res;
}
SpareseMatrixD_Row row_stack(SpareseMatrixD_Row origin, SpareseMatrixD_Row diag) {
	SpareseMatrixD_Row res(origin.rows() + diag.rows(), origin.cols());
	res.topRows(origin.rows()) = origin;
	res.bottomRows(diag.rows()) = diag;
	return res;
}

MatrixXd row_stack(MatrixXd mat1, MatrixXd mat2) {
	MatrixXd res(mat1.rows() + mat2.rows(), mat1.cols());
	res.topRows(mat1.rows()) = mat1;
	res.bottomRows(mat2.rows()) = mat2;
	return res;
}

MatrixXd col_stack(MatrixXd mat1, MatrixXd mat2) {
	MatrixXd res(mat1.rows(), mat1.cols() + mat2.cols());
	res.leftCols(mat1.cols()) = mat1;
	res.rightCols(mat2.cols()) = mat2;
	return res;
}

void DrawLine(CVMat& img, CoordinateDouble coordstart, CoordinateDouble coordend) {
	cv::Point start((int)coordstart.col, (int)coordstart.row);
	cv::Point end((int)coordend.col, (int)coordend.row);
	int thickness = 1;
	int lineType = 1;
	cv::line(img, start, end, cv::Scalar(0, 255, 0), thickness, lineType);
}

void DrawLine(CVMat& img, LineD line) {
	cv::Point start((int)line.col1, (int)line.row1);
	cv::Point end((int)line.col2, (int)line.row2);
	int thickness = 1;
	int lineType = 1;
	cv::line(img, start, end, cv::Scalar(0, 255, 0), thickness, lineType);
}

CVMat drawmesh(CVMat src, vector<vector<CoordinateDouble>> mesh, Config config) {
	int meshNumRow = config.meshNumRow;
	int meshNumCol = config.meshNumCol;

	for (int row = 0; row < meshNumRow; row++) {
		for (int col = 0; col < meshNumCol; col++) {
			CoordinateDouble now = mesh[row][col];
			if (row == meshNumRow - 1 && col<meshNumCol - 1) {
				CoordinateDouble right = mesh[row][col + 1];
				DrawLine(src, now, right);
			}
			else if (row < meshNumRow - 1 && col == meshNumCol - 1) {
				CoordinateDouble down = mesh[row + 1][col];
				DrawLine(src, now, down);
			}
			else if (row < meshNumRow - 1 && col < meshNumCol - 1) {
				CoordinateDouble right = mesh[row][col + 1];
				DrawLine(src, now, right);
				CoordinateDouble down = mesh[row + 1][col];
				DrawLine(src, now, down);
			}
		}
	}
	cv::namedWindow("Mesh", CV_WINDOW_AUTOSIZE);
	cv::imshow("Mesh", src);
	cv::waitKey(0);
	return src;
}

void enlarge_mesh(vector<vector<CoordinateDouble>>& mesh, double enlarge_factor, Config config) {
	int numMeshRow = config.meshNumRow;
	int numMeshCol = config.meshNumCol;
	for (int row = 0; row < numMeshRow; row++) {
		for (int col = 0; col < numMeshCol; col++) {
			CoordinateDouble &coord = mesh[row][col];

			coord.row = coord.row * enlarge_factor;
			coord.col = coord.col * enlarge_factor;
		}
	}
};

