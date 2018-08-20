#include"Common.h"
#include"LocalWrapping.h"
#include"GlobalWrapping.h"
#include "lsd.h"
#include "GL/glut.h"
#include<cmath>
#define WindowTitle  "OpenGLŒ∆¿Ì≤‚ ‘"
GLuint texGround;
vector<vector<CoordinateDouble>> outputmesh;
vector<vector<CoordinateDouble>> mesh;
CVMat img;
GLuint matToTexture(cv::Mat mat, GLenum minFilter = GL_LINEAR,
	GLenum magFilter = GL_LINEAR, GLenum wrapFilter = GL_REPEAT) {
	//cv::flip(mat, mat, 0);
	// Generate a number for our textureID's unique handle
	GLuint textureID;
	glGenTextures(1, &textureID);

	// Bind to our texture handle
	glBindTexture(GL_TEXTURE_2D, textureID);

	// Catch silly-mistake texture interpolation method for magnification
	if (magFilter == GL_LINEAR_MIPMAP_LINEAR ||
		magFilter == GL_LINEAR_MIPMAP_NEAREST ||
		magFilter == GL_NEAREST_MIPMAP_LINEAR ||
		magFilter == GL_NEAREST_MIPMAP_NEAREST)
	{
		//cout << "You can't use MIPMAPs for magnification - setting filter to GL_LINEAR" << endl;
		magFilter = GL_LINEAR;
	}

	// Set texture interpolation methods for minification and magnification
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);

	// Set texture clamping method
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapFilter);

	// Set incoming texture format to:
	// GL_BGR       for CV_CAP_OPENNI_BGR_IMAGE,
	// GL_LUMINANCE for CV_CAP_OPENNI_DISPARITY_MAP,
	// Work out other mappings as required ( there's a list in comments in main() )
	GLenum inputColourFormat = GL_BGR_EXT;
	if (mat.channels() == 1)
	{
		inputColourFormat = GL_LUMINANCE;
	}

	// Create the texture
	glTexImage2D(GL_TEXTURE_2D,     // Type of texture
		0,                 // Pyramid level (for mip-mapping) - 0 is the top level
		GL_RGB,            // Internal colour format to convert to
		mat.cols,          // Image width  i.e. 640 for Kinect in standard mode
		mat.rows,          // Image height i.e. 480 for Kinect in standard mode
		0,                 // Border width in pixels (can either be 1 or 0)
		inputColourFormat, // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
		GL_UNSIGNED_BYTE,  // Image data type
		mat.ptr());        // The actual image data itself

						   // If we're using mipmaps then generate them. Note: This requires OpenGL 3.0 or higher

	return textureID;
}
//
void display(void){
	glLoadIdentity();
	// «Â≥˝∆¡ƒª
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindTexture(GL_TEXTURE_2D, texGround);
	for (int row = 0; row < 20; row++) {
		for (int col = 0; col < 20; col++) {
			CoordinateDouble &coord = outputmesh[row][col];
			CoordinateDouble &localcoord = mesh[row][col];
			coord.row /= img.rows;
			coord.col /= img.cols;
			coord.row -= 0.5;
			coord.col -= 0.5;
			coord.row *= 2;
			coord.col *= 2;
			coord.row=clamp(coord.row, -1, 1);
			coord.col=clamp(coord.col, -1, 1);
			//cout << coord << " ";

			localcoord.row /= img.rows;
			localcoord.col /= img.cols;
			localcoord.row = clamp(localcoord.row, 0, 1);
			localcoord.col = clamp(localcoord.col, 0, 1);
		}
		// cout << endl;
	}
	//system("pause");
	
	for (int row = 0; row < 19; row++) {
		for (int col = 0; col < 19; col++) {
			CoordinateDouble local_left_top = mesh[row][col];
			CoordinateDouble local_right_top = mesh[row][col+1];
			CoordinateDouble local_left_bottom = mesh[row+1][col];
			CoordinateDouble local_right_bottom = mesh[row+1][col+1];
			

			CoordinateDouble global_left_top = outputmesh[row][col];
			CoordinateDouble global_right_top = outputmesh[row][col + 1];
			CoordinateDouble global_left_bottom = outputmesh[row + 1][col];
			CoordinateDouble global_right_bottom = outputmesh[row + 1][col + 1];
			
			//
			glBegin(GL_QUADS);
				glTexCoord2d(local_right_top.col, local_right_top.row); glVertex3d(global_right_top.col, -1*global_right_top.row, 0.0f);
				glTexCoord2d(local_right_bottom.col, local_right_bottom.row); glVertex3d(global_right_bottom.col, -1*global_right_bottom.row, 0.0f);
				glTexCoord2d(local_left_bottom.col, local_left_bottom.row);	glVertex3d(global_left_bottom.col, -1*global_left_bottom.row, 0.0f);
				glTexCoord2d(local_left_top.col, local_left_top.row);glVertex3d(global_left_top.col, -1*global_left_top.row, 0.0f);
			glEnd();
			
		}
	}
	/*
	int row = 18;
	int col = 18;
	CoordinateDouble local_left_top = mesh[row][col];
	CoordinateDouble local_right_top = mesh[row][col + 1];
	CoordinateDouble local_left_bottom = mesh[row + 1][col];
	CoordinateDouble local_right_bottom = mesh[row + 1][col + 1];


	CoordinateDouble global_left_top = outputmesh[row][col];
	CoordinateDouble global_right_top = outputmesh[row][col + 1];
	CoordinateDouble global_left_bottom = outputmesh[row + 1][col];
	CoordinateDouble global_right_bottom = outputmesh[row + 1][col + 1];
	glBegin(GL_QUADS);
	glTexCoord2d(local_right_top.col, local_right_top.row);
	glVertex3d(1, 0.7, 0.0f);
	glTexCoord2d(local_right_bottom.col, local_right_bottom.row);
	glVertex3d(1, -0.7, 0.0f);
	glTexCoord2d(local_left_bottom.col, local_left_bottom.row); 	
	glVertex3d(-1, -1, 0.0f);
	glTexCoord2d(local_left_top.col, local_left_top.row); 
	glVertex3d(-1, 1, 0.0f);
	glEnd();
	*/
	glutSwapBuffers();
}
int main(int argc, char* argv[]) {
	img = cv::imread("C:\\Users\\guyuchao\\source\\repos\\Rectangling_Panoramic\\Rectangling_Panoramic\\testimg\\1.jpg");
	

	double Time = (double)cvGetTickCount();
	//cv::resize(img, img, cv::Size(0, 0), 0.5, 0.5);
	CVMat scaled_img;
	cv::resize(img, scaled_img, cv::Size(0, 0), 0.5, 0.5);

	Config config(scaled_img.rows,scaled_img.cols,20,20);
	CVMat mask = Mask_contour(scaled_img);
	CVMat tmpmask;
	mask.copyTo(tmpmask);
	CVMat wrapped_img = CVMat::zeros(scaled_img.size(), CV_8UC3);

	vector<vector<Coordinate>> displacementMap = Local_wrap(scaled_img,wrapped_img,tmpmask);
	mesh = get_rectangle_mesh(scaled_img,config);
	//drawmesh(wrapped_img, mesh, config);
	//system("pasue");
	wrap_mesh_back(mesh,displacementMap,config);
	
	cout << "wrap back"<<endl;

	SpareseMatrixD_Row shape_energy = get_shape_mat(mesh,config);
	cout << "get shape energy"<<endl;
	SpareseMatrixD_Row Q = get_vertex_to_shape_mat(mesh,config);
	pair<SpareseMatrixD_Row, VectorXd> pair_dvec_B = get_boundary_mat(scaled_img, mesh, config);
	cout << "get border constraint" << endl;
	
	vector<pair<int, double>>id_theta;
	vector < LineD > line_flatten;
	vector<double> rotate_theta;
	vector<vector<vector<LineD>>> LineSeg = init_line_seg(scaled_img, mask, config, line_flatten, mesh, id_theta,rotate_theta);
	


	
	for (int iter = 1; iter <= 10; iter++) {
		cout << iter << endl;
		int Nl = 0;
		vector<pair<MatrixXd, MatrixXd>> BilinearVec;//need to update
		vector<bool> bad;
		SpareseMatrixD_Row line_energy = get_line_mat(scaled_img, mask, mesh, rotate_theta, LineSeg, BilinearVec, config, Nl, bad);
		cout << "get line energy" << "  " << Nl << endl;
		//combine
		double Nq = config.meshQuadRow*config.meshQuadCol;
		double lambdaB = INF;
		double lambdaL = 100;
		SpareseMatrixD_Row shape = (1 / Nq)*(shape_energy*Q);
		SpareseMatrixD_Row boundary = lambdaB * pair_dvec_B.first;
		SpareseMatrixD_Row line = (lambdaL / Nl)*(line_energy*Q);

		SpareseMatrixD_Row K = row_stack(shape, line);
		SpareseMatrixD_Row K2 = row_stack(K, boundary);

		//print_sparse_mat_row(shape_energy,1);
		//system("pause");
		VectorXd B = pair_dvec_B.second;
		VectorXd BA = VectorXd::Zero(K2.rows());
		BA.tail(B.size()) = lambdaB * B;

		SparseMatrixD K2_trans = K2.transpose();
		SparseMatrixD A = K2_trans * K2;
		VectorXd b = K2_trans * BA;
		VectorXd x;

		CSolve *p_A = new CSolve(A);
		x = p_A->solve(b);
		//update theta
		outputmesh = vector_to_mesh(x, config);
		int tmplinenum = -1;
		VectorXd thetagroup = VectorXd::Zero(50);
		VectorXd thetagroupcnt = VectorXd::Zero(50);
		for (int row = 0; row < config.meshQuadRow; row++) {
			for (int col = 0; col < config.meshQuadCol; col++) {
				vector<LineD> linesegInquad = LineSeg[row][col];
				int QuadID = row * config.meshQuadCol + col;
				if (linesegInquad.size() == 0) {
					continue;
				}
				else {
					VectorXd S = get_vertice(row, col, outputmesh);
					for (int k = 0; k < linesegInquad.size(); k++) {
						tmplinenum++;
						//cout << tmplinenum<<endl;
						if (bad[tmplinenum] == true) {
							continue;
						}
						//cout << tmplinenum;
						pair<MatrixXd, MatrixXd> Bstartend = BilinearVec[tmplinenum];
						MatrixXd start_W_mat = Bstartend.first;
						MatrixXd end_W_mat = Bstartend.second;
						Vector2d newstart = start_W_mat * S;
						Vector2d newend = end_W_mat * S;

						double theta = atan((newstart(1) - newend(1)) / (newstart(0) - newend(0)));
						double deltatheta = theta - id_theta[tmplinenum].second;
						if (isnan(id_theta[tmplinenum].second) || isnan(deltatheta)) {
							continue;
						}

						if (deltatheta > (PI / 2)) {
							deltatheta -= PI;
						}
						if (deltatheta < (-PI / 2)) {
							deltatheta += PI;
						}
						thetagroup(id_theta[tmplinenum].first) += deltatheta;
						thetagroupcnt(id_theta[tmplinenum].first) += 1;
						//cout << newstart << endl << endl << newend;
					}
				}
			}
		}

		//cal mean theta
		for (int ii = 0; ii < thetagroup.size(); ii++) {
			thetagroup(ii) /= thetagroupcnt(ii);

		}
		//update rotate_theta
		for (int ii = 0; ii < rotate_theta.size(); ii++) {
			rotate_theta[ii] = thetagroup[id_theta[ii].first];
		}
	}
	//cout << x;
	//system("pause");
	//vector<vector<CoordinateDouble>> outputmesh = vector_to_mesh(x,config);
	//drawmesh(scaled_img, outputmesh, config);
	//system("pause");
	enlarge_mesh(mesh, 2, config);
	enlarge_mesh(outputmesh, 2, config);
	//CVMat outputimg = CVMat::zeros(img.size(), CV_32FC3);
	//CVMat ouputcnt = CVMat::zeros(img.size(), CV_32FC3);


	/*
	for (int row = 0; row < config.meshQuadRow; row++) {
		cout << row << endl;
		for (int col = 0; col < config.meshQuadCol; col++) {

			VectorXd Vq = get_vertice(row, col, outputmesh);//x0,y0,x1,y1
			VectorXd Vo = get_vertice(row, col, mesh);//x0,y0
			double col_len = max(Vq(0), max(Vq(2), max(Vq(4), Vq(6)))) - min(Vq(0), min(Vq(2), min(Vq(4), Vq(6))));
			double row_len = max(Vq(1), max(Vq(3), max(Vq(5), Vq(7)))) - min(Vq(1), min(Vq(3), min(Vq(5), Vq(7))));
			double col_step = 1 / (4 * col_len);
			double row_step = 1 / (4 * row_len);
			//system("pause");
			for (double i = 0; i < 1; i += row_step) {
				for (double j = 0; j < 1; j += col_step) {
					double v1w = 1 - i - j + i * j;
					double v2w = j - i * j;
					double v3w = i - i * j;
					double v4w = i * j;
					MatrixXd matt(2, 8);
					matt << v1w, 0, v2w, 0, v3w, 0, v4w, 0,
						0, v1w, 0, v2w, 0, v3w, 0, v4w;
					VectorXd pout = matt * Vq;
					VectorXd pref = matt * Vo;
					if (int(pout(1)) >= 0 && int(pout(0)) >= 0 && int(pout(1)) < img.rows&&int(pout(0)) < img.cols) {
						colorPixel pixel = img.at<colorPixel>(int(pref(1)), int(pref(0)));
						cv::Vec3f pixelf = cv::Vec3f(float(pixel[0]), float(pixel[1]), float(pixel[2]));
						outputimg.at<cv::Vec3f>(int(pout(1)), int(pout(0))) = outputimg.at<cv::Vec3f>(int(pout(1)), int(pout(0))) + pixelf;
						ouputcnt.at<cv::Vec3f>(int(pout(1)), int(pout(0))) += cv::Vec3f(1, 1, 1);
					}
					else {
						//cout << "unfill";
					}
				}
			}
		}
	}
	*/
	//CVMat finaloutput = outputimg / (255 * ouputcnt);
	
	//drawmesh(finaloutput, outputmesh, config);
	//drawmesh(finaloutput, outpu, config);
	//fill_image(finaloutput);
	/*cv::namedWindow("Border", CV_WINDOW_AUTOSIZE);
	cv::imshow("Border", finaloutput);
	cv::waitKey(0);*/


	//glut
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(img.cols, img.rows);
	glutCreateWindow(WindowTitle);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);    // ∆Ù”√Œ∆¿Ì
	texGround = matToTexture(img);
	glutDisplayFunc(&display);   //◊¢≤·∫Ø ˝ 
	Time = (double)cvGetTickCount() - Time;

	printf("run time = %gms\n", Time / (cvGetTickFrequency() * 1000));//∫¡√Î
	glutMainLoop(); //—≠ª∑µ˜”√
					//
	system("pause");
	return 0;
}