#include "..\PSINSCore\PSINS.h"
#include "..\UserDef\point_calib.h"

enum {

	NoCalib = 0,

	FirstCalibRecing = 1,
	SecondCalibRecing = 2,
	CalibFinished = 3,

};

CVect3 pos_last_rec = O31, vn_last_rec = O31;
double kCalibInterval = 3000.0;
double t1 = 0.0, t2 = 0.0, t3 = 0.0;
double dPE[3] = { 0.0 }, dPN[3] = { 0.0 };
double dVE[3] = { 0.0 }, dVN[3] = { 0.0 };
double g_epsZ = 0.0;

CVect ZME = CVect(5, 1), ZMN = CVect(5, 1);
CVect XE = CVect(4, 1), XN = CVect(4, 1);

int PSINS_SINS_PointCalib1(CSINS& SINS, CVect3& pos, CVect3& vn, int& iter) {

	double lat = pos.i, cl = cos(lat);
	double phiE = 0.0, phiN = 0.0;
	double ws = sqrt(glv.g0 / glv.Re);
	double c1 = 0.0, c2 = 0.0, c3 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;

	//SINS.pos.j -= g_epsZ * SINS.nts;

	//如果用户发送一次外部位置和速度后，且在1小时内没有再次发送新的外部信息时，默认使用上一组数据；
	if (IsZero(pos)) {

		if (iter < 0 || iter > 2) return NoCalib;

		if (iter == 1 && (SINS.tk - t1 >= kCalibInterval)) {     //距离接受校正信息已过3000s

			pos = pos_last_rec;
			vn = vn_last_rec;

		}
		else if (iter == 1 && (SINS.tk - t1 < kCalibInterval)) {

			return FirstCalibRecing;

		}
		if (iter == 2 && (SINS.tk - t2 >= kCalibInterval)) { //距离接受校正信息已过6000s

			pos = pos_last_rec;
			vn = vn_last_rec;

		}
		else if (iter == 2 && (SINS.tk - t2 < kCalibInterval)) {

			return SecondCalibRecing;

		}

	}

	if (!IsZero(pos)) {

		pos_last_rec = pos;
		vn_last_rec = vn;

		if (iter == 0) {

			t1 = SINS.tk;

			dPE[iter] = (SINS.pos.i - pos.i) * SINS.eth.RMh; //北向位置误差
			dPN[iter] = (SINS.pos.j - pos.j) * SINS.eth.RNh * cl; //东向位置误差

			dVE[iter] = SINS.vn.i - vn.i;
			dVN[iter] = SINS.vn.j - vn.j;

			iter++;
			return FirstCalibRecing;
		}
		else if (iter == 1) {

			t2 = SINS.tk;

			if ((t2 - t1) < kCalibInterval) {

				iter = 0; // 重新开始第一个
				t1 = t2;
				dPE[iter] = (SINS.pos.i - pos.i) * SINS.eth.RMh; //北向位置误差
				dPN[iter] = (SINS.pos.j - pos.j) * SINS.eth.RNh * cl; //东向位置误差

				dVE[iter] = SINS.vn.i - vn.i;
				dVN[iter] = SINS.vn.j - vn.j;

				iter++;
				return FirstCalibRecing;

			}

			dPE[iter] = (SINS.pos.i - pos.i) * SINS.eth.RMh; //北向位置误差
			dPN[iter] = (SINS.pos.j - pos.j) * SINS.eth.RNh * cl; //东向位置误差

			dVE[iter] = SINS.vn.i - vn.i;
			dVN[iter] = SINS.vn.j - vn.j;

			iter++;
			return SecondCalibRecing;

		}
		else if (iter == 2) {

			t3 = SINS.tk;

			if ((t3 - t2) < kCalibInterval) {

				iter = 0; // 重新开始第一个
				t1 = t3;
				dPE[iter] = (SINS.pos.i - pos.i) * SINS.eth.RMh; //北向位置误差
				dPN[iter] = (SINS.pos.j - pos.j) * SINS.eth.RNh * cl; //东向位置误差

				dVE[iter] = SINS.vn.i - vn.i;
				dVN[iter] = SINS.vn.j - vn.j;

				iter++;
				return FirstCalibRecing;

			}

			c1 = cos(ws * t1); s1 = sin(ws * t1);
			c2 = cos(ws * t2); s2 = sin(ws * t2);
			c3 = cos(ws * t3); s3 = sin(ws * t3);

			double t21 = t2 - t1, t31 = t3 - t1;
			double t11 = t1 * t1, t22 = t2 * t2, t33 = t3 * t3;

			dPE[iter] = (SINS.pos.i - pos.i) * SINS.eth.RMh; //北向位置误差
			dPN[iter] = (SINS.pos.j - pos.j) * SINS.eth.RNh * cl; //东向位置误差

			dVE[iter] = SINS.vn.i - vn.i;
			dVN[iter] = SINS.vn.j - vn.j;

			CMat H = CMat(5, 4);
			H(0, 0) = 1; H(0, 1) = t1; H(0, 2) = c1; H(0, 3) = s1;
			H(1, 0) = 1; H(1, 1) = t2; H(1, 2) = c2; H(1, 3) = s2;
			H(2, 0) = 1; H(2, 1) = t3; H(2, 2) = c3; H(2, 3) = s3;
			H(3, 0) = t21; H(3, 1) = 0.5 * (t22 - t11); H(3, 2) = (s2 - s1) / ws; H(3, 3) = -(c2 - c1) / ws;
			H(4, 0) = t31; H(4, 1) = 0.5 * (t33 - t11); H(4, 2) = (s3 - s1) / ws; H(4, 3) = -(c3 - c1) / ws;


			ZME(0) = dVE[0]; ZME(1) = dVE[1]; ZME(2) = dVE[2];
			ZME(3) = dPE[1] - dPE[0]; ZME(4) = dPE[2] - dPE[0];

			ZMN(0) = dVN[0]; ZMN(1) = dVN[1]; ZMN(2) = dVN[2];
			ZMN(3) = dPN[1] - dPN[0]; ZMN(4) = dPN[2] - dPN[0];

			XE = lss(H, ZME); XN = lss(H, ZMN);

			//		    double XE0 = 0.0, XN0 = 0.0;
			//			XE0 = dPE[2] - XE(0) * t3 - 0.5 * XE(1) * t33 - XE(2) * s3 / ws - XE(3) * (1 - c3) / ws;
			//			XN0 = dPN[2] - XN(0) * t3 - 0.5 * XN(1) * t33 - XN(2) * s3 / ws - XN(3) * (1 - c3) / ws;

			//			CVect temp(3,1);
			//			temp(0) = (dPE[0] - (XE0 + XE(0) * t1 + 0.5 * XE(1) * t11 + XE(2) * s1 / ws + XE(3) * (1 - c1) / ws));
			//			temp(1) = (dPE[1] - (XE0 + XE(0) * t2 + 0.5 * XE(1) * t22 + XE(2) * s2 / ws + XE(3) * (1 - c2) / ws));
			//			temp(2) = (dPE[2] - (XE0 + XE(0) * t3 + 0.5 * XE(1) * t33 + XE(2) * s3 / ws + XE(3) * (1 - c3) / ws));
			//
			//			g_epsZ = -(temp(2) - temp(1)) / (t31 * glv.Re * cl);

			phiE = (XN(1) - XN(2) * ws * s3 + XN(3) * ws * c3) / glv.g0;
			phiN = -(XE(1) - XE(2) * ws * s3 + XE(3) * ws * c3) / glv.g0;
			//phiU = -ws * ws * (XN(2) * cos(ws * t3) + XN(3) * sin(ws * t3)) / (glv.g0 * wie * cos(lat)) + tan(lat) * phiN;

			printf("phiE:%f, phiN:%f, g_epsZ:%f\r\n", phiE / glv.sec, phiN / glv.sec,g_epsZ/glv.dph);

			//校准
			qdelphi(SINS.qnb, CVect3(phiE, phiN, 0.0));

			SINS.vn.i -= dVE[iter];
			SINS.vn.j -= dVN[iter];

			SINS.pos.i = pos.i;
			SINS.pos.j = pos.j;

			iter = 0;
			return CalibFinished;

		}

	}
}

