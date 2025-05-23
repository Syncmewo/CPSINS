#include "..\CPSINS\PSINSCore\PSINS.h"
#include "..\CPSINS\UserDef\point_calib.h"
#include "..\eigen-3.4.0\Eigen\Eigen"

#define CalibTime  2 * 3600
#define CalibiTime 0.2
#define DataLen (int)(CalibTime/CalibiTime)

enum {

	NoCalib = 0,

	CalibDataRecing = 1,
	CalibFinished = 2,

	Waiting = 4,

};

double t = 0.0;
double SsSf, SsCf, CsCf, CsSf, Rce, Rse;
double Cs, Ss, Ce, Se, Cf, Sf;
Eigen::MatrixXd H(DataLen * 2, 7);
Eigen::VectorXd X(7,1), Y(DataLen * 2, 1);

double phiE, PhiE0, PhiN, PhiN0, PhiU, PhiU0, EbN, EbU, dVnE0, dVnN0;

int PSINS_SINS_PointCalib(CSINS& SINS, CVect3& pos, CVect3& vn, int& iter) {

	double lat = pos.i, cl = cos(lat), sl = sin(lat), el = 1.0 / cl, tl = tan(lat);
	double ws = sqrt(glv.g0 / glv.Re), wie = glv.wie, wf = glv.wie * sin(lat), SeuratT = 2 * PI / ws;
	double R = glv.Re, VI = sqrt(glv.g0 * glv.Re);
	double wn = wie * cl, wu = wie * sl;

	if (!IsZero(pos)) {

		if (iter == 0) {
		
			t = SINS.tk;

			int iter2 = 2 * iter, iter2p1 = iter2 + 1;

			Cs = cos(ws * t), Ss = sin(ws * t);
			Ce = cos(wie * t), Se = sin(wie * t);
			Cf = cos(wf * t), Sf = sin(wf * t);

			SsSf = Ss * Sf; SsCf = Ss * Cf;CsCf = Cs * Cf; CsSf = Cs * Sf;
			Rce = R * Se; Rse = R * Se;

			H(iter2,0)=VI*SsSf;H(iter2,1)=-VI*SsCf;H(iter2,2)=-wn*sl*Rse;H(iter2,3)=CsCf;H(iter2,4)=CsSf;H(iter2,5)=cl*R;H(iter2,6)=-sl*Rce;
			H(iter2p1,0)=VI*SsCf;H(iter2p1, 1)=VI*SsSf;H(iter2p1,2)=-wn*Rce;H(iter2p1,3)=-CsSf;H(iter2p1,4)=CsCf;H(iter2p1,5)=0.0;H(iter2p1,6)=Rse;

			Y(iter2) = SINS.vn.i;   // dVnE
			Y(iter2p1) = SINS.vn.j;//  dVnN

			iter++;

			return CalibDataRecing;
		
		}

		if (iter > 0 && iter < DataLen) {
		
			if ((SINS.tk - t) < CalibiTime) return Waiting;
			
			t = SINS.tk;

			int iter2 = 2*iter, iter2p1 = iter2 + 1;

			Cs = cos(ws * t), Ss = sin(ws * t);
			Ce = cos(wie * t), Se = sin(wie * t);
			Cf = cos(wf * t), Sf = sin(wf * t);

			SsSf = Ss * Sf; SsCf = Ss * Cf; CsCf = Cs * Cf; CsSf = Cs * Sf;
			Rce = R * Se; Rse = R * Se;

			H(iter2, 0) = VI * SsSf; H(iter2, 1) = -VI * SsCf; H(iter2, 2) = -wn * sl * Rse; H(iter2, 3) = CsCf; H(iter2, 4) = CsSf; H(iter2, 5) = cl * R; H(iter2, 6) = -sl * Rce;
			H(iter2p1, 0) = VI * SsCf; H(iter2p1, 1) = VI * SsSf; H(iter2p1, 2) = -wn * Rce; H(iter2p1, 3) = -CsSf; H(iter2p1, 4) = CsCf; H(iter2p1, 5) = 0.0; H(iter2p1, 6) = Rse;

			Y(iter2) = SINS.vn.i;   // dVnE
			Y(iter2p1) = SINS.vn.j;//  dVnN

			iter++;

			if (iter < DataLen) return CalibDataRecing;
		}

		if (iter == DataLen) {
			
			//LSEfit
			Eigen::MatrixXd HT = H.transpose();
			Eigen::MatrixXd iHTH = (HT * H);
			X = iHTH.inverse() * HT *Y;

		
			PhiE0 = X(0); PhiN0 = X(1); PhiU0 = X(2);
			EbN = X(5) * cl + X(6) * sl;
			EbU = X(5) * sl - X(6) * cl;

			dVnE0 = X(3); dVnN0 = X(5) - PhiU0 * wu * R;

			phiE = PhiE0 * CsCf + PhiN0 * CsSf - PhiU0 * wn * SsSf / ws + dVnE0 * SsSf / VI - dVnN0 * SsCf / VI - EbN * SsSf / ws;
			PhiN = (-PhiE0 * CsSf + PhiN0 * CsCf + PhiU0 * wn * SsCf / ws + dVnE0 * SsCf / VI - dVnN0 * SsSf / VI - EbN * SsCf / ws) * -1.0;
			PhiU = (PhiE0 * el * (Se - sl * CsCf) + PhiN0 * tl * (CsCf - Ce) + PhiU0 * (Ce + wu * SsSf / ws) +
				tl * dVnE0 * SsCf / VI + tl * dVnN0 * SsSf / VI +
				EbN * tl * (Se / wie - SsCf / ws) - EbU * Se / wie) * -1.0;
			printf("PhiE:%f,PhiN:%f,PhiU:%f\r\n", phiE / glv.sec, PhiN / glv.sec, PhiU / glv.sec);

			//У׼
			qdelphi(SINS.qnb, -CVect3(phiE, PhiN, PhiU));

			SINS.lvr.i = phiE; SINS.lvr.j = PhiN; SINS.lvr.k = PhiU;

			SINS.pos.i = pos.i;
			SINS.pos.j = pos.j;

			iter++;
			return CalibFinished;

		}
		
		if (iter > DataLen) {
			
			t = SINS.tk;

			Cs = cos(ws * t), Ss = sin(ws * t);
			Ce = cos(wie * t), Se = sin(wie * t);
			Cf = cos(wf * t), Sf = sin(wf * t);

			SsSf = Ss * Sf; SsCf = Ss * Cf; CsCf = Cs * Cf; CsSf = Cs * Sf;
			Rce = R * Se; Rse = R * Se;

			phiE = PhiE0 * CsCf + PhiN0 * CsSf - PhiU0 * wn * SsSf / ws + dVnE0 * SsSf / VI - dVnN0 * SsCf / VI - EbN * SsSf / ws;
			PhiN = (-PhiE0 * CsSf + PhiN0 * CsCf + PhiU0 * wn * SsCf / ws + dVnE0 * SsCf / VI - dVnN0 * SsSf / VI - EbN * SsCf / ws)* -1.0;
			PhiU = (PhiE0 * el * (Se - sl * CsCf) + PhiN0 * tl * (CsCf - Ce) + PhiU0 * (Ce + wu * SsSf / ws) +
				tl * dVnE0 * SsCf / VI + tl * dVnN0 * SsSf / VI +
				EbN * tl * (Se / wie - SsCf / ws) - EbU * Se / wie) * -1.0;
		
			SINS.lvr.i = phiE; SINS.lvr.j = PhiN; SINS.lvr.k = PhiU;
		
			return CalibFinished;
		}

	}
}
