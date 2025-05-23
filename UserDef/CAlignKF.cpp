#include "CAlignKF.h" 

CAlignVn::CAlignVn():CKalman(15, 3)
{

}

void CAlignVn::Init(const CSINS &sins0)
{
	pos0 = sins0.pos;  qnb = sins0.qnb;

	Pk.SetDiag2(fXXZU(10, 30, MIN), fXXX(0.1), fdPOS(10.0), fDPH3(0.1), fUG3(200.0));
	Qt.Set2(fDPSH3(0.001), fUGPSHZ3(1.0), fOO9);
	//Xmax.Set(fINF9, fDPH3(0.1), fMG3(1.0));
	Rt.Set2(fXXX(0.01));
	FBTau.Set(fXXX(0.1), fXXX(0.1), fINF3, fXXX(0.1), fXXX(0.1));
}

void CAlignVn::SetFtSins(CSINS &sins) {

	sins.etm();
	Ft.SetMat3(0, 0, sins.Maa, sins.Mav, sins.Map), Ft.SetMat3(0, 9, -sins.Cnb);
	Ft.SetMat3(3, 0, sins.Mva, sins.Mvv, sins.Mvp), Ft.SetMat3(3, 12, sins.Cnb);
	Ft.SetMat3(6, 3, sins.Mpv, sins.Mpp);
	// 0-14 phi,dvn,dpos,eb,db

}

void CAlignVn::SetKFHk() {

	Hk(0, 3) = Hk(1, 4) = Hk(2, 5) = 1.0;    // dvn
	SetMeasMask(007);

}

void CAlignVn::TimeUpdate(double ts) {
	
	CKalman::TimeUpdate(ts);

}

void CAlignVn::MeasUpdate() {


	CKalman::MeasUpdate();

}

void CAlignVn::SetMeasVn(CSINS &sins) {

		*(CVect3*)&Zk.dd[0] = sins.mvn;
		SetMeasFlag(0007);
		sins.pos = pos0;
}

void CAlignVn::FeedBack(CSINS& sins) {

	int nq = CKalman::nq; double fbts = sins.nts;
	CKalman::Feedback(nq, fbts);
	//FBXk = Xk; Xk -= FBXk;
	sins.qnb -= *(CVect3*)&FBXk.dd[0];  sins.vn -= *(CVect3*)&FBXk.dd[3];  sins.pos -= *(CVect3*)&FBXk.dd[6];
	sins.eb += *(CVect3*)&FBXk.dd[9];	sins.db += *(CVect3*)&FBXk.dd[12];  // 0-14 phi,dvn,dpos,eb,db

	qnb = sins.qnb;
}


