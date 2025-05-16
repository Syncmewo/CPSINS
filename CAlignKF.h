#include "..\CPSINS\PSINSCore\PSINS.h"


class CAlignVn :public CKalman
{
public:
	CVect3 pos0;
	CQuat qnb;

	CAlignVn();
	void Init(const CSINS& sins0);

	void TimeUpdate(double ts);
	void MeasUpdate();

	void SetFtSins(CSINS& sins);
	void SetKFHk();
	void SetMeasVn(CSINS& sins);
	void FeedBack(CSINS& sins);


	virtual void SetFt(int nnq) {};			// process matrix setting
	virtual void SetHk(int nnq) {};			// measurement matrix setting
	virtual void SetMeas(void) {};				// set measurement
	virtual void Feedback(int nnq, double fbts) {};	// state feedback
};
