/*
 * CD.h
 *
 *  Created on: Aug 11, 2016
 *      Author: lien
 */

#ifndef CD_H_
#define CD_H_

#include "CD3d_model.h"
#include "measure/CD3d_hull.h"
#include "RAPID.H"
#include "util/mathtool/Matrix.h"
#include <list>

class CD_RAPID
{
public:
    CD_RAPID();
    virtual ~CD_RAPID();

	int buildModel(const cd_m& m);
	int buildModel(const cd_hull& m);
	bool isInCollision(int m1, int m2);
  bool isInCollision(int m1, int m2, list< pair<int,int> >& collisions);

	//clear up everything
	void destroy();

protected:
    vector<RAPID_model*> m_models;
};

#endif /* CD_H_ */
