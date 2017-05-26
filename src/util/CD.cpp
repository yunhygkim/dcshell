/*
 * CD.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: jml
 */

#include "CD.h"


CD_RAPID::CD_RAPID()
{
    //nothing to do here
}

CD_RAPID::~CD_RAPID()
{
	destroy();
}

//clear up everything
void CD_RAPID::destroy()
{
	for (auto cd_model : this->m_models)
		delete cd_model;

	this->m_models.clear();
}

int CD_RAPID::buildModel(const cd_m& m)
{
  auto rapid_model = new RAPID_model();

	const int tsize = m.f().size();

  rapid_model->BeginModel();
	int tid = 0;
	for (auto tri : m.f())
	{
		if(tri->area()>1e-10)
		{
			const auto& pt1 = tri->v[0]->pos;
			const auto& pt2 = tri->v[1]->pos;
			const auto& pt3 = tri->v[2]->pos;

	    double p1[3]={pt1[0],pt1[1],pt1[2]};
	    double p2[3]={pt2[0],pt2[1],pt2[2]};
	    double p3[3]={pt3[0],pt3[1],pt3[2]};

	    rapid_model->AddTri(p1,p2,p3,tid);
	  }
		tid++;
  }
  //end RAPID model

  rapid_model->EndModel();


  this->m_models.push_back(rapid_model);

  // return the cd id
  return m_models.size()-1;
}

int CD_RAPID::buildModel(const cd_hull& m)
{
	auto rapid_model = new RAPID_model();

	const int tsize = m.getFaceSize();

	rapid_model->BeginModel();
	//for (auto tri : m.f())
	for (int tid = 0; tid < tsize; tid++)
	{
		int vid1 = m.vid(tid, 0);
		int vid2 = m.vid(tid, 1);
		int vid3 = m.vid(tid, 2);

		const auto& pt1 = m.getVertex(vid1);
		const auto& pt2 = m.getVertex(vid2);
		const auto& pt3 = m.getVertex(vid3);

		double p1[3] = { pt1[0], pt1[1], pt1[2] };
		double p2[3] = { pt2[0], pt2[1], pt2[2] };
		double p3[3] = { pt3[0], pt3[1], pt3[2] };

		rapid_model->AddTri(p1, p2, p3, tid);
	}
	//end RAPID model

	rapid_model->EndModel();


	this->m_models.push_back(rapid_model);

	// return the cd id
	return m_models.size() - 1;
}

bool CD_RAPID::isInCollision(int m1, int m2)
{

	double T[3] = { 0, 0, 0 };
	double eye[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

	RAPIDres res;
	RAPID_Collide(res,
		eye, T, 1, m_models[m1],
		eye, T, 1, m_models[m2],
		RAPID_FIRST_CONTACT);

	return res.RAPID_num_contacts>0;
}

bool CD_RAPID::isInCollision(int m1, int m2, list< pair<int,int> >& collisions)
{

	double T[3] = { 0, 0, 0 };
	double eye[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

	RAPIDres res;
	RAPID_Collide(res,
		eye, T, 1, m_models[m1],
		eye, T, 1, m_models[m2],
		RAPID_ALL_CONTACTS);

  for(int i=0;i<res.RAPID_num_contacts;i++)
	{
		auto& c=res.RAPID_contact[i];
		collisions.push_back(make_pair(c.id1, c.id2));
	}

	return res.RAPID_num_contacts>0;
}
