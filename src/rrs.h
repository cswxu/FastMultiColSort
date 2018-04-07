/*
 * RRS.h
 *
 *  Created on: 28 Oct 2015
 *      Author: neo
 */

#ifndef SRC_RRS_H_
#define SRC_RRS_H_

#include	<algorithm>
#include	<cassert>
#include	<cmath>
#include	<ctime>
#include	<random>
#include	<vector>
#include	<iostream>

namespace multiAttrSort {

/**
 * A naive implementation of RRS black-box optimization
 * in Ye et. al, 2003
 */
template <class PointType, class CostType>
class RRS {
public:
	/** Read the paper for meaning for all these crab! */
	RRS(double p, double r, double q, double v, double c, double st);
	virtual ~RRS();

	/** Invoke the RRS algorithm */
	void Optimize();
	/** Retrieve the optimal point */
	PointType GetOptimalPoint();
	/** Retrieve the optimal cost */
	CostType GetOptimalCost();

protected:

	/** Generate a random sample from the global space */
	virtual PointType GetRandomSample() = 0;
	/** Cost function */
	virtual CostType  GetCost(const PointType& point) = 0;
	/** Stop criterion: should I stop? */
	virtual bool StopCritera() = 0;
	// Get a random neighbor of a given point within distance ro.
	// (ro is a fraction of the global distance range)
	virtual PointType GetRandomNeighbor(const PointType& point, double ro) = 0;

	/** A dice() between 0 and 1. */
	double Dice();

private:
	// parameters
	double p_,r_, q_, v_, c_, st_; // all 0~1 ?

	CostType opt_cost_;
	PointType opt_point_;

	std::default_random_engine generator_;
	std::uniform_real_distribution<double> distribution_;
};

template <class P, class C>
RRS<P,C>::RRS(double p, double r, double q, double v, double c, double st)
	: p_(p), r_(r), q_(q), v_(v), c_(c), st_(st),
	  generator_(std::time(0)),
	  distribution_(0.0, 1.0){

}

template <class P, class C>
RRS<P,C>::~RRS() {
	// TODO Auto-generated destructor stub
}



template <class P, class C>
void RRS<P,C>::Optimize(){
	// Define variables
	std::vector<C> F;

	// Calculate n and l
	const size_t n = std::ceil(std::log(1-p_)/std::log(1-r_));
	const size_t l = std::ceil(std::log(1-q_)/std::log(1-v_));
	assert(n>=1);
	assert(l>=1);

	std::cout << "n=" << n << ", l=" << l << std::endl;

	// Generate n random samples, take the one with min cost
	// Set it to x0 and yr, add cost to F
	opt_point_ = GetRandomSample();
	opt_cost_ = GetCost(opt_point_);
	for(size_t i=1; i < n; i++){
		auto x = GetRandomSample();
		auto cost = GetCost(x);
		if(cost < opt_cost_){
			opt_cost_ = cost;
			opt_point_ = x;
		}
	}

	auto x0 = opt_point_;
	auto cost0 = opt_cost_;
	auto yr = opt_cost_;
	F.push_back(opt_cost_);

	bool exploit_flag = true;
	size_t i = 0;

	// Main loop body
	while(!StopCritera()){

		// Exploitation process
		if(exploit_flag){
			size_t j=0;
			auto fc = GetCost(x0);
			auto xl = x0;
			auto ro = r_;
			while(ro > st_){
				auto x1 = GetRandomNeighbor(xl, ro);
				auto cost1 = GetCost(x1);
				if(cost1 < fc){	// find a better point, re-align
					xl = x1;
					fc = cost1;
					j = 0;
				}
				else{
					j++;
				}

				if(j == l){	// Fail to find a better point, shrink the sample space
					ro = c_ * ro;
					j = 0;
				}
			}
			exploit_flag = false;
			if(fc < opt_cost_){
				opt_cost_ = fc;
				opt_point_ = xl;
			}
		}	// end exploitation process


		// Exploration Process
		assert(!exploit_flag);
		auto x1 = GetRandomSample();
		auto cost1 = GetCost(x1);
		if(cost1 < cost0){
			x0 = x1;
			cost0 = cost1;
		}
		if(cost1 < yr){	// Find a promising point, set the flag to exploit
			exploit_flag = true;
		}

		if(i == n){
			// Update the exploitation threshold every n samples in the parameter space
			F.push_back(cost0);
			yr = std::accumulate(F.begin(), F.end(), C(0)) / F.size();
			i = 0;
		}

		i++;
	}

}

template <class P, class C>
P RRS<P,C>::GetOptimalPoint(){
	return opt_point_;
}

template <class P, class C>
C RRS<P,C>::GetOptimalCost(){
	return opt_cost_;
}

template<class PointType, class CostType>
double RRS<PointType, CostType>::Dice() {
	return distribution_(generator_);
}

} /* namespace multiAttrSort */



#endif /* SRC_RRS_H_ */
