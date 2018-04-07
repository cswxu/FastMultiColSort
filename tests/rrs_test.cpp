/*
 * rss_test.cpp
 *
 *  Created on: 28 Oct 2015
 *      Author: neo
 */

#include    "gtest/gtest.h"

#include	"src/rrs.h"
#include	<cassert>
#include	<cmath>
#include	<functional>
#include	<vector>

namespace multiAttrSort{

/**
 * Let's find the minimum of a function f(x,y)
 * in the area centered at origin with radius 10.
 * For simplicity, use Hamming distance here (not Euclidean)
 */
class Foo2DRRS : public RRS<std::vector<double>, double> {

public:
	typedef std::vector<double> PointType;
	typedef double	CostType;

	Foo2DRRS(std::function<double(double, double)> fn)
	: RRS<std::vector<double>, double>(0.99,0.1,0.99,0.8,0.5,0.001),
	  f_(fn){
		last_opt_.clear();
	}

	PointType GetRandomSample() override {
		double x = Dice() * 2 * radius_ - radius_;
		double y = Dice() * 2 * radius_ - radius_;
		return PointType({x,y});
	}

	CostType  GetCost(const PointType& point) override {
		assert(2 == point.size());
		return f_(point[0], point[1]);
	}

	bool StopCritera() override {
		std::cout << GetOptimalPoint()[0] << "," << GetOptimalPoint()[1] << ":" << GetOptimalCost() << std::endl;
		counter++;
		return counter > 1000;

	}

	PointType GetRandomNeighbor(const PointType& point, double ro) override {
		// Well ... this doesn't guarantee we are in the safe domain
		auto x1 = point[0] + (Dice()-0.5) * radius_ * ro;
		auto y1 = point[1] + (Dice()-0.5) * radius_ * ro;

		std::cout << "Neighbor:"
				<< point[0] << "," << point[1] << "=>"
				<< x1 << "," << y1 << std::endl;

		return PointType({x1,y1});
	}

protected:
	double radius_ = 10;
	std::function<double(double, double)> f_;

	PointType last_opt_;
	size_t counter = 0;

	double HammingDistance(PointType p1, PointType p2){
		return std::abs(p1[0] - p2[0]) + std::abs(p1[1] - p2[1]);
	}

};

TEST(RRSTest, Bowl){
	// f(x,y) = (x-1)^2 + (y-2)^2
	std::function<double(double, double)> bowl = [](double x, double y){return std::pow(x-1, 2) + std::pow(y-2,2);};
	Foo2DRRS foo(bowl);
	foo.Optimize();

	EXPECT_TRUE( std::abs(foo.GetOptimalCost()-(0)) < 1e-3);

}

TEST(RRSTest, Bowl2){
	// f(x,y) = (x-1)^2 + (y-2)^2 when x > 0
	// 		  = (x+5)^2 + (y+1)^2 when x < 0
	std::function<double(double, double)> bowl
			= [](double x, double y){
			return x > 0.0 ? std::pow(x-1, 2) + std::pow(y-2,2)
			: std::pow(x+5,2) + std::pow(y+1,2) - 1;};
	Foo2DRRS foo(bowl);
	foo.Optimize();

	EXPECT_TRUE( std::abs(foo.GetOptimalCost()-(-1)) < 1e-3);

}

}	// namespace
