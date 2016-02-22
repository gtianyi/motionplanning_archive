#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

#include <ompl/tools/benchmark/Benchmark.h>
#include <ompl/base/goals/GoalState.h>
#include "QuadrotorPlanning.hpp"
#include "SE3RigidBodyPlanning.hpp"
#include "config.hpp"

class QuadrotorSpatialGoal : public ompl::base::GoalState {
public:
	QuadrotorSpatialGoal(const ompl::base::SpaceInformationPtr &si, const ompl::base::State *state) : ompl::base::GoalState(si) {
		setState(state);
		se3State = state_->as<ompl::base::CompoundStateSpace::StateType>()->as<ompl::base::SE3StateSpace::StateType>(0);
	}

	virtual double distanceGoal(const ompl::base::State *state) const {
		auto s = state->as<ompl::base::CompoundStateSpace::StateType>()->as<ompl::base::SE3StateSpace::StateType>(0);
		double dx = s->getX() - se3State->getX();
		double dy = s->getY() - se3State->getY();
		double dz = s->getZ() - se3State->getZ();
		return sqrt(dx*dx + dy*dy + dz*dz);
	}
protected:
	const ompl::base::SE3StateSpace::StateType *se3State = NULL;
};

BenchmarkData quadrotorBenchmark(const FileMap &params) {
	ompl::app::QuadrotorPlanning *quadrotor = new ompl::app::QuadrotorPlanning();

	globalParameters.globalAppBaseControl = quadrotor;

	ompl::app::SE3RigidBodyPlanning *abstract = new ompl::app::SE3RigidBodyPlanning();

	globalParameters.globalAbstractAppBaseGeometric = abstract;

	ompl::control::SimpleSetupPtr quadrotorPtr(quadrotor);

	ompl::base::StateSpacePtr stateSpace(quadrotorPtr->getStateSpace());

	if(params.exists("EnvironmentBounds")) {
		auto boundsVals = params.doubleList("EnvironmentBounds");
		std::vector<double> low, high;
		for(unsigned int i = 0; i < boundsVals.size(); i++) {
			low.emplace_back(boundsVals[i++]);
			high.emplace_back(boundsVals[i]);
		}

		// set the bounds for the R^3 part of SE(3)
		ompl::base::RealVectorBounds bounds(3);
		for(unsigned int i = 0; i < 3; i++) {
			bounds.setLow(i, low[i]);
			bounds.setHigh(i, high[i]);
		}

		stateSpace->as<ompl::base::CompoundStateSpace>()->as<ompl::base::SE3StateSpace>(0)->setBounds(bounds);

		abstract->getStateSpace()->as<ompl::base::SE3StateSpace>()->setBounds(bounds);
	} else {
		OMPL_WARN("using default environment bounds");
	}

	// define start state
	ompl::base::ScopedState<ompl::base::SE3StateSpace> start(quadrotor->getGeometricComponentStateSpace());
	auto startLoc = params.doubleList("Start");

	start->setX(startLoc[0]);
	start->setY(startLoc[1]);
	start->setZ(startLoc[2]);
	start->rotation().setIdentity();

	// define goal state
	ompl::base::ScopedState<ompl::base::SE3StateSpace> goal(quadrotor->getGeometricComponentStateSpace());
	auto goalLoc = params.doubleList("Goal");
	goal->setX(goalLoc[0]);
	goal->setY(goalLoc[1]);
	goal->setZ(goalLoc[2]);
	goal->rotation().setIdentity();

	double goalRadius = params.doubleVal("GoalRadius");

	// set the start & goal states
	auto myGoal = new BlimpSpatialGoal(
		quadrotor->getSpaceInformation(),
		quadrotor->getFullStateFromGeometricComponent(goal).get());
	myGoal->setThreshold(goalRadius);
	auto goalPtr = ompl::base::GoalPtr(myGoal);
	quadrotorPtr->setGoal(goalPtr);

	abstract->setStartAndGoalStates(start, goal, goalRadius);

	struct passwd *pw = getpwuid(getuid());
	const char *homedir = pw->pw_dir;
	std::string homeDirString(homedir);

	auto agentMesh = params.stringVal("AgentMesh");
	auto environmentMesh = params.stringVal("EnvironmentMesh");

	quadrotor->setRobotMesh(homeDirString + "/gopath/src/github.com/skiesel/moremotionplanning/models/" + agentMesh);
	quadrotor->setEnvironmentMesh(homeDirString + "/gopath/src/github.com/skiesel/moremotionplanning/models/" + environmentMesh);

	abstract->setRobotMesh(homeDirString + "/gopath/src/github.com/skiesel/moremotionplanning/models/" + agentMesh);
	abstract->setEnvironmentMesh(homeDirString + "/gopath/src/github.com/skiesel/moremotionplanning/models/" + environmentMesh);

	if(params.exists("MinControlDuration")) {
		quadrotorPtr->getSpaceInformation()->setMinMaxControlDuration(params.integerVal("MinControlDuration"), params.integerVal("MaxControlDuration"));
	}

	if(params.exists("PropagationStepSize")) {
		quadrotorPtr->getSpaceInformation()->setPropagationStepSize(params.doubleVal("PropagationStepSize"));
	}

	if(!quadrotorPtr->getSpaceInformation()->getStatePropagator()->canSteer()) //if it can steer, leave it alone!
		quadrotorPtr->getSpaceInformation()->setDirectedControlSamplerAllocator(directedControlSamplerAllocator);

	quadrotorPtr->setup();
	abstract->setup();

	globalParameters.abstractBounds = abstract->getStateSpace()->as<ompl::base::SE3StateSpace>()->getBounds();
	globalParameters.abstractBounds.resize(6);
	for(unsigned int i = 3; i < 6; i++) {
		globalParameters.abstractBounds.setLow(i, -M_PI);
		globalParameters.abstractBounds.setHigh(i, M_PI);
	}

	globalParameters.copyVectorToAbstractState = [](ompl::base::State *s, const std::vector<double> &values) {
		ompl::base::SE3StateSpace::StateType *state = s->as<ompl::base::SE3StateSpace::StateType>();
		state->setXYZ(values[0], values[1], values[2]);
		double c1 = cos(values[5] / 2); //yaw
		double s1 = sin(values[5] / 2);
		double c2 = cos(values[4] / 2); //pitch 
		double s2 = sin(values[4] / 2);
		double c3 = cos(values[3] / 2); //roll
		double s3 = sin(values[3] / 2);
		double c1c2 = c1*c2;
		double s1s2 = s1*s2;
		state->rotation().w =c1c2 * c3 - s1s2 * s3;
		state->rotation().x =c1c2 * s3 + s1s2 * c3;
		state->rotation().y =s1 * c2 * c3 + c1 * s2 * s3;
		state->rotation().z =c1 * s2 * c3 - s1 * c2 * s3;
	};

	globalParameters.copyAbstractStateToVector = [](std::vector<double> &values, const ompl::base::State *s) {
		const ompl::base::SE3StateSpace::StateType *state = s->as<ompl::base::SE3StateSpace::StateType>();
		values.resize(6);
		values[0] = state->getX();
		values[1] = state->getY();
		values[2] = state->getZ();

		double x = state->rotation().x;
		double y = state->rotation().y;
		double z = state->rotation().z;
		double w = state->rotation().w;
		
		values[3] = atan2(2 * (x * y + w * z), w * w + x * x - y * y - z * z);
		values[4] = atan2(2 * (y * z + w * x), w * w - x * x - y * y + z * z);
		values[5] = asin(-2 * (x * z - w * y));
	};



	BenchmarkData data;
	data.benchmark = new ompl::tools::Benchmark(*quadrotorPtr, quadrotor->getName());
	data.simplesetup = quadrotorPtr;
	data.decomposition = quadrotor->allocDecomposition();

	return data;
}