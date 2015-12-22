#include <unordered_set>
#include <unordered_map>

#include <ompl/base/samplers/UniformValidStateSampler.h>

#include <ompl/datastructures/NearestNeighborsSqrtApprox.h>
#include <ompl/datastructures/NearestNeighborsGNAT.h>
#include <ompl/datastructures/NearestNeighborsGNATNoThreadSafety.h>

#include "../domains/geometry/detail/FCLStateValidityChecker.hpp"

#include "../structs/probabilitydensityfunction.hpp"
#include "../structs/inplacebinaryheap.hpp"

namespace ompl {

namespace base {

class PlakuStateSampler : public ompl::base::UniformValidStateSampler {

protected:
	struct Vertex {
		Vertex(unsigned int id) : heapIndex(std::numeric_limits<unsigned int>::max()), id(id),
			numSelections(0), heuristic(std::numeric_limits<double>::infinity()), weight(0), onOpen(false) {}

		~Vertex() {}

		void setHeuristicAndPath(double heur, const std::vector<unsigned int> &path) {
			heuristic = heur;
			regionPath = path;
			regionPath.push_back(id);
		}

		void addPathCandidate(double heur, const std::vector<unsigned int> &path) {
			if(heur < heuristic) {
				bestIndex = heuristicCandidates.size();
				heuristic = heur;
				regionPath = path;
				regionPath.push_back(id);
			}
			heuristicCandidates.push_back(heur);
			regionPathCandidates.emplace_back(path.begin(), path.end());
			regionPathCandidates.back().emplace_back(id);
		}

		void currentInvalidFindNextBestPath() {
			regionPathCandidates.erase(regionPathCandidates.begin() + bestIndex);
			heuristicCandidates.erase(heuristicCandidates.begin() + bestIndex);

			if(heuristicCandidates.size() == 0) {
				bestIndex = 0;
				heuristic = std::numeric_limits<double>::infinity();
				regionPath.resize(0);
				return;
			}

			heuristic = std::numeric_limits<double>::infinity();
			for(unsigned int i = 0; i < heuristicCandidates.size(); ++i) {
				if(heuristicCandidates[i] < heuristic) {
					heuristic = heuristicCandidates[i];
					bestIndex = i;
				}
			}
			regionPath = regionPathCandidates[bestIndex];
		}

		void selected(double alpha) {
			numSelections++;
			weight = pow(alpha, numSelections) / (std::numeric_limits<double>::epsilon() + heuristic);
		}

		unsigned int getRandomRegionAlongPathToGoal(ompl::RNG &randomNumbers) const {
			unsigned int randomIndex = (unsigned int)(randomNumbers.uniform01() * regionPath.size());
			return regionPath[randomIndex];
		}

		unsigned int heapIndex;
		static bool pred(const Vertex *a, const Vertex *b) {
			return a->heuristic < b->heuristic;
		}
		static unsigned int getHeapIndex(const Vertex *r) {
			return r->heapIndex;
		}
		static void setHeapIndex(Vertex *r, unsigned int i) {
			r->heapIndex = i;
		}


		static bool HeapCompare(const Vertex *r1, const Vertex *r2) {
			return r1->weight < r2->weight;
		}

		unsigned int id, numSelections;
		double heuristic, weight;
		std::vector<unsigned int> regionPath;
		ompl::base::State *state;

		unsigned int bestIndex;
		std::vector<double> heuristicCandidates;
		std::vector<std::vector<unsigned int>> regionPathCandidates;
		bool onOpen;
	};

	struct Edge {
		Edge() {}

		Edge(unsigned int endpoint, double weight) : endpoint(endpoint), weight(weight) {

		}

		std::size_t operator()(const Edge &e) const {
			return e.endpoint;
		}

		bool operator==(const Edge &e) const {
			return e.endpoint == endpoint;
		}

		unsigned int endpoint;
		double weight;
	};

	double abstractDistanceFunction(const Vertex *a, const Vertex *b) const {
		assert(a->state != NULL);
		assert(b->state != NULL);
		return globalParameters.globalAbstractAppBaseGeometric->getStateSpace()->distance(a->state, b->state);
	}

	const base::State *stateIsAlreadyGeometric(const base::State *state, unsigned int /*index*/) const {
		return state;
	}

public:
	PlakuStateSampler(ompl::base::SpaceInformation *base, ompl::base::State *start, const ompl::base::GoalPtr &goal, double alpha, double b, double stateRadius)
		: UniformValidStateSampler(base), fullStateSampler(base->allocStateSampler()), alpha(alpha), b(b), stateRadius(stateRadius), activeRegion(NULL) {

		//Stolen from tools::SelfConfig::getDefaultNearestNeighbors
		if(base->getStateSpace()->isMetricSpace()) {
			// if (specs.multithreaded)
			//  nn.reset(new NearestNeighborsGNAT<Vertex*>());
			//else
			nn.reset(new NearestNeighborsGNATNoThreadSafety<Vertex *>());
		} else {
			nn.reset(new NearestNeighborsSqrtApprox<Vertex *>());
		}

		nn->setDistanceFunction(boost::bind(&PlakuStateSampler::abstractDistanceFunction, this, _1, _2));

		generateVertices(10000);
		generateEdges(10);

		unsigned int bestId = 0;
		double minDistance = std::numeric_limits<double>::infinity();
		double distance;
		for(auto vertex : vertices) {
			ompl::base::ScopedState<> vertexState(globalParameters.globalAppBaseControl->getGeometricComponentStateSpace());
			vertexState = vertex->state;
			ompl::base::ScopedState<> fullState = globalParameters.globalAppBaseControl->getFullStateFromGeometricComponent(vertexState);

			goal->isSatisfied(fullState.get(), &distance);
			if(distance <  minDistance) {
				bestId = vertex->id;
				minDistance = distance;
			}
		}

    Vertex startVertex(0);
    startVertex.state = start;
    startRegionId = nn->nearest(&startVertex)->id;
    
    goalRegionId = bestId;
		dijkstra(vertices[goalRegionId]);

		reached(start);
	}

	virtual ~PlakuStateSampler() {}

	std::vector<double> getColor(double min, double max, double value) const {
		std::vector<double> color(3);

		value = ((value - min) / (max - min)) * 765;

		if(value < 255) {
			color[0] = 0;
			color[1] = value / 2;
			color[2] = 255 - value;
		} else if(value < 510) {
			double relVal = value - 255;
			color[0] = relVal;
			color[1] = (relVal + 255) / 2;
			color[2] = 0;
		} else {
			double relVal = value - 510;
			color[0] = 255;
			color[1] = 255 - relVal;
			color[2] = 0;
		}

		for(unsigned int i = 0; i < 3; ++i) {
			color[i] /= 255;
		}

		return color;
	}

	void dumpToStderr() const {
		double min = std::numeric_limits<double>::infinity();
		double max = -std::numeric_limits<double>::infinity();
		for(unsigned int i = 0; i < vertices.size(); ++i) {
			if(vertices[i]->weight < min) min = vertices[i]->weight;
			if(vertices[i]->weight > max) max = vertices[i]->weight;
		}

		for(unsigned int i = 0; i < vertices.size(); ++i) {
			auto state = vertices[i]->state->as<ompl::base::SE3StateSpace::StateType>();
			auto color = getColor(min, max, vertices[i]->weight);

			fprintf(stderr, "point %g %g %g %g %g %g 1\n", state->getX(), state->getY(), state->getZ(), color[0], color[1], color[2]);
		}
	}

	void generatePythonPlotting() const {
		FILE *f = fopen("prm.py", "w");

		double min = std::numeric_limits<double>::infinity();
		double max = -std::numeric_limits<double>::infinity();
		for(unsigned int i = 0; i < vertices.size(); ++i) {
			if(vertices[i]->weight < min) min = vertices[i]->weight;
			if(vertices[i]->weight > max) max = vertices[i]->weight;
		}

		fprintf(f, "import numpy as np\nfrom mpl_toolkits.mplot3d import Axes3D\nimport matplotlib.pyplot as plt\n");
		fprintf(f, "fig = plt.figure()\nax = fig.add_subplot(111, projection='3d')\nax.scatter(");

		for(unsigned int coord = 0; coord < 4; ++coord) {
			if(coord == 3) {
				fprintf(f, "c=");
			}
			fprintf(f, "[");
			for(unsigned int i = 0; i < vertices.size(); ++i) {
				auto state = vertices[i]->state->as<ompl::base::SE3StateSpace::StateType>();
				if(coord == 0) {
					fprintf(f, (i < vertices.size()-1) ? "%g, ": "%g", state->getX());
				} else if(coord == 1) {
					fprintf(f, (i < vertices.size()-1) ? "%g, ": "%g", state->getY());
				} else if(coord == 2) {
					fprintf(f, (i < vertices.size()-1) ? "%g, ": "%g", state->getZ());
				} else if(coord == 3) {
					auto color = getColor(min, max, vertices[i]->weight);
					fprintf(f, (i < vertices.size()-1) ? "(%g, %g, %g), ": "(%g, %g, %g)", color[0], color[1], color[2]);
				}
			}
			fprintf(f, (coord < 3) ? "]," : "], depthshade=False");
		}
		fprintf(f, ")\nplt.show()\n");
		fclose(f);
	}

	virtual bool sample(ompl::base::State *state) {
		if(activeRegion != NULL) {
			if(!activeRegion->onOpen) {
				regionHeap.push_back(activeRegion);
				std::push_heap(regionHeap.begin(), regionHeap.end(), Vertex::HeapCompare);
				activeRegion->onOpen = true;
			}
		}

		if(randomNumbers.uniform01() < b) {
			assert(!regionHeap.empty());

			activeRegion = regionHeap.front();
			std::pop_heap(regionHeap.begin(), regionHeap.end(), Vertex::HeapCompare);
			regionHeap.pop_back();

			activeRegion->selected(alpha);
			activeRegion->onOpen = false;

			unsigned int regionAlongPath = activeRegion->getRandomRegionAlongPathToGoal(randomNumbers);

			ompl::base::ScopedState<> vertexState(globalParameters.globalAppBaseControl->getGeometricComponentStateSpace());
			vertexState = vertices[regionAlongPath]->state;

			ompl::base::ScopedState<> fullState = globalParameters.globalAppBaseControl->getFullStateFromGeometricComponent(vertexState);

			fullStateSampler->sampleUniformNear(state, fullState.get(), stateRadius);

			return true;
		} else {
			fullStateSampler->sampleUniform(state);
			return true;
		}
	}

	virtual bool sampleNear(ompl::base::State *, const ompl::base::State *, const double) {
		throw ompl::Exception("PlakuStateSampler::sampleNear", "not implemented");
		return false;
	}

	void reached(ompl::base::State *state) {
		ompl::base::ScopedState<> incomingState(si_->getStateSpace());
		incomingState = state;

		Vertex v(0);
		auto ss = globalParameters.globalAppBaseControl->getGeometricComponentState(incomingState, 0);
		v.state = ss.get();
		unsigned int newCellId = nn->nearest(&v)->id;

		if(!vertices[newCellId]->onOpen) {
			vertices[newCellId]->selected(alpha);
			regionHeap.push_back(vertices[newCellId]);
			std::push_heap(regionHeap.begin(), regionHeap.end(), Vertex::HeapCompare);
			vertices[newCellId]->onOpen = true;
		}
	}

protected:
	void generateVertices(unsigned int howMany) {
		ompl::base::StateSpacePtr abstractSpace = globalParameters.globalAbstractAppBaseGeometric->getStateSpace();
		ompl::base::ValidStateSamplerPtr abstractSampler = globalParameters.globalAbstractAppBaseGeometric->getSpaceInformation()->allocValidStateSampler();

		vertices.resize(howMany);

		for(unsigned int i = 0; i < howMany; ++i) {
			vertices[i] = new Vertex(i);
			vertices[i]->state = abstractSpace->allocState();
			abstractSampler->sample(vertices[i]->state);
			nn->add(vertices[i]);
		}
	}

	void generateEdges(unsigned int howManyConnections) {
		ompl::base::MotionValidatorPtr motionValidator = globalParameters.globalAbstractAppBaseGeometric->getSpaceInformation()->getMotionValidator();;

		auto distanceFunc = nn->getDistanceFunction();

		std::vector<Vertex *> neighbors;
		for(auto vertex : vertices) {
			neighbors.clear();
			nn->nearestK(vertex, howManyConnections+1, neighbors);

			assert(howManyConnections+1 >= neighbors.size());

			for(auto neighbor : neighbors) {

				if(edges.find(vertex->id) != edges.end() &&
				        edges[vertex->id].find(neighbor->id) != edges[vertex->id].end()) {
					continue;
				}

				double distance = distanceFunc(vertex, neighbor);
				if(distance == 0) continue;

				if(motionValidator->checkMotion(vertex->state, neighbor->state)) {
					edges[vertex->id][neighbor->id] = Edge(neighbor->id, distance);
					edges[neighbor->id][vertex->id] = Edge(vertex->id, distance);
				}
			}
		}
	}

	std::vector<unsigned int> getNeighboringCells(unsigned int index) const {
		auto vertexAndEdges = edges.find(index);
		assert(vertexAndEdges != edges.end());

		const auto &edges = vertexAndEdges->second;

		std::vector<unsigned int> ids;
		ids.reserve(edges.size());
		for(const auto &edge : edges) {
			ids.push_back(edge.second.endpoint);
		}

		return ids;
	}

	double getEdgeCostBetweenCells(unsigned int c1, unsigned int c2) const {
		auto vertexAndEdges = edges.find(c1);
		assert(vertexAndEdges != edges.end());

		const auto &edges = vertexAndEdges->second;

		auto edge = edges.find(c2);
		assert(edge != edges.end());

		return edge->second.weight;
	}

	void dijkstra(Vertex *start) {
		InPlaceBinaryHeap<Vertex, Vertex> open;
		std::unordered_set<unsigned int> closed;
		start->setHeuristicAndPath(0, std::vector<unsigned int>());
		open.push(start);

		closed.insert(start->id);

		while(!open.isEmpty()) {
			Vertex *current = open.pop();

			closed.insert(current->id);

			std::vector<unsigned int> kids = getNeighboringCells(current->id);
			for(unsigned int kid : kids) {
				if(closed.find(kid) != closed.end()) continue;

				double newHeuristic = current->heuristic + getEdgeCostBetweenCells(current->id, kid);
				Vertex *kidPtr = vertices[kid];

				if(newHeuristic < kidPtr->heuristic) {
					kidPtr->setHeuristicAndPath(newHeuristic, current->regionPath);

					if(open.inHeap(kidPtr)) {
						open.siftFromItem(kidPtr);
					} else {
						open.push(kidPtr);
					}
				}
			}
		}
	}

	StateSamplerPtr fullStateSampler;
	boost::shared_ptr< NearestNeighbors<Vertex *> > nn;
	std::vector<Vertex *> vertices, regionHeap;
	std::unordered_map<unsigned int, std::unordered_map<unsigned int, Edge>> edges;
	ompl::RNG randomNumbers;

	unsigned int startRegionId, goalRegionId;
	Vertex *activeRegion;

	double alpha, b, stateRadius;
};

}

}