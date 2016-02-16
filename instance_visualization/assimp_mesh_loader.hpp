#pragma once

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <vector>

class AssimpMeshLoader {
public:
	struct Vertex {
		Vertex(double x, double y, double z) {
			vals.push_back(x);
			vals.push_back(y);
			vals.push_back(z);
		}

		double operator[](unsigned int index) const {
			return vals[index];
		}
		std::vector<double> vals;
	};

	struct Triangle {
		Triangle(unsigned int v0, unsigned int v1, unsigned int v2) {
			indices.push_back(v0);
			indices.push_back(v1);
			indices.push_back(v2);
		}

		unsigned int operator[](unsigned int index) const {
			return indices[index];
		}
		std::vector<unsigned int> indices;
	};


	AssimpMeshLoader(const std::string &pFile) : error(false) {
		scene = importer.ReadFile(pFile, aiProcess_FixInfacingNormals | aiProcess_GenNormals |  aiProcess_Triangulate);

		if(!scene) {
			fprintf(stderr, "%s\n", importer.GetErrorString());
			error = true;
		}
	}

	void get(std::vector< std::vector<Vertex> > &vertices,
	         std::vector< std::vector<Triangle> > &triangles,
	         std::vector< std::vector<double> > &normals) const {

		aiMatrix4x4 identity;
		walkSceneHierarchy(scene, scene->mRootNode, identity, vertices, triangles, normals);
	}

	void walkSceneHierarchy(const aiScene *scene, const aiNode *current, const aiMatrix4x4 &parentTransform,
	                        std::vector< std::vector<Vertex> > &vertices,
	                        std::vector< std::vector<Triangle> > &triangles,
	                        std::vector< std::vector<double> > &normals) const {

		const aiMatrix4x4 &transform = parentTransform * current->mTransformation;

		for(unsigned int j = 0; j < current->mNumMeshes; ++j) {
			assert(current->mMeshes[j] < scene->mNumMeshes);
			unsigned int meshIndex = current->mMeshes[j];

			auto mesh = scene->mMeshes[meshIndex];

			vertices.emplace_back();
			triangles.emplace_back();

			for(unsigned int j = 0; j < mesh->mNumVertices; ++j) {
				auto vertex = mesh->mVertices[j];
				double x = vertex.x;
				double y = vertex.y;
				double z = vertex.z;

				aiVector3D vec(x, y, z);
				vec *= transform;

				//apparently these are packed on linux and
				//we can't pass them in directly?
				x = vec.x;
				y = vec.y;
				z = vec.z;

				vertices.back().emplace_back(x,y,z);

				normals.emplace_back(3);
				if(mesh->HasNormals()) {
					aiVector3D norm(mesh->mNormals[j][0], mesh->mNormals[j][1], mesh->mNormals[j][2]);
					// norm *= transform;
					normals.back()[0] = norm.x;
					normals.back()[1] = norm.y;
					normals.back()[2] = norm.z;
				} else {
					normals.back()[0] = 1;
					normals.back()[1] = 0;
					normals.back()[2] = 0;
				}
			}

			for(unsigned int j = 0; j < mesh->mNumFaces; ++j) {
				auto face = mesh->mFaces[j];

				if(face.mNumIndices != 3) {
					fprintf(stderr, "face does not have 3 vertices: %d\n", face.mNumIndices);
					continue;
				}

				triangles.back().emplace_back(face.mIndices[0],
				                              face.mIndices[1],
				                              face.mIndices[2]);
			}
		}

		for(unsigned int i = 0; i < current->mNumChildren; ++i) {
			walkSceneHierarchy(scene, current->mChildren[i], transform, vertices, triangles, normals);
		}
	}

	//don't let the importer get destructed if you want to use scene (really stupid...)
	Assimp::Importer importer;
	const aiScene *scene;
	bool error;
};
