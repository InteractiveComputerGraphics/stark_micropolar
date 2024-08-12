#include "include.h"

#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <random>

#define TINYOBJLOADER_IMPLEMENTATION
#include "tinyobjloader.h"

#include "unordered_array_set_and_map.h"
#include "mesh_utils.h"


// Tools
void push_back_if_not_present(std::vector<int>& v, const int value)
{
	if (std::find(v.begin(), v.end(), value) == v.end()) {
		v.push_back(value);
	}
}
bool is_outward_facing(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& tetCenter) {
	Eigen::Vector3d normal = (v1 - v0).cross(v2 - v0);
	Eigen::Vector3d faceCenter = (v0 + v1 + v2) / 3.0;
	Eigen::Vector3d toCenter = tetCenter - faceCenter;
	return normal.dot(toCenter) < 0; // outward if dot product is negative
}
// ========================================================================================================

double stark::deg2rad(const double deg)
{
	return 2.0 * M_PI * (deg / 360.0);
}
double stark::rad2deg(const double rad)
{
	return rad * 180.0 / M_PI;
}

std::vector<stark::Mesh<3>> stark::load_obj(const std::string& path)
{
    std::vector<stark::Mesh<3>> tri_meshes;

    // Initialise tinyobjloader objects and read the file
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn;
    std::string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, path.c_str());

    // Input checking
    if (!warn.empty()) {
        std::cout << warn << std::endl;
    }

    if (!err.empty()) {
        std::cerr << err << std::endl;
    }

    if (!ret) {
        exit(1);
    }

    // Global information
    const int total_n_vertices = (int)attrib.vertices.size() / 3;

    // Write the geometric information into individual triangular meshes
    // Loop over meshes
    for (int shape_i = 0; shape_i < shapes.size(); shape_i++) {

        // Initialize individual triangular mesh
        Mesh<3> tri_mesh;
        tri_mesh.conn.resize(shapes[shape_i].mesh.num_face_vertices.size());
        std::vector<bool> global_nodes_present(total_n_vertices, false);

        // Loop over triangles
        int index_offset = 0;
        for (int tri_i = 0; tri_i < shapes[shape_i].mesh.num_face_vertices.size(); tri_i++) {
            if (shapes[shape_i].mesh.num_face_vertices[tri_i] != 3) {
                std::cout << "learnSPH error: readTriMeshesFromObj can only read triangle meshes." << std::endl;
            }

            // Gather triangle global indices
            std::array<int, 3> triangle_global_indices;
            for (int vertex_i = 0; vertex_i < 3; vertex_i++) {
                tinyobj::index_t idx = shapes[shape_i].mesh.indices[(int)(3 * tri_i + vertex_i)];
                const int global_vertex_index = idx.vertex_index;
                triangle_global_indices[vertex_i] = global_vertex_index;
                global_nodes_present[global_vertex_index] = true;
            }
            tri_mesh.conn[tri_i] = triangle_global_indices;
        }

        // Reduce global indexes to local indexes
        std::vector<int> global_to_local_vertex_idx(total_n_vertices, -1);
        int local_vertices_count = 0;
        for (int global_vertex_i = 0; global_vertex_i < total_n_vertices; global_vertex_i++) {
            if (global_nodes_present[global_vertex_i]) {
                // Map global -> local
                global_to_local_vertex_idx[global_vertex_i] = local_vertices_count;
                local_vertices_count++;

                // Add vertex to the local mesh vertex vector
                tinyobj::real_t vx = attrib.vertices[(int)(3 * global_vertex_i + 0)];
                tinyobj::real_t vy = attrib.vertices[(int)(3 * global_vertex_i + 1)];
                tinyobj::real_t vz = attrib.vertices[(int)(3 * global_vertex_i + 2)];
                tri_mesh.vertices.push_back({ vx, vy, vz });
            }
        }

        // Change triangle indices
        for (int tri_i = 0; tri_i < tri_mesh.conn.size(); tri_i++) {
            for (int vertex_i = 0; vertex_i < 3; vertex_i++) {
                tri_mesh.conn[tri_i][vertex_i] = global_to_local_vertex_idx[tri_mesh.conn[tri_i][vertex_i]];
            }
        }

        tri_meshes.push_back(tri_mesh);
    }

	return tri_meshes;
}
void stark::write_VTK(const std::string& path, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 4>>& tets)
{
	vtkio::VTKFile vtk_file;

	if (vertices.size() == 0) {
		vtk_file.write_empty(path);
	}
	else {
		vtk_file.set_points_from_twice_indexable(vertices);
		vtk_file.set_cells_from_twice_indexable(tets, vtkio::CellType::Tetra);
		vtk_file.write(path);
	}
}
void stark::write_VTK(const std::string& path, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles, const bool generate_normals)
{
	vtkio::VTKFile vtk_file;

	if (vertices.size() == 0) {
		vtk_file.write_empty(path);
	}
	else {
		vtk_file.set_points_from_twice_indexable(vertices);
		vtk_file.set_cells_from_twice_indexable(triangles, vtkio::CellType::Triangle);
		if (generate_normals) {
			std::vector<Eigen::Vector3d> normals;
			compute_node_normals(normals, vertices, triangles);
			vtk_file.set_point_data_from_twice_indexable("normals", normals, vtkio::AttributeType::Vectors);
		}
		vtk_file.write(path);
	}
}
void stark::write_VTK(const std::string& path, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 2>>& edges)
{
	vtkio::VTKFile vtk_file;

	if (vertices.size() == 0) {
		vtk_file.write_empty(path);
	}
	else {
		vtk_file.set_points_from_twice_indexable(vertices);
		vtk_file.set_cells_from_twice_indexable(edges, vtkio::CellType::Line);
		vtk_file.write(path);
	}
}
void stark::write_VTK(const std::string& path, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 1>>& points)
{
	vtkio::VTKFile vtk_file;

	if (vertices.size() == 0) {
		vtk_file.write_empty(path);
	}
	else {
		vtk_file.set_points_from_twice_indexable(vertices);
		vtk_file.set_cells_as_particles(vertices.size());
		vtk_file.write(path);
	}
}

Eigen::Vector3d stark::triangle_normal(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2)
{
	return (p0 - p2).cross(p1 - p2).normalized();
}
double stark::triangle_area(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2)
{
	return 0.5 * (p0 - p2).cross(p1 - p2).norm();
}
double stark::unsigned_tetra_volume(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3)
{
	return (1.0 / 6.0) * ((p1 - p0).cross(p2 - p0)).dot(p3 - p0);
}
double stark::signed_tetra_volume(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3)
{
	return std::abs(unsigned_tetra_volume(p0, p1, p2, p3));
}

void stark::find_node_node_map_simplex(std::vector<std::vector<int>>& output, const int32_t* connectivity, const int32_t n_simplices, const int32_t n_nodes_per_simplex, const int32_t n_nodes)
{
	output.resize(n_nodes);
	for (int simplex_i = 0; simplex_i < n_simplices; simplex_i++) {
		const int begin = n_nodes_per_simplex * simplex_i;
		for (int i = 0; i < n_nodes_per_simplex; i++) {
			const int node_i = connectivity[begin + i];
		for (int j = i + 1; j < n_nodes_per_simplex; j++) {
			const int node_j = connectivity[begin + j];
			push_back_if_not_present(output[node_i], node_j);
			push_back_if_not_present(output[node_j], node_i);
		}
		}
	}
}
void stark::find_internal_angles(std::vector<std::array<int, 4>>& internal_angles, const std::vector<std::array<int, 3>>& triangles, const int n_nodes)
{
	/*
		For each edge in the triangle mesh, find the 2 nodes that are common neighbors of both edge end-points.
	*/
	internal_angles.clear();
	if (triangles.size() == 0) {
		return;
	}

	std::vector<std::vector<int>> node_node_map;
	find_node_node_map_simplex(node_node_map, &triangles[0][0], (int)triangles.size(), 3, n_nodes);
	for (std::vector<int>& nodes : node_node_map) {
		std::sort(nodes.begin(), nodes.end());  // std::set_intersection assumes sorted
	}

	std::vector<std::vector<int>> node_face_map;
	node_face_map.resize(n_nodes);
	for (int face_i = 0; face_i < (int)triangles.size(); face_i++) {
		const std::array<int, 3>& face = triangles[face_i];
		for (int node_i : face) {
			push_back_if_not_present(node_face_map[node_i], face_i);
		}
	}

	std::vector<std::array<int, 2>> edges;
	find_edges_from_simplices(edges, triangles, n_nodes);

	std::vector<int> buffer;
	std::vector<int> buffer2;
	internal_angles.reserve(edges.size());
	for (int edge_i = 0; edge_i < (int)edges.size(); edge_i++) {
		const std::array<int, 2>& edge = edges[edge_i];
		const std::vector<int>& neighs_i = node_node_map[edge[0]];
		const std::vector<int>& neighs_j = node_node_map[edge[1]];
		buffer.clear();
		// Check which neighbors are shared by the two nodes
		std::set_intersection(neighs_i.begin(), neighs_i.end(), neighs_j.begin(), neighs_j.end(), std::back_inserter(buffer));
		if (buffer.size() == 2) {
			// If the edge has exactly two shared neighbors, the internal angle is unambiguous...
			internal_angles.push_back({ edge[0], edge[1], buffer[0], buffer[1] });
		} else {
			// ...otherwise, we have to check the face connectivity of the shared neighbors
			buffer2.clear();
			// Check which shared neighbors are part of a face that also contains the edge
			for (const int i : buffer) {
				for (const int face_i : node_face_map[i]) {
					const auto& face = triangles[face_i];
					if (std::find(face.begin(), face.end(), edge[0]) != face.end() && std::find(face.begin(), face.end(), edge[1]) != face.end()) {
						buffer2.push_back(i);
					}
				}
			}
			if (buffer2.size() == 2) {
				internal_angles.push_back({ edge[0], edge[1], buffer2[0], buffer2[1] });
			} else {
				std::cout << "Stark error: triangle mesh has edges with more than two incident triangles." << std::endl;
				exit(-1);
			}
		}
	}
}
void stark::find_perimeter_edges(std::vector<std::array<int, 2>>& out_perimeter_edges, std::vector<int>& out_edge_to_triangle_node_map, const std::vector<std::array<int, 3>>& triangles, const int n_nodes)
{
	unordered_array_map<int, 2, int> edge_count;
	for (const auto& triangle : triangles) {
		edge_count[{ std::min(triangle[0], triangle[1]), std::max(triangle[0], triangle[1])} ]++;
		edge_count[{ std::min(triangle[1], triangle[2]), std::max(triangle[1], triangle[2])} ]++;
		edge_count[{ std::min(triangle[2], triangle[0]), std::max(triangle[2], triangle[0])} ]++;
	}

	std::vector<std::array<int, 2>> perimeter_edges;
	for (const auto& [edge, count] : edge_count) {
		if (count == 1) { // Edge appears only once, so it's a perimeter edge
			perimeter_edges.push_back(edge);
		}
	}

	reduce_connectivity(out_perimeter_edges, out_edge_to_triangle_node_map, perimeter_edges, n_nodes);
}
std::tuple<std::vector<std::array<int, 2>>, std::vector<int>> stark::find_perimeter_edges(const std::vector<std::array<int, 3>>& triangles, const int n_nodes)
{
	std::vector<std::array<int, 2>> edges;
	std::vector<int> edge_to_triangle_node_map;
	find_perimeter_edges(edges, edge_to_triangle_node_map, triangles, n_nodes);
	return std::make_tuple(edges, edge_to_triangle_node_map);
}
void stark::find_surface(std::vector<std::array<int, 3>>& out_triangles, std::vector<int>& out_triangle_to_tet_node_map, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 4>>& tets)
{
	out_triangles.clear();
	out_triangle_to_tet_node_map.clear();

	// Store the faces (with sorted connectivity) that occur only once and map them to their corresponding tet index
	unordered_array_map<int, 3, int> unique_face_tet_map;

	for (int tet_i = 0; tet_i < (int)tets.size(); tet_i++) {
		const std::array<int, 4>& tet = tets[tet_i];
		std::array<std::array<int, 3>, 4> faces = { {
			{{tet[0], tet[1], tet[2]}},
			{{tet[0], tet[1], tet[3]}},
			{{tet[0], tet[2], tet[3]}},
			{{tet[1], tet[2], tet[3]}}
		} };

		for (auto& face : faces) {
			std::sort(face.begin(), face.end());

			// If face exist, remove it from the list. Otherwise add it with the tet index.
			auto it = unique_face_tet_map.find(face);
			if (it == unique_face_tet_map.end()) {
				unique_face_tet_map[face] = tet_i;
			}
			else {
				unique_face_tet_map.erase(face);
			}
		}
	}

	// New connectivity + Face winding
	std::vector<std::array<int, 3>> unique_triangles; 
	unique_triangles.reserve(unique_face_tet_map.size());
	for (const auto& it : unique_face_tet_map) {
		std::array<int, 3> face = it.first;

		// Face winding to point outwards from the tet
		const std::array<int, 4>& tet = tets[it.second];
		const Eigen::Vector3d center = (vertices[tet[0]] + vertices[tet[1]] + vertices[tet[2]] + vertices[tet[3]]) / 4.0;
		if (is_outward_facing(vertices[face[0]], vertices[face[1]], vertices[face[2]], center)) {
			std::swap(face[0], face[1]);
		}

		unique_triangles.push_back(face);
	}

	// Reduce connectivity to a full, smaller mesh
	reduce_connectivity(out_triangles, out_triangle_to_tet_node_map, unique_triangles, (int)vertices.size());
}
std::tuple<std::vector<std::array<int, 3>>, std::vector<int>> stark::find_surface(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 4>>& tets)
{
	std::vector<std::array<int, 3>> out_triangles;
	std::vector<int> out_triangle_to_tet_node_map;
	find_surface(out_triangles, out_triangle_to_tet_node_map, vertices, tets);
	return { out_triangles, out_triangle_to_tet_node_map };
}
void stark::clean_triangle_mesh(std::vector<Eigen::Vector3d>& out_vertices, std::vector<std::array<int, 3>>& out_triangles, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles, const double merge_by_distance)
{
	// Quantize all vertices to a grid with cell size equal to merge_by_distance
	if (merge_by_distance > 0.0) {

		// AABB
		Eigen::AlignedBox3d bbox;
		for (const Eigen::Vector3d& v : vertices) {
			bbox.extend(v);
		}

		// Check dimensions
		double max_n_cells_double = bbox.diagonal().maxCoeff() / merge_by_distance;
		if (max_n_cells_double > (double)std::numeric_limits<int>::max()) {
			std::cout << "Stark error: clean_triangle_mesh merge_by_distance is too small." << std::endl;
			exit(-1);
		}

		// Quantize
		unordered_array_map<int, 3, int> unique_vertices;  // {cell_ijk: new_vertex_idx}
		std::vector<int> old_to_new_map(vertices.size(), -1);
		for (size_t old_i = 0; old_i < vertices.size(); ++old_i) {

			const std::array<int, 3> ijk = { 
				(int)((vertices[old_i].x() - bbox.min().x())/merge_by_distance),
				(int)((vertices[old_i].y() - bbox.min().y())/merge_by_distance),
				(int)((vertices[old_i].z() - bbox.min().z())/merge_by_distance)
			};

			int new_idx = -1;
			auto it = unique_vertices.find(ijk);
			if (it == unique_vertices.end()) {
				new_idx = out_vertices.size();
				out_vertices.push_back(vertices[old_i]);
				unique_vertices[ijk] = new_idx;
			}
			else {
				new_idx = it->second;
			}
			old_to_new_map[old_i] = new_idx;
		}

		// Update triangles
		out_triangles = apply_map(triangles, old_to_new_map);
	}
	else {
		out_vertices = vertices;
		out_triangles = triangles;
	}

	// Remove duplicated and degenerated triangles
	unordered_array_map<int, 3, int> unique_triangles;  //  {sorted_idx_triangle: old_triangle_idx}  We need this to conserve winding
	for (int tri_i = 0; tri_i < (int)out_triangles.size(); tri_i++) {
		std::array<int, 3> tri = out_triangles[tri_i];
		if (tri[0] != tri[1] && tri[0] != tri[2] && tri[1] != tri[2]) {
			std::sort(tri.begin(), tri.end());
			unique_triangles[tri] = tri_i;
		}
	}
	const std::vector<std::array<int, 3>> prev_triangles = out_triangles;
	out_triangles.clear();
	for (const auto& unique_tri : unique_triangles) {
		out_triangles.push_back(prev_triangles[unique_tri.second]);
	}
}
std::tuple<std::vector<Eigen::Vector3d>, std::vector<std::array<int, 3>>> stark::clean_triangle_mesh(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles, const double merge_by_distance)
{
	std::vector<Eigen::Vector3d> out_vertices;
	std::vector<std::array<int, 3>> out_triangles;
	clean_triangle_mesh(out_vertices, out_triangles, vertices, triangles, merge_by_distance);
	return { out_vertices, out_triangles };
}
void stark::find_sharp_edges(std::vector<std::array<int, 2>>& out_edges, std::vector<int>& out_old_to_new_map, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles, double angle_deg_threshold)
{
	std::vector<std::array<int, 2>> edges;
	std::vector<std::array<int, 4>> internal_angles;
	find_internal_angles(internal_angles, triangles, (int)vertices.size());

	const double cos_angle_rad_threshold = std::cos(deg2rad(angle_deg_threshold));
	edges.reserve(internal_angles.size());
	for (const std::array<int, 4>& angle : internal_angles) {
		const Eigen::Vector3d& p0 = vertices[angle[0]];
		const Eigen::Vector3d& p1 = vertices[angle[1]];
		const Eigen::Vector3d& p2 = vertices[angle[2]];
		const Eigen::Vector3d& p3 = vertices[angle[3]];

		const Eigen::Vector3d n0 = triangle_normal(p0, p1, p2);
		const Eigen::Vector3d n1 = triangle_normal(p1, p0, p3);
		if (n0.dot(n1) < cos_angle_rad_threshold) {
			edges.push_back({ angle[0], angle[1] });
		}
	}

	reduce_connectivity(out_edges, out_old_to_new_map, edges, (int)vertices.size());
}
std::tuple<std::vector<std::array<int, 2>>, std::vector<int>> stark::find_sharp_edges(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles, double angle_deg_threshold)
{
	std::vector<std::array<int, 2>> out_edges;
	std::vector<int> out_old_to_new_map;
	find_sharp_edges(out_edges, out_old_to_new_map, vertices, triangles, angle_deg_threshold);
	return { out_edges, out_old_to_new_map };
}


double stark::total_volume(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 4>>& tets)
{
	double volume = 0.0;
	for (const std::array<int, 4>&tet : tets) {
		volume += unsigned_tetra_volume(vertices[tet[0]], vertices[tet[1]], vertices[tet[2]], vertices[tet[3]]);
	}
	return volume;
}
void stark::compute_node_normals(std::vector<Eigen::Vector3d>& output, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles)
{
	output.resize(vertices.size(), Eigen::Vector3d::Zero());
	for (const std::array<int, 3> &triangle : triangles) {
		const Eigen::Vector3d normal = triangle_normal(vertices[triangle[0]], vertices[triangle[1]], vertices[triangle[2]]);
		const double area = triangle_area(vertices[triangle[0]], vertices[triangle[1]], vertices[triangle[2]]);
		output[triangle[0]] += area*normal;
		output[triangle[1]] += area*normal;
		output[triangle[2]] += area*normal;
	}
	for (Eigen::Vector3d& normal : output) {
		normal.normalize();
	}
}


void stark::move(std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& translation)
{
	for (Eigen::Vector3d& point : points) {
		point += translation;
	}
}
void stark::rotate_deg(std::vector<Eigen::Vector3d>& points, const double angle, const Eigen::Vector3d& axis)
{
	Eigen::Matrix3d R = Eigen::AngleAxis<double>(deg2rad(angle), axis.normalized()).toRotationMatrix();
	for (Eigen::Vector3d& point : points) {
		point = R * point;
	}
}
void stark::rotate_deg(std::vector<Eigen::Vector3d>& points, const double angle, const Eigen::Vector3d& axis, const Eigen::Vector3d& pivot)
{
	move(points, -pivot);
	rotate_deg(points, angle, axis);
	move(points, pivot);
}
void stark::scale(std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& scale)
{
	for (Eigen::Vector3d& point : points) {
		point = scale.cwiseProduct(point);
	}
}
void stark::scale(std::vector<Eigen::Vector3d>& points, const double s)
{
	scale(points, { s, s, s });
}
void stark::mirror(std::vector<Eigen::Vector3d>& points, const int dim, const double pivot)
{
	for (Eigen::Vector3d& point : points) {
		const double dist = point[dim] - pivot;
		point[dim] = pivot - dist;
	}
}

Eigen::Vector3d stark::rotate_deg(const Eigen::Vector3d& point, const Eigen::Matrix3d& R, const Eigen::Vector3d& pivot)
{
	const Eigen::Vector3d p_shifted = point - pivot;
	const Eigen::Vector3d p_shifted_rotated = R*p_shifted;
	const Eigen::Vector3d p_rotated = p_shifted_rotated + pivot;
	return p_rotated;
}

stark::Mesh<6> stark::tri3_to_tri6(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int32_t, 3>>& triangles)
{
    Mesh<6> mesh;

    auto& quadratic_vertices = mesh.vertices;
    auto& quadratic_triangles = mesh.conn;

    quadratic_vertices.reserve(vertices.size());
    quadratic_triangles.reserve(triangles.size());

    struct Edge {
        int32_t idx_a;
        int32_t idx_b;

        static Edge sorted(int32_t a, int32_t b) { return Edge { std::min(a, b),std::max(a, b) }; }
    };

    struct HashEdge
    {
        std::size_t operator()(const Edge& e) const
        {
            static_assert(sizeof(Edge) == sizeof(std::size_t));
            std::size_t h;
            std::memcpy(&h, &e, sizeof(std::size_t));
            return h;
        }
    };

    struct EqualEdge
    {
        bool operator()(const Edge& lhs, const Edge& rhs) const
        {
            return lhs.idx_a == rhs.idx_a && lhs.idx_b == rhs.idx_b;
        }
    };

    std::unordered_map<Edge, int, HashEdge, EqualEdge> edge_to_midpoint;

    // Copy all linear vertices into new storage for vertices of quadratic mesh
    quadratic_vertices.insert(quadratic_vertices.begin(), vertices.begin(), vertices.end());

    const int x = std::numeric_limits<int>::max();
    for (const auto& tri : triangles) {
        std::array<int, 6> new_triangle {x,x,x,x,x,x};
        for (int i = 0; i < 3; ++i) {
            const int vert_a_idx = tri[i];
            const int vert_b_idx = tri[(i + 1) % 3];
            const Edge edge = Edge::sorted(vert_a_idx, vert_b_idx);

            auto it = edge_to_midpoint.find(edge);
            int midpoint_idx = -1;
            if (it != edge_to_midpoint.end()) {
                midpoint_idx = it->second;
            } else {
                midpoint_idx = quadratic_vertices.size();
                edge_to_midpoint[edge] = midpoint_idx;

                const Eigen::Vector3d midpoint = (quadratic_vertices[vert_a_idx] + quadratic_vertices[vert_b_idx]) * 0.5;
                quadratic_vertices.push_back(midpoint);
            }

            new_triangle[i*2] = vert_a_idx;
            new_triangle[i*2 + 1] = midpoint_idx;
        }
        quadratic_triangles.push_back(new_triangle);
    }

    return mesh;
}
stark::Mesh<3> stark::tri6_to_tri3_coarsen(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int32_t, 6>>& triangles)
{
    Mesh<3> mesh;
    mesh.vertices = vertices;

    auto& linear_triangles = mesh.conn;

    for (int tri_i = 0; tri_i < triangles.size(); ++tri_i) {
        const auto& tri = triangles[tri_i];
        linear_triangles.push_back({tri[0], tri[2], tri[4]});
    }

    return mesh;
}
stark::Mesh<3> stark::tri6_to_tri3_subdivide(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int32_t, 6>>& triangles)
{
    Mesh<3> mesh;
    mesh.vertices = vertices;

    auto& linear_triangles = mesh.conn;

    for (int tri_i = 0; tri_i < triangles.size(); ++tri_i) {
        const auto& tri = triangles[tri_i];
        linear_triangles.push_back({tri[0], tri[1], tri[5]});
        linear_triangles.push_back({tri[1], tri[2], tri[3]});
        linear_triangles.push_back({tri[3], tri[4], tri[5]});
        linear_triangles.push_back({tri[1], tri[3], tri[5]});
    }

    return mesh;
}
stark::Mesh<3> stark::tri6_to_tri3_subdivide_n(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int32_t, 6>>& triangles, int subdivision_levels)
{
    Mesh<3> mesh_tri3;
    if (subdivision_levels == 0) {
        mesh_tri3 = tri6_to_tri3_coarsen(vertices, triangles);
    } else {
        mesh_tri3 = tri6_to_tri3_subdivide(vertices, triangles);
        for (int i = 1; i < subdivision_levels; ++i) {
            auto mesh_tri6 = tri3_to_tri6(mesh_tri3.vertices, mesh_tri3.conn);
            mesh_tri3 = tri6_to_tri3_subdivide(mesh_tri6.vertices, mesh_tri6.conn);
        }
    }
    return mesh_tri3;
}

stark::Mesh<10> stark::tri3_to_tri10(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int32_t, 3>>& triangles)
{
    Mesh<10> mesh;

    auto& cubic_vertices = mesh.vertices;
    auto& cubic_triangles = mesh.conn;

    cubic_vertices.reserve(vertices.size());
    cubic_triangles.reserve(triangles.size());

    struct Edge {
        int32_t idx_a;
        int32_t idx_b;

        static Edge sorted(int32_t a, int32_t b) { return Edge { std::min(a, b),std::max(a, b) }; }
    };

    struct HashEdge
    {
        std::size_t operator()(const Edge& e) const
        {
            static_assert(sizeof(Edge) == sizeof(std::size_t));
            std::size_t h;
            std::memcpy(&h, &e, sizeof(std::size_t));
            return h;
        }
    };

    struct EqualEdge
    {
        bool operator()(const Edge& lhs, const Edge& rhs) const
        {
            return lhs.idx_a == rhs.idx_a && lhs.idx_b == rhs.idx_b;
        }
    };

    struct InteriorVertices {
        int32_t idx_a_mid;  // Vertex closer to vertex a
        int32_t idx_b_mid;  // Vertex closer to vertex b
    };

    std::unordered_map<Edge, InteriorVertices, HashEdge, EqualEdge> edge_to_midpoint;

    // Copy all linear vertices into new storage for vertices of quadratic mesh
    cubic_vertices.insert(cubic_vertices.begin(), vertices.begin(), vertices.end());

    const int x = std::numeric_limits<int>::min();
    for (const auto& tri : triangles) {
        std::array<int, 10> new_triangle {x,x,x,x,x,x,x,x,x,x};
        for (int i = 0; i < 3; ++i) {
            const int vert_a_idx = tri[i];
            const int vert_b_idx = tri[(i + 1) % 3];
            const Edge edge = Edge::sorted(vert_a_idx, vert_b_idx);

            auto it = edge_to_midpoint.find(edge);
            int midpoint_a_idx = -1;
            int midpoint_b_idx = -1;
            if (it != edge_to_midpoint.end()) {
                const auto& interior = it->second;
                if (edge.idx_a == vert_a_idx) {
                    midpoint_a_idx = interior.idx_a_mid;
                    midpoint_b_idx = interior.idx_b_mid;
                } else {
                    midpoint_a_idx = interior.idx_b_mid;
                    midpoint_b_idx = interior.idx_a_mid;
                }
            } else {
                midpoint_a_idx = cubic_vertices.size();
                midpoint_b_idx = cubic_vertices.size() + 1;

                InteriorVertices interior;
                if (edge.idx_a == vert_a_idx) {
                    interior.idx_a_mid = midpoint_a_idx;
                    interior.idx_b_mid = midpoint_b_idx;
                } else {
                    interior.idx_a_mid = midpoint_b_idx;
                    interior.idx_b_mid = midpoint_a_idx;
                }
                edge_to_midpoint[edge] = interior;

                const Eigen::Vector3d midpoint_a = (2.0/3.0) * cubic_vertices[vert_a_idx] + (1.0/3.0) * cubic_vertices[vert_b_idx];
                const Eigen::Vector3d midpoint_b = (1.0/3.0) * cubic_vertices[vert_a_idx] + (2.0/3.0) * cubic_vertices[vert_b_idx];
                cubic_vertices.push_back(midpoint_a);
                cubic_vertices.push_back(midpoint_b);
            }

            new_triangle[i*3] = vert_a_idx;
            new_triangle[i*3 + 1] = midpoint_a_idx;
            new_triangle[i*3 + 2] = midpoint_b_idx;
        }

        // Add center vertex of triangle
        const auto& vert_a = vertices[tri[0]];
        const auto& vert_b = vertices[tri[1]];
        const auto& vert_c = vertices[tri[2]];
        const Eigen::Vector3d center = (1.0/3.0)*vert_a + (1.0/3.0)*vert_b + (1.0/3.0)*vert_c;
        const int center_idx = cubic_vertices.size();
        cubic_vertices.push_back(center);
        new_triangle[9] = center_idx;

        cubic_triangles.push_back(new_triangle);
    }

    return mesh;
}
stark::Mesh<3> stark::tri10_to_tri3_subdivide(const std::vector<Eigen::Vector3d>& vertices,
                                              const std::vector<std::array<int32_t, 10>>& triangles)
{
    Mesh<3> mesh;
    mesh.vertices = vertices;

    auto& linear_triangles = mesh.conn;
    linear_triangles.reserve(triangles.size() * 9);

    for (int tri_i = 0; tri_i < triangles.size(); ++tri_i) {
        const auto& tri = triangles[tri_i];
        linear_triangles.push_back({tri[0], tri[1], tri[8]});
        linear_triangles.push_back({tri[1], tri[9], tri[8]});
        linear_triangles.push_back({tri[1], tri[2], tri[9]});
        linear_triangles.push_back({tri[2], tri[4], tri[9]});
        linear_triangles.push_back({tri[2], tri[3], tri[4]});
        linear_triangles.push_back({tri[8], tri[9], tri[7]});
        linear_triangles.push_back({tri[9], tri[5], tri[7]});
        linear_triangles.push_back({tri[9], tri[4], tri[5]});
        linear_triangles.push_back({tri[7], tri[5], tri[6]});
    }

    return mesh;
}

template <std::size_t N_SUBDIV_TRIS, std::size_t N_SUBDIV_VERTS>
std::pair<std::vector<std::array<int32_t, 3>>, std::vector<stark::VertexLocalCoords>> tri_to_tri3_subdivide_with_local_coords(
    const int num_triangles,
    int subdivision_levels,
    const std::array<std::array<int, 3>, N_SUBDIV_TRIS>& subdiv_tris,
    const std::array<std::array<double, 3>, N_SUBDIV_VERTS>& subdiv_verts)
{
    using namespace stark;

    Mesh<3> mesh_tri3;
    std::vector<VertexLocalCoords> data;

    // Initialize triangle soup
    std::vector<std::array<int, 3>> triangle_soup;
    triangle_soup.reserve(num_triangles);

    std::vector<VertexLocalCoords> vertex_data;
    vertex_data.reserve(num_triangles * 3);

    const std::array<Eigen::Vector2d, 3> tri3_ref = {{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}};

    for (int tri_i = 0; tri_i < num_triangles; ++tri_i) {
        const int offset = vertex_data.size();

        triangle_soup.push_back({offset + 0, offset + 1, offset + 2});
        for (int i = 0; i < 3; ++i) {
            vertex_data.push_back(VertexLocalCoords{
                    tri_i,
                    tri3_ref[i]
            });
        }
    }

    if (subdivision_levels == 0) {
        return std::make_pair(triangle_soup, vertex_data);
    } else {
        const auto subdivide_once = [&subdiv_tris, &subdiv_verts](const std::vector<std::array<int, 3>>& triangle_soup, const std::vector<VertexLocalCoords>& data) {
            std::vector<std::array<int, 3>> new_triangles;
            new_triangles.reserve(triangle_soup.size() * N_SUBDIV_TRIS);

            std::vector<VertexLocalCoords> new_data;
            new_data.reserve(triangle_soup.size() * N_SUBDIV_TRIS * 3);

            for (int tri_i = 0; tri_i < triangle_soup.size(); ++tri_i) {
                const auto& tri = triangle_soup[tri_i];
                const auto& some_data = data[tri[0]];

                for (const auto& subdiv_tri : subdiv_tris) {
                    const int offset = new_data.size();
                    new_triangles.push_back({offset + 0, offset + 1, offset + 2});

                    for (int i = 0; i < 3; ++i) {
                        const int subdiv_vert = subdiv_tri[i];
                        const auto& c = subdiv_verts[subdiv_vert];

                        new_data.push_back(VertexLocalCoords{
                                some_data.origin_tri,
                                c[0] * data[tri[0]].local_coords + c[1] * data[tri[1]].local_coords + c[2] * data[tri[2]].local_coords,
                        });
                    }
                }
            }

            return std::make_pair(new_triangles, new_data);
        };

        for (int i = 0; i < subdivision_levels; ++i) {
            auto [new_triangle_soup, new_vertex_data] = subdivide_once(triangle_soup, vertex_data);
            triangle_soup = new_triangle_soup;
            vertex_data = new_vertex_data;
        }

        return std::make_pair(triangle_soup, vertex_data);
    }
}

std::pair<std::vector<std::array<int32_t, 3>>, std::vector<stark::VertexLocalCoords>> stark::tri6_to_tri3_subdivide_with_local_coords(const std::vector<std::array<int32_t, 6>>& triangles, int subdivision_levels)
{
    std::array<std::array<int, 3>, 4> subdiv_tris {{
        {0, 1, 5},
        {1, 2, 3},
        {5, 3, 4},
        {1, 3, 5},
    }};

    std::array<std::array<double, 3>, 6> linear_combinations {{
        {1.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.5, 0.5},
        {0.0, 0.0, 1.0},
        {0.5, 0.0, 0.5},
    }};

    return tri_to_tri3_subdivide_with_local_coords(triangles.size(), subdivision_levels, subdiv_tris, linear_combinations);
}
std::pair<std::vector<std::array<int32_t, 3>>, std::vector<stark::VertexLocalCoords>> stark::tri10_to_tri3_subdivide_with_local_coords(const std::vector<std::array<int32_t, 10>>& triangles, int subdivision_levels)
{
    std::array<std::array<int, 3>, 9> subdiv_tris {{
        {0, 1, 8},
        {1, 9, 8},
        {1, 2, 9},
        {2, 4, 9},
        {2, 3, 4},
        {8, 9, 7},
        {9, 5, 7},
        {9, 4, 5},
        {7, 5, 6},
    }};

    const double one_third = 1.0/3.0;
    const double two_third = 2.0/3.0;
    std::array<std::array<double, 3>, 10> linear_combinations {{
        {1.0, 0.0, 0.0},
        {two_third, one_third, 0.0},
        {one_third, two_third, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, two_third, one_third},
        {0.0, one_third, two_third},
        {0.0, 0.0, 1.0},
        {one_third, 0.0, two_third},
        {two_third, 0.0, one_third},
        {one_third, one_third, one_third},
    }};

    return tri_to_tri3_subdivide_with_local_coords(triangles.size(), subdivision_levels, subdiv_tris, linear_combinations);
}

stark::Mesh<3> stark::quad4_to_tri3(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int32_t, 4>>& quads)
{
    Mesh<3> mesh;
    mesh.vertices = vertices;

    auto& linear_triangles = mesh.conn;

    for (int quad_i = 0; quad_i < quads.size(); ++quad_i) {
        const auto& quad = quads[quad_i];
        linear_triangles.push_back({quad[0], quad[1], quad[2]});
        linear_triangles.push_back({quad[0], quad[2], quad[3]});
    }

    return mesh;
}
stark::Mesh<9> stark::quad4_to_quad9(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int32_t, 4>>& quads)
{
    Mesh<9> mesh;

    auto& quadratic_vertices = mesh.vertices;
    auto& quadratic_quads = mesh.conn;

    quadratic_vertices.reserve(vertices.size());
    quadratic_quads.reserve(quads.size());

    struct Edge {
        int32_t idx_a;
        int32_t idx_b;

        static Edge sorted(int32_t a, int32_t b) { return Edge { std::min(a, b),std::max(a, b) }; }
    };

    struct HashEdge
    {
        std::size_t operator()(const Edge& e) const
        {
            static_assert(sizeof(Edge) == sizeof(std::size_t));
            std::size_t h;
            std::memcpy(&h, &e, sizeof(std::size_t));
            return h;
        }
    };

    struct EqualEdge
    {
        bool operator()(const Edge& lhs, const Edge& rhs) const
        {
            return lhs.idx_a == rhs.idx_a && lhs.idx_b == rhs.idx_b;
        }
    };

    std::unordered_map<Edge, int, HashEdge, EqualEdge> edge_to_midpoint;

    // Copy all linear vertices into new storage for vertices of quadratic mesh
    quadratic_vertices.insert(quadratic_vertices.begin(), vertices.begin(), vertices.end());

    const int x = std::numeric_limits<int>::max();
    for (const auto& quad : quads) {
        std::array<int, 9> new_quad {x,x,x,x,x,x,x,x,x};
        for (int i = 0; i < 4; ++i) {
            const int vert_a_idx = quad[i];
            const int vert_b_idx = quad[(i + 1) % 4];
            const Edge edge = Edge::sorted(vert_a_idx, vert_b_idx);

            auto it = edge_to_midpoint.find(edge);
            int midpoint_idx = -1;
            if (it != edge_to_midpoint.end()) {
                midpoint_idx = it->second;
            } else {
                midpoint_idx = quadratic_vertices.size();
                edge_to_midpoint[edge] = midpoint_idx;

                const Eigen::Vector3d midpoint = 0.5 * (quadratic_vertices[vert_a_idx] + quadratic_vertices[vert_b_idx]);
                quadratic_vertices.push_back(midpoint);
            }

            new_quad[i*2] = vert_a_idx;
            new_quad[i*2 + 1] = midpoint_idx;
        }

        new_quad[8] = quadratic_vertices.size();
        quadratic_vertices.push_back(0.25 * (quadratic_vertices[quad[0]] + quadratic_vertices[quad[1]] + quadratic_vertices[quad[2]] + quadratic_vertices[quad[3]]));

        quadratic_quads.push_back(new_quad);
    }

    return mesh;
}
stark::Mesh<4> stark::quad9_to_quad4_subdivide(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int32_t, 9>>& quads)
{
    Mesh<4> mesh;
    mesh.vertices = vertices;

    auto& linear_quads = mesh.conn;
    linear_quads.reserve(quads.size() * 4);

    for (int quad_i = 0; quad_i < quads.size(); ++quad_i) {
        const auto& quad = quads[quad_i];
        linear_quads.push_back({quad[0], quad[1], quad[8], quad[7]});
        linear_quads.push_back({quad[1], quad[2], quad[3], quad[8]});
        linear_quads.push_back({quad[8], quad[3], quad[4], quad[5]});
        linear_quads.push_back({quad[7], quad[8], quad[5], quad[6]});
    }

    return mesh;
}
std::pair<std::vector<std::array<int32_t, 3>>, std::vector<stark::VertexLocalCoords>> quad_to_tri3_subdivide_with_local_coords(const int n_quads, int subdivision_levels)
{
    using namespace stark;

    Mesh<3> mesh_tri3;
    std::vector<VertexLocalCoords> data;

    // Initialize triangle soup
    std::vector<std::array<int, 3>> triangle_soup;
    triangle_soup.reserve(n_quads * 2);

    std::vector<VertexLocalCoords> vertex_data;
    vertex_data.reserve(n_quads * 6);

    const std::array<Eigen::Vector2d, 4> quad4_ref = {{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};

    for (int quad_i = 0; quad_i < n_quads; ++quad_i) {
        const int offset = vertex_data.size();

        triangle_soup.push_back({offset + 0, offset + 1, offset + 2});
        triangle_soup.push_back({offset + 3, offset + 4, offset + 5});
        vertex_data.push_back(VertexLocalCoords{ quad_i, quad4_ref[0] });
        vertex_data.push_back(VertexLocalCoords{ quad_i, quad4_ref[1] });
        vertex_data.push_back(VertexLocalCoords{ quad_i, quad4_ref[3] });
        vertex_data.push_back(VertexLocalCoords{ quad_i, quad4_ref[1] });
        vertex_data.push_back(VertexLocalCoords{ quad_i, quad4_ref[2] });
        vertex_data.push_back(VertexLocalCoords{ quad_i, quad4_ref[3] });
    }

    if (subdivision_levels == 0) {
        return std::make_pair(triangle_soup, vertex_data);
    } else {
        const auto subdivide_once = [](const std::vector<std::array<int, 3>>& triangle_soup, const std::vector<VertexLocalCoords>& data) {
            std::vector<std::array<int, 3>> new_triangles;
            new_triangles.reserve(triangle_soup.size() * 4);

            std::vector<VertexLocalCoords> new_data;
            new_data.reserve(triangle_soup.size() * 4 * 3);

            std::array<std::array<int, 3>, 4> subdiv_tris {{
                                                                   {0, 1, 5},
                                                                   {1, 3, 5},
                                                                   {1, 2, 3},
                                                                   {5, 3, 4},
                                                           }};

            std::array<std::array<double, 3>, 6> linear_combinations {{
                                                                              {1.0, 0.0, 0.0},
                                                                              {0.5, 0.5, 0.0},
                                                                              {0.0, 1.0, 0.0},
                                                                              {0.0, 0.5, 0.5},
                                                                              {0.0, 0.0, 1.0},
                                                                              {0.5, 0.0, 0.5},
                                                                      }};

            for (int tri_i = 0; tri_i < triangle_soup.size(); ++tri_i) {
                const auto& tri = triangle_soup[tri_i];
                const auto& some_data = data[tri[0]];

                for (const auto& subdiv_tri : subdiv_tris) {
                    const int offset = new_data.size();
                    new_triangles.push_back({offset + 0, offset + 1, offset + 2});

                    for (int i = 0; i < 3; ++i) {
                        const int subdiv_vert = subdiv_tri[i];
                        const auto& c = linear_combinations[subdiv_vert];

                        new_data.push_back(VertexLocalCoords{
                                some_data.origin_tri,
                                c[0] * data[tri[0]].local_coords + c[1] * data[tri[1]].local_coords + c[2] * data[tri[2]].local_coords,
                        });
                    }
                }
            }

            return std::make_pair(new_triangles, new_data);
        };

        for (int i = 0; i < subdivision_levels; ++i) {
            auto [new_triangle_soup, new_vertex_data] = subdivide_once(triangle_soup, vertex_data);
            triangle_soup = new_triangle_soup;
            vertex_data = new_vertex_data;
        }

        return std::make_pair(triangle_soup, vertex_data);
    }
}
std::pair<std::vector<std::array<int32_t, 3>>, std::vector<stark::VertexLocalCoords>> stark::quad4_to_tri3_subdivide_with_local_coords(const std::vector<std::array<int32_t, 4>>& quads, int subdivision_levels)
{
    return quad_to_tri3_subdivide_with_local_coords(quads.size(), subdivision_levels);
}
std::pair<std::vector<std::array<int32_t, 3>>, std::vector<stark::VertexLocalCoords>> stark::quad9_to_tri3_subdivide_with_local_coords(const std::vector<std::array<int32_t, 9>>& quads, int subdivision_levels)
{
    return quad_to_tri3_subdivide_with_local_coords(quads.size(), subdivision_levels);
}
