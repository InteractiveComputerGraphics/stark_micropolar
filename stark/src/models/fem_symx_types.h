#pragma once

#include <array>

#include <symx>

#include "fem_quadrature.h"

namespace stark
{
namespace Fem
{
/// Connectivity for an FEM-based energy with quadrature
template<typename Element>
struct FemEnergyConnectivity {
    /// Connectivity: idx, group, DOFs, local quadrature index (for weights, points), global quadrature index (for quadrature point data)
    symx::LabelledConnectivity<2 + Element::num_nodes + 2> connectivity;

    Quadrature<Element::ref_dim> quadrature;
    std::vector<std::array<double, Element::ref_dim>> quadrature_points;
    std::vector<double> quadrature_weights;

    FemEnergyConnectivity() : connectivity({}) {
        // Initialize labels
        connectivity.labels.at(0) = "idx";
        connectivity.labels.at(1) = "group";
        for (int i = 0; i < Element::num_nodes; i++) {
            connectivity.labels.at(2 + i) = fmt::format("v{}", i);
        }
        connectivity.labels.at(2 + Element::num_nodes) = "q_idx_local";
        connectivity.labels.at(2 + Element::num_nodes + 1) = "q_idx_global";
    }

    FemEnergyConnectivity(ElementQuadrature<Element::ref_shape, Element::ref_dim> quad) : FemEnergyConnectivity() {
        set_quadrature(quad);
    }

    static std::vector<symx::Index> element_dof_slice(const symx::Element& ELE_WITH_QUAD) {
        return ELE_WITH_QUAD.slice(2, 2 + Element::num_nodes);
    };

    void set_quadrature(const Quadrature<Element::ref_dim>& quad) {
        if (this->connectivity.size() > 0) {
            throw std::runtime_error("FemEnergyConnectivity::set_quadrature: cannot safely change quadrature after elements have been added");
        }

        this->quadrature = quad;
        const int num_quad_points = this->quadrature.points.size();

        this->quadrature_points.resize(num_quad_points);
        for (int quad_i = 0; quad_i < num_quad_points; quad_i++) {
            const auto& p = this->quadrature.points[quad_i];
            for (int j = 0; j < Element::ref_dim; j++) {
                this->quadrature_points[quad_i][j] = p[j];
            }
        }
        quadrature_weights = quad.weights;
    }

    void push_element(const int idx, const int group, const std::array<int, Element::num_nodes>& dofs, const int q_idx_local, const int q_idx_global) {
        std::array<int, 4 + Element::num_nodes> conn;
        conn.at(0) = idx;
        conn.at(1) = group;
        for (int i = 0; i < Element::num_nodes; i++) {
            conn.at(2 + i) = dofs.at(i);
        }
        conn.at(2 + Element::num_nodes) = q_idx_local;
        conn.at(2 + Element::num_nodes + 1) = q_idx_global;
        connectivity.push_back(conn);
    }
};

/// Connectivity for an FEM-based energy with quadrature for problems with two DOF fields
template<typename PrimaryElement, typename SecondaryElement>
struct MixedFemEnergyConnectivity {
    static_assert(PrimaryElement::ref_dim == SecondaryElement::ref_dim, "Primary and secondary elements must have the same reference dimension");
    static_assert(PrimaryElement::ref_shape == SecondaryElement::ref_shape, "Primary and secondary elements must have the same reference shape");

    /// Connectivity: idx, group, DOFs, local quadrature index (for weights, points), global quadrature index (for quadrature point data)
    symx::LabelledConnectivity<2 + PrimaryElement::num_nodes + SecondaryElement::num_nodes + 2> connectivity;

    Quadrature<PrimaryElement::ref_dim> quadrature;
    std::vector<std::array<double, PrimaryElement::ref_dim>> quadrature_points;
    std::vector<double> quadrature_weights;

    MixedFemEnergyConnectivity() : connectivity({}) {
        // Initialize labels
        connectivity.labels.at(0) = "idx";
        connectivity.labels.at(1) = "group";
        for (int i = 0; i < PrimaryElement::num_nodes; i++) {
            connectivity.labels.at(2 + i) = fmt::format("v{}", i);
        }
        for (int i = 0; i < SecondaryElement::num_nodes; i++) {
            connectivity.labels.at(2 + PrimaryElement::num_nodes + i) = fmt::format("w{}", i);
        }
        connectivity.labels.at(2 + PrimaryElement::num_nodes + SecondaryElement::num_nodes + 0) = "q_idx_local";
        connectivity.labels.at(2 + PrimaryElement::num_nodes + SecondaryElement::num_nodes + 1) = "q_idx_global";
    }

    MixedFemEnergyConnectivity(ElementQuadrature<PrimaryElement::ref_shape, PrimaryElement::ref_dim> quad) : MixedFemEnergyConnectivity() {
        set_quadrature(quad);
    }

    static std::vector<symx::Index> element_primary_dof_slice(const symx::Element& ELE_WITH_QUAD) {
        return ELE_WITH_QUAD.slice(2, 2 + PrimaryElement::num_nodes);
    };

    static std::vector<symx::Index> element_secondary_dof_slice(const symx::Element& ELE_WITH_QUAD) {
        return ELE_WITH_QUAD.slice(2 + PrimaryElement::num_nodes, 2 + PrimaryElement::num_nodes + SecondaryElement::num_nodes);
    };

    void set_quadrature(const Quadrature<PrimaryElement::ref_dim>& quad) {
        if (this->connectivity.size() > 0) {
            throw std::runtime_error("MixedFemEnergyConnectivity::set_quadrature: cannot safely change quadrature after elements have been added");
        }

        this->quadrature = quad;
        const int num_quad_points = this->quadrature.points.size();

        this->quadrature_points.resize(num_quad_points);
        for (int quad_i = 0; quad_i < num_quad_points; quad_i++) {
            const auto& p = this->quadrature.points[quad_i];
            for (int j = 0; j < PrimaryElement::ref_dim; j++) {
                this->quadrature_points[quad_i][j] = p[j];
            }
        }
        quadrature_weights = quad.weights;
    }

    void push_element(const int idx, const int group, const std::array<int, PrimaryElement::num_nodes>& primary_dofs, const std::array<int, SecondaryElement::num_nodes>& secondary_dofs, const int q_idx_local, const int q_idx_global) {
        std::array<int, 4 + PrimaryElement::num_nodes + SecondaryElement::num_nodes> conn;
        conn.at(0) = idx;
        conn.at(1) = group;
        for (int i = 0; i < PrimaryElement::num_nodes; i++) {
            conn.at(2 + i) = primary_dofs.at(i);
        }
        for (int i = 0; i < SecondaryElement::num_nodes; i++) {
            conn.at(2 + PrimaryElement::num_nodes + i) = secondary_dofs.at(i);
        }
        conn.at(2 + PrimaryElement::num_nodes + SecondaryElement::num_nodes + 0) = q_idx_local;
        conn.at(2 + PrimaryElement::num_nodes + SecondaryElement::num_nodes + 1) = q_idx_global;
        connectivity.push_back(conn);
    }
};

struct FemSurfaceMeshes {
    /// Element type for each added mesh
    std::vector<Fem::ElementType> element_type;
    /// Index into element specific mesh vector for each added mesh
    std::vector<std::size_t> mesh_index;

    std::vector<std::vector<std::array<int, 3>>> mesh_tri3;
    std::vector<std::vector<std::array<int, 6>>> mesh_tri6;
    std::vector<std::vector<std::array<int, 10>>> mesh_tri10;
    std::vector<std::vector<std::array<int, 4>>> mesh_quad4;
    std::vector<std::vector<std::array<int, 9>>> mesh_quad9;

    template <typename ElementT>
    void add_mesh(const std::vector<std::array<int, ElementT::num_nodes>>& mesh) {
        if constexpr (ElementT::element_type == Fem::ElementType::Tri3) {
            add_tri3(mesh);
        } else if constexpr (ElementT::element_type == Fem::ElementType::Tri6) {
            add_tri6(mesh);
        } else if constexpr (ElementT::element_type == Fem::ElementType::Tri10) {
            add_tri10(mesh);
        } else if constexpr (ElementT::element_type == Fem::ElementType::Quad4) {
            add_quad4(mesh);
        } else if constexpr (ElementT::element_type == Fem::ElementType::Quad9) {
            add_quad9(mesh);
        } else {
            static_assert(!std::is_same_v<ElementT, ElementT>, "Unsupported element type");
        }
    }

    template <typename ElementT>
    std::vector<std::array<int, ElementT::num_nodes>>& get_mesh(int mesh_index) {
        if (this->element_type.at(mesh_index) != ElementT::element_type) {
            throw std::runtime_error("FemSurfaceMeshes::get_mesh: element type mismatch");
        }

        const int local_mesh_index = this->mesh_index.at(mesh_index);
        if constexpr(ElementT::element_type == Fem::ElementType::Tri3) {
            return this->mesh_tri3.at(local_mesh_index);
        } else if constexpr(ElementT::element_type == Fem::ElementType::Tri6) {
            return this->mesh_tri6.at(local_mesh_index);
        } else if constexpr(ElementT::element_type == Fem::ElementType::Tri10) {
            return this->mesh_tri10.at(local_mesh_index);
        } else if constexpr(ElementT::element_type == Fem::ElementType::Quad4) {
            return this->mesh_quad4.at(local_mesh_index);
        } else if constexpr(ElementT::element_type == Fem::ElementType::Quad9) {
            return this->mesh_quad9.at(local_mesh_index);
        } else {
            static_assert(!std::is_same_v<ElementT, ElementT>, "Unsupported element type");
        }
    }

    void add_tri3(const std::vector<std::array<int, 3>>& mesh) {
        this->element_type.push_back(Fem::ElementType::Tri3);
        this->mesh_index.push_back(this->mesh_tri3.size());
        this->mesh_tri3.push_back(mesh);
    }

    void add_tri6(const std::vector<std::array<int, 6>>& mesh) {
        this->element_type.push_back(Fem::ElementType::Tri6);
        this->mesh_index.push_back(this->mesh_tri6.size());
        this->mesh_tri6.push_back(mesh);
    }

    void add_tri10(const std::vector<std::array<int, 10>>& mesh) {
        this->element_type.push_back(Fem::ElementType::Tri10);
        this->mesh_index.push_back(this->mesh_tri10.size());
        this->mesh_tri10.push_back(mesh);
    }

    void add_quad4(const std::vector<std::array<int, 4>>& mesh) {
        this->element_type.push_back(Fem::ElementType::Quad4);
        this->mesh_index.push_back(this->mesh_quad4.size());
        this->mesh_quad4.push_back(mesh);
    }

    void add_quad9(const std::vector<std::array<int, 9>>& mesh) {
        this->element_type.push_back(Fem::ElementType::Quad9);
        this->mesh_index.push_back(this->mesh_quad9.size());
        this->mesh_quad9.push_back(mesh);
    }

    std::size_t size() const {
        return this->element_type.size();
    }

    int num_elements(int mesh_index) {
        int n_elements = -1;
        switch (this->element_type.at(mesh_index)) {
            case Fem::ElementType::Tri3: n_elements = this->mesh_tri3.at(this->mesh_index.at(mesh_index)).size(); break;
            case Fem::ElementType::Tri6: n_elements = this->mesh_tri6.at(this->mesh_index.at(mesh_index)).size(); break;
            case Fem::ElementType::Tri10: n_elements = this->mesh_tri10.at(this->mesh_index.at(mesh_index)).size(); break;
            case Fem::ElementType::Quad4: n_elements = this->mesh_quad4.at(this->mesh_index.at(mesh_index)).size(); break;
            case Fem::ElementType::Quad9: n_elements = this->mesh_quad9.at(this->mesh_index.at(mesh_index)).size(); break;
            default: throw std::runtime_error("Unknown element type");
        }
        return n_elements;
    }
};
}
}
