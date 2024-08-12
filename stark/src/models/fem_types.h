#pragma once

#include <fmt/format.h>

namespace stark {
namespace Fem {
    enum class ElementType : int {
        Tri3 = 0,
        Tri6 = 1,
        Tri10 = 2,
        Quad4 = 3,
        Quad9 = 4,
    };

    enum class ElementShape : int {
        Tri2d = 0,
        Quad2d = 1,
        Tet3d = 2,
    };

    static std::string element_type_to_string(const ElementType element_type);
    static ElementShape element_type_to_shape(const ElementType element_type);

    static std::string element_type_to_string(const ElementType element_type) {
        switch (element_type) {
            case ElementType::Tri3:
                return "Tri3";
            case ElementType::Tri6:
                return "Tri6";
            case ElementType::Tri10:
                return "Tri10";
            case ElementType::Quad4:
                return "Quad4";
            case ElementType::Quad9:
                return "Quad9";
            default:
                throw std::runtime_error(fmt::format("element_type_to_string: Unknown/unimplemented element type (int){}", static_cast<int>(element_type)));
        }
    }

    static ElementShape element_type_to_shape(const ElementType element_type) {
        switch (element_type) {
            case ElementType::Tri3:
            case ElementType::Tri6:
            case ElementType::Tri10:
                return ElementShape::Tri2d;
            case ElementType::Quad4:
            case ElementType::Quad9:
                return ElementShape::Quad2d;
            default:
                throw std::runtime_error(fmt::format("element_type_to_shape: Unknown/unimplemented element type {}", element_type_to_string(element_type)));
        }
    }
}
}
