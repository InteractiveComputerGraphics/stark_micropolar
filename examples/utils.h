#pragma once

#include <optional>

/// Wraps std::optional and redefines * and -> operators to throw if no value is contained.
template<typename T>
struct Option : std::optional<T> {
    Option() = default;
    Option(std::nullopt_t value) : std::optional<T>(value) {}
    Option(std::optional<T> value) : std::optional<T>(value) {}
    Option(T value) : std::optional<T>(value) {}

    operator T() const { return this->value(); }

    T& operator*() & { return this->value(); }
    const T& operator*() const& { return this->value(); }
    T&& operator*() && { return this->value(); }
    const T&& operator*() const&& { return this->value(); }
    T* operator->() { return &this->value(); }
    const T* operator->() const { return &this->value(); }

    void apply_to(T& target) const {
        if (this->has_value()) {
            target = this->value();
        }
    }

    void apply_to(Option<T>& target) const {
        if (this->has_value()) {
            target = this->value();
        }
    }
};
