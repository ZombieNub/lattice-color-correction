#pragma once

#include <functional>
#include "ipm/ipx/basis.h"

struct EdgeInfo {
    int n1;
    int n2;

    EdgeInfo(const int in1, const int in2) {
        n1 = in1 < in2 ? in1 : in2;
        n2 = in1 < in2 ? in2 : in1;
    }

    bool operator==(const EdgeInfo &e) const {
        return (n1 == e.n1 && n2 == e.n2) || (n1 == e.n2 && n2 == e.n1);
    }

    bool operator<(const EdgeInfo &e) const {
        int ln1 = n1 < n2 ? n1 : n2;
        int ln2 = n1 < n2 ? n2 : n1;

        int rn1 = e.n1 < e.n2 ? e.n1 : e.n2;
        int rn2 = e.n1 < e.n2 ? e.n2 : e.n1;

        if (ln1 == rn1) {
            return ln2 < rn2;
        }
        return ln1 < rn1;
    }
};

template <>
struct std::hash<EdgeInfo> {
    std::size_t operator()(const EdgeInfo &e) const noexcept {
        const int n1 = e.n1 < e.n2 ? e.n1 : e.n2;
        const int n2 = e.n1 < e.n2 ? e.n2 : e.n1;
        const std::size_t h1 = std::hash<int>()(n1);
        const std::size_t h2 = std::hash<int>()(n2);
        return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};