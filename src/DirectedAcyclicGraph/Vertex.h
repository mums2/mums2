//
// Created by Gregory Johnson on 5/15/25.
//

#ifndef VERTEX_H
#define VERTEX_H
struct Vertex {
    size_t indexParentNode{}, indexChildNode{};
    double score{};
};

struct CompareVertexes {
    bool operator()(Vertex const& s1, Vertex const & s2) const {
        return s1.score > s2.score;
    }
};
#endif //VERTEX_H
