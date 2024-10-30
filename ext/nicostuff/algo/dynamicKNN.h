#ifndef ___DYNAMIC_KNN__H__
#define ___DYNAMIC_KNN__H__

#include <ultimaille/all.h>
#include <iostream>




template<int D> struct DynamicKNN { // do not try anything but D=2 or D=3

    const std::vector<vec<D>>& pts;
    int n;
    std::vector<int> tree;

    DynamicKNN(const std::vector<vec<D>>& points) : pts(points), n(points.size()), tree(points.size()) {
        update_knn();
    }

    void build(const int L, const int R, const int dim = 0) { // build is O(n log n) using divide and conquer
        if (L >= R) return;
        int M = (L + R) / 2; // get median in O(n), split dim coordinate
        std::nth_element(tree.begin() + L, tree.begin() + M, tree.begin() + R, [this, dim](int a, int b) { return pts[a][dim] < pts[b][dim]; });
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp task
#endif
        build(L, M, (dim + 1) % D);
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp task
#endif
        build(M + 1, R, (dim + 1) % D);
    }


    void update_knn() {
        n = pts.size();
        if (n == 0) return;
        tree.resize(n);

        std::iota(tree.begin(), tree.end(), 0);
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel
#pragma omp single nowait
#endif
        build(0, n);
    }



    // k-nearest neighbor query, O(k log(k) log(n)) on average
    std::vector<int> query(const vec<D>& p, const int k = 1) {
        if (n + 1000 < pts.size()) update_knn();

        std::priority_queue<std::pair<double, int> > pq; // priority queue for KNN, keep the K nearest
        if (n > 0) {
            perform_query(pq, p, k, 0, n, false); // recursion
        }
        // START get extra points
        for (int i = n; i < pts.size(); i++) {
            vec3 I = pts[i];
            double d = (p - I).norm2();
            if (n > 0 && d > pq.top().first) continue;

            pq.push(std::make_pair(d, i));
            if (static_cast<int>(pq.size()) > k) pq.pop();
        }
        // END get extra points


        std::vector<int> neighbors;
        while (!pq.empty()) { // collect points
            neighbors.push_back(pq.top().second);
            pq.pop();
        }
        std::reverse(neighbors.begin(), neighbors.end());
        return neighbors;
    }

    void perform_query(std::priority_queue<std::pair<double, int> >& pq, const vec<D>& p, const int k, const int L, const int R, const int dim) const {
        if (L >= R) return;
        int M = (L + R) / 2;
        vec<D> d = p - pts[tree[M]];
        if (static_cast<int>(pq.size()) < k || d * d < pq.top().first) { // if point is nearer to the kth farthest, put point in queue
            pq.push(std::make_pair(d * d, tree[M]));
            if (static_cast<int>(pq.size()) > k) pq.pop(); // keep k elements only
        }
        int nearL = L, nearR = M, farL = M + 1, farR = R;
        if (d[dim] > 0) { // right is nearer
            std::swap(nearL, farL);
            std::swap(nearR, farR);
        }
        // query the nearer child
        perform_query(pq, p, k, nearL, nearR, (dim + 1) % D);

        if (static_cast<int>(pq.size()) < k || d[dim] * d[dim] < pq.top().first) // query the farther child if there might be candidates
            perform_query(pq, p, k, farL, farR, (dim + 1) % D);
    }
};


 inline void test_DynamicKNN() {
		PointSet pts;
		pts.create_points(1000);
		FOR(v, pts.size()) pts[v] = vec3(rand_range(0, 1), rand_range(0, 1), rand_range(0, 1));
        DynamicKNN<3> knn(*pts.data);

        FOR (iter,20){
            PointAttribute<int> attr(pts, 0);
            int off_v = pts.create_points(10000);
            for (int v = off_v; v < pts.size();v++) pts[v] = vec3(rand_range(0, 1), rand_range(0, 1), rand_range(0, 1));
            for (int v = off_v; v < pts.size(); v++)  attr[v] = 1;

            FOR(v,pts.size())                       knn.query(pts[v], 10);

            PointAttribute<int> nn(pts,0);
            auto res = knn.query(vec3(0, 0, 0), 10);
            for (int i : res) attr[i] +=2;
            Drop(pts, attr).apply("attr");
        }
    }

#endif