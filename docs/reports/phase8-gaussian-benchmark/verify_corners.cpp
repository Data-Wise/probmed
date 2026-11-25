#include <algorithm>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

struct Bin {
    int start;
    std::int64_t weight;
};

bool covered(int u, int v, int z1, int z2,
             const std::vector<std::pair<int, int>>& offsets) {
    bool found_below = false;
    bool found_above = false;
    for (auto [dx, dy] : offsets) {
        bool below = true;
        bool above = true;
        for (int mu1 : {u, u + 10}) {
            for (int mu2 : {v, v + 10}) {
                for (int err1 : {z1, z1 + 10}) {
                    for (int err2 : {z2, z2 + 10}) {
                        int residual = (mu1 + err1 + dx) * (mu2 + err2 + dy) - mu1 * mu2;
                        below &= residual <= 0;
                        above &= residual >= 0;
                    }
                }
            }
        }
        found_below |= below;
        found_above |= above;
        if (found_below && found_above) return true;
    }
    return false;
}

int main() {
    int bin_count, offset_count;
    std::int64_t tail_weight;
    std::cin >> bin_count >> offset_count >> tail_weight;
    std::vector<Bin> bins(bin_count);
    for (auto& bin : bins) std::cin >> bin.start >> bin.weight;
    std::vector<std::pair<int, int>> offsets(offset_count);
    for (auto& offset : offsets) std::cin >> offset.first >> offset.second;
    std::int64_t largest = 0;
    int count = 0;
    for (int u = 0; u * u <= 20000; u += 10) {
        for (int v = 0; v <= u && u * u + v * v <= 20000; v += 10) {
            ++count;
            std::int64_t upper = tail_weight;
            for (const auto& first : bins) {
                for (const auto& second : bins) {
                    if (!covered(u, v, first.start, second.start, offsets)) {
                        upper += first.weight * second.weight;
                    }
                }
            }
            largest = std::max(largest, upper);
        }
    }
    std::cout << count << " parameter boxes, max numerator " << largest << '\n';
    return largest == 264990180241388LL && count == 92 ? 0 : 1;
}
