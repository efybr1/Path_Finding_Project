#include <iostream>
#include <vector>
#include <algorithm> // For std::make_heap, std::push_heap, std::pop_heap

class Wrapper {
private:
    int number; // Unique identifier
    int value; // Value for testing
    int index; // Pointer into vector

public:
    // Constructor
    Wrapper(int num, int val) : number(num), value(val), index(-1) {}

    // Getter and setter functions for index
    int getIndex() const {
        return index;
    }
    void setIndex(int idx) {
        index = idx;
    }

    // Getter functions for number and value
    int getNumber() const {
        return number;
    }
    int getValue() const {
        return value;
    }
};

auto compareByValue = [](const Wrapper& a, const Wrapper& b) {
    return a.getValue() > b.getValue();
};

void customPush(std::vector<Wrapper>& pq, Wrapper& wrapper) {
    wrapper.setIndex(pq.size());
    pq.push_back(wrapper);
    std::push_heap(pq.begin(), pq.end(), compareByValue);
}



int main() {
    std::vector<Wrapper> pq;

    Wrapper w1(1, 5);
    Wrapper w2(2, 4);
    Wrapper w3(3, 3);
    Wrapper w4(4, 2);
    Wrapper w5(5, 1);

    customPush(pq, w3);
    customPush(pq, w2);
    customPush(pq, w5);
    customPush(pq, w4);
    customPush(pq, w1);

    std::cout << "Before make_heap:" << std::endl;
    for (size_t i = 0; i < pq.size(); ++i) {
    const auto& wrapper = pq[i];
    std::cout << "Number: " << wrapper.getNumber() << ", Value: " << wrapper.getValue() << std::endl;
    }

    // Convert vector into a priority queue
    std::make_heap(pq.begin(), pq.end(), compareByValue);

    std::cout << "\nAfter make_heap:" << std::endl;

    while (!pq.empty()) {
        std::cout << "Number: " << pq.front().getNumber() << ", Value: " << pq.front().getValue() << std::endl;
        std::pop_heap(pq.begin(), pq.end(), compareByValue);
        pq.pop_back();
    }

    return 0;
}
