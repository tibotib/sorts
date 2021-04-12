#include <iostream>
#include <unistd.h>
#include <tuple>
#include <limits>
#include <vector>
#include <cmath>
#include <list>
#include <queue>
#include <algorithm>

namespace Util {
        template <typename T>
        T min(const T& a, const T& b) {
                return (a < b) ? a : b;
        }

        template <typename T>
        T max(const T& a, const T& b) {
                return (a > b) ? a : b;
        }

        template <typename T>
        void swap(T &a, T &b) {
                T tmp = a;
                a = b;
                b = tmp;
        }

        template <typename T>
        void display(T *a, int s) {
                for(int i = 0; i < s; i++)
                        std::cout << a[i] << " ";
                std::cout << std::endl;
        }

        template <typename T>
        void display(const std::vector<T> &a) {
                for(int i = 0; i < a.size(); i++)
                        std::cout << a[i] << " ";
                std::cout << std::endl;
        }

        template <typename T>
        void display(const std::vector<std::vector<T>> &a) {
                for(int i = 0; i < a.size(); i++)
                        Util::display(a[i]);
        }

        template <typename T>
        void generate(T *tab, int n) {
                srand(getpid());
                for(int i = 0; i < n; i++)
                        tab[i] = rand() % 100;
        }

        template <typename T>
        void generate(std::vector<T> &tab) {
                srand(getpid());
                for(int i = 0; i < tab.size(); i++)
                        tab[i] = rand() % 100;
        }

        template <typename T>
        T min(T *tab, int n) {
                T min = tab[0];
                for(int i = 1; i < n; i++) {
                        if(tab[i] < min)
                                min = tab[i];
                }
                return min;
        }

        template <typename T>
        T max(T *tab, int n) {
                T max = tab[0];
                for(int i = 1; i < n; i++) {
                        if(tab[i] > max)
                                max = tab[i];
                }
                return max;
        }

        template <typename T>
        std::pair<T, T> max_min(T *a, int n) {
                std::pair<T, T> ret(a[0], a[0]);
                for(int i = 1; i < n; i+=2) {
                        int max = i;
                        int min = i;
                        if(a[i] > a[i-1])
                                min = i-1;
                        if(a[max] > ret.first)
                                ret.first = a[max];
                        if(a[min] < ret.second)
                                ret.second = a[min];
                }
                if(n % 2 == 1) {
                        if(a[n-1] > ret.first)
                                ret.first = a[n-1];
                        else if(a[n-1] < ret.second)
                                ret.second = a[n-1];
                }
                return ret;
        }

        template <typename T>
        void display2D(const std::vector<T> &vec) {
                int n = sqrt(vec.size());
                for(int i = 0; i < n; i++) {
                        for(int j = 0; j < n; j++)
                                std::cout << vec[i*n+j] << " ";
                        std::cout << std::endl;
                }
        }

        template <typename T>
        void display(const std::vector<std::list<T>> &gr) {
                for(int i = 0; i < gr.size(); i++) {
                        std::cout << "Vertices " << i + 1 << " : " ;
                        for(auto e : gr[i])
                                std::cout << e << " ";
                        std::cout << std::endl;
                }
        }

        template <typename T>
        void display(const std::vector<std::pair<T, T>> &vec) {
                for(int i = 0; i < vec.size(); i++) {
                        std::cout << "NODE " << i+1 << " : " << vec[i].first << " and " << vec[i].second << std::endl;
                }
        }

        void display(const std::vector<std::list<std::pair<int, int>>> &gr) {
                for(int i = 0; i < gr.size(); i++) {
                        std::cout << "Vertices : " << i + 1 << " : " ;
                        for(auto e : gr[i])
                                std::cout << e.first << " value : " << e.second << std::endl;;
                        std::cout << std::endl;
                }
        }

        struct triplet{
                int x;
                int y;
                int z;
        };
        typedef struct triplet triplet;
};

#define LEFT(i) (i << 1)
#define RIGHT(i) ((i << 1) + 1)
#define PARENT(i) (i >> 1)

template <typename T>
class Heap {
public :
        Heap(const std::vector<T> &vec) : tab_(vec.size()) {
                std::copy(vec.begin(), vec.end(), this->tab_.begin());
                for(int i = (this->tab_.size() >> 1); i >= 0; i--) {
                        this->insert(i);
                }
        }

        Heap(T *tab, int s) : tab_(s){
                std::copy(tab, tab + s, this->tab_.begin());
                for(int i = (this->tab_.size() >> 1); i >= 0; i--) {
                        this->insert(i);
                }
        }

        void insert(T ele, int p) {//to insert element at the end and get it up
                if(p >= this->tab_.size())
                        this->tab_.push_back(p);
                else
                        this->tab_[p] = ele;

                auto p2 = PARENT(p);
                while(this->tab_[p] > this->tab_[p2] && p2 < p) {
                        Util::swap(this->tab_[p], this->tab_[p2]);
                        p = p2;
                        p2 = PARENT(p);
                }
        }

        void insert(int i) {//insert the element at index i
                if(i >= this->tab_.size() || i < 0)
                        return ;
                int l   = LEFT(i);
                int r   = RIGHT(i);
                int max = i;
                if(l < this->tab_.size() && this->tab_[i] < this->tab_[l])
                        max = l;
                if(r < this->tab_.size() && this->tab_[max] < this->tab_[r])
                        max = r;
                if(i == max)
                        return ;

                Util::swap(this->tab_[i], this->tab_[max]);
                this->insert(max);
        }

        T pop() {
                T ret = this->tab_[0];
                this->tab_[0] = this->tab_.back();
                this->tab_.pop_back();
                this->insert(0);
                return ret;
        }

        T max() {
                return this->tab_[0];
        }

        void display() {
                Util::display(this->tab_);
        }

private :
        std::vector<T> tab_;
};

namespace Sort {
        template <typename T>
        void merge(T *tab, int b, int m, int e) {
                int s1 = m - b + 1;
                int s2 = e - m;
                T *t1 = new T[s1];
                T *t2 = new T[s2];
                for(int i = 0; i < s1; i++)
                        t1[i] = tab[i + b];
                for(int i = 0; i < s2; i++)
                        t2[i] = tab[i + m + 1];

                int index_t1 = 0;
                int index_t2 = 0;
                int i        = 0;

                for(; index_t1 < s1 && index_t2 < s2; i++) {
                        if(t1[index_t1] < t2[index_t2]) {
                                tab[b + i] = t1[index_t1];
                                ++index_t1;
                                continue;
                        }
                        tab[b + i] = t2[index_t2];
                        ++index_t2;
                }

                for(; index_t1 < s1; index_t1++, i++)
                        tab[b + i] = t1[index_t1];
                for(; index_t2 < s2; index_t2++, i++)
                        tab[b + i] = t2[index_t2];

                delete t1;
                delete t2;
        }

        template <typename T>
        void merge_sort(T *tab, int f, int s) {//first call(tab, 0, tab.size() - 1)
                if(f < s) {
                        int mid = (f + s) >> 1;
                        merge_sort(tab, f, mid);
                        merge_sort(tab, mid + 1, s);
                        merge(tab, f, mid, s);
                }
        }

        template <typename T>
        void insertion_sort(T *a, int s) {
                for(int ii = 1; ii < s; ii++) {
                        int n = 1;
                        int i = ii;
                        while(a[i] < a[ii - n] && n <= ii) {
                                Util::swap(a[ii - n], a[i]);
                                --i;
                                ++n;
                        }
                }
        }

        template <typename T>
        void insertion_sort(std::vector<T> &a) {//for selection function
                for(int ii = 1; ii < a.size(); ii++) {
                        int n = 1;
                        int i = ii;
                        while(a[i] < a[ii - n] && n <= ii) {
                                Util::swap(a[ii - n], a[i]);
                                --i;
                                ++n;
                        }
                }
        }

        template <typename T>
        void bubble_sort(T *tab, int n) {
                for(int i = 0; i < n-1; i++) {
                        for(int j = 1; j < n-i; j++) {
                                if(tab[j-1] > tab[j])
                                        Util::swap(tab[j-1], tab[j]);
                        }
                }
        }

        template <typename T>
        void heap_sort(T *tab, int s) {
                Heap<T> hp = Heap<T>(tab, s);
                for(int i = s-1; i >= 0; i--) {
                        tab[i] = hp.pop();
                }
        }

        template <typename T>
        int partition(T *tab, int p, int r) {
                T pivot = tab[r];
                int n  = 0;
                for(int i = p; i < r; i++) {
                        if(tab[i] > pivot) {
                                ++n;
                        }
                        else{
                                Util::swap(tab[i - n], tab[i]);
                        }
                }
                Util::swap(tab[r], tab[r-n]);
                return r-n;
        }

        template <typename T>
        int random_partition(T *tab, int p, int r) {
                srand(getpid());
                auto i = p + (rand() % (p-r));
                Util::swap(tab[i], tab[r]);
                return partition(tab, p, r);
        }

        template <typename T>
        void quick_sort(T *tab, int b, int e) {//first call(tab, 0, tab.size() - 1)
                if(b < e) {
                        int mid = random_partition(tab, b, e);
                        quick_sort(tab, b, mid-1);
                        quick_sort(tab, mid+1, e);
                }
        }

        template <typename T = uint>
        void counting_sort(T *a, int size_a, int k) {//first call(tab, tab.size(), range of integers)
                uint *tab = new uint[k];
                for(int i = 0; i < size_a; i++) {
                        ++tab[a[i]];
                }
                int index_a = 0;
                for(int i = 0; i < k; i++) {
                        for(int j = 0; j < tab[i]; j++) {
                                a[index_a] = i;
                                ++index_a;
                        }
                }
                delete []tab;
        }
};
