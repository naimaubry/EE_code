#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>  // For timing
#include <fstream>




// MyPoint class definition
class MyPoint : public std::array<double, 3> {
public:
    static const int DIM = 3;

    MyPoint() {}
    MyPoint(double x, double y, double z) {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
    }

    double distance(const MyPoint& other) const {
        double dist = 0.0;
        for (int i = 0; i < DIM; ++i) {
            dist += ((*this)[i] - other[i]) * ((*this)[i] - other[i]);
        }
        return std::sqrt(dist);
    }
};

// KDTree class definition
class KDTree {
private:
    struct Node {
        MyPoint point;
        Node* left;
        Node* right;

        Node(const MyPoint& pt) : point(pt), left(nullptr), right(nullptr) {}
    };

    Node* root;
    int depth;

    Node* build(std::vector<MyPoint>& points, size_t start, size_t end, int depth) {
        if (start >= end) return nullptr;
        int axis = depth % MyPoint::DIM;
        size_t mid = (start + end) / 2;
        std::nth_element(points.begin() + start, points.begin() + mid, points.begin() + end,
            [axis](const MyPoint& a, const MyPoint& b) {
                return a[axis] < b[axis];
            });

        Node* node = new Node(points[mid]);
        node->left = build(points, start, mid, depth + 1);
        node->right = build(points, mid + 1, end, depth + 1);
        return node;
    }

public:
    KDTree(const std::vector<MyPoint>& points) {
        std::vector<MyPoint> pts = points;
        root = build(pts, 0, pts.size(), 0);
    }

    MyPoint nearest(const MyPoint& query, Node* node, int depth, MyPoint& best) const {
        if (node == nullptr) return best; //!node

        if (query.distance(node->point) < query.distance(best)) {
            best = node->point;
        }

        int axis = depth % MyPoint::DIM;
        Node* nextNode = (query[axis] < node->point[axis]) ? node->left : node->right;
        Node* otherNode = (query[axis] < node->point[axis]) ? node->right : node->left;

        best = nearest(query, nextNode, depth + 1, best);
        if (std::abs(query[axis] - node->point[axis]) < query.distance(best)) {
            best = nearest(query, otherNode, depth + 1, best);
        }

        return best;
    }

    MyPoint nearest(const MyPoint& query) const {
        return nearest(query, root, 0, root->point);
    }
};

// Octree class definition
class Octree {
private:
    struct OctreeNode {
        MyPoint point;
        OctreeNode* children[8];

        OctreeNode(const MyPoint& pt) : point(pt) {
            std::fill(std::begin(children), std::end(children), nullptr);
        }
    };

    OctreeNode* root;
    MyPoint minBounds, maxBounds;

    int getOctant(const MyPoint& point, const MyPoint& mid) const {
        int octant = 0;
        if (point[0] >= mid[0]) octant |= 4;
        if (point[1] >= mid[1]) octant |= 2;
        if (point[2] >= mid[2]) octant |= 1;
        return octant;
    }

    void insert(OctreeNode*& node, const MyPoint& point, const MyPoint& min, const MyPoint& max) {
        if (!node) {
            node = new OctreeNode(point);
            return;
        }

        MyPoint mid;
        for (int i = 0; i < MyPoint::DIM; ++i) {
            mid[i] = (min[i] + max[i]) / 2.0;
        }
//Find out which octant the query point is in
        int octant = getOctant(point, mid);
        MyPoint newMin = min, newMax = max;
        for (int i = 0; i < MyPoint::DIM; ++i) {
            if (octant & (4 >> i)) {
                newMin[i] = mid[i];
            } else {
                newMax[i] = mid[i];
            }
        }

        insert(node->children[octant], point, newMin, newMax);
    }

    void nearest(const MyPoint& query, OctreeNode* node, const MyPoint& min, const MyPoint& max, MyPoint& best) const {
        if (!node) return;

        const double coef=1.1;
        if (query.distance(node->point) < query.distance(best)) {
            best = node->point;
        }

        MyPoint mid;
        for (int i = 0; i < MyPoint::DIM; ++i) {
            mid[i] = (min[i] + max[i]) / 2.0;
        }

        int octant = getOctant(query, mid);
        MyPoint newMin = min, newMax = max;
        for (int i = 0; i < MyPoint::DIM; ++i) {
            if (octant & (4 >> i)) {
                newMin[i] = mid[i];
            } else {
                newMax[i] = mid[i];
            }
        }

//look in the same octant
        nearest(query, node->children[octant], newMin, newMax, best);
        
        double bestDist=query.distance(best)*coef;
        
        MyPoint bestMin,bestMax;
        bestMin[0]=query[0]-bestDist;
        bestMin[1]=query[1]-bestDist;
        bestMin[2]=query[2]-bestDist;
        
        bestMax[0]=query[0]+bestDist;
        bestMax[1]=query[1]+bestDist;
        bestMax[2]=query[2]+bestDist;
        
//look at neighboring
        for (int i = 0; i < 8; ++i) {
            if (i != octant && node->children[i]) {
                MyPoint neighborMin = min, neighborMax = max; //childMin and childMax, we are checking children above for condition
                for (int j = 0; j < MyPoint::DIM; ++j) {
                    if (i & (4 >> j)) {
                        neighborMin[j] = mid[j];
                    } else {
                        neighborMax[j] = mid[j];
                    }
                }
                
                if(neighborMin[0]>bestMax[0])continue; //condition (if true then no way we have a better distance)
                if(neighborMin[1]>bestMax[1])continue;
                if(neighborMin[2]>bestMax[2])continue;
                if(neighborMax[0]<bestMin[0])continue;
                if(neighborMax[1]<bestMin[1])continue;
                if(neighborMax[2]<bestMin[2])continue;
                
                nearest(query, node->children[i], neighborMin, neighborMax, best);
                
                double bestDistloc = query.distance(best);
                if(bestDistloc<bestDist){
                    bestDist=bestDistloc;
                    bestMin[0]=query[0]-bestDist;
                    bestMin[1]=query[1]-bestDist;
                    bestMin[2]=query[2]-bestDist;
                    
                    bestMax[0]=query[0]+bestDist;
                    bestMax[1]=query[1]+bestDist;
                    bestMax[2]=query[2]+bestDist;
                }
            }
        }
    }

public:
    Octree(const std::vector<MyPoint>& points) {
        minBounds = MyPoint(0, 0, 0);
        maxBounds = MyPoint(1, 1, 1);

        root = nullptr;
        for (const auto& point : points) {
            insert(root, point, minBounds, maxBounds);
        }
    }

    MyPoint nearest(const MyPoint& query) const {
        MyPoint best = root->point;
        nearest(query, root, minBounds, maxBounds, best);
        return best;
    }
};

// UniformGrid class definition
class UniformGrid {
private:
    std::vector<std::vector<MyPoint>> grid;
    double cellSize;
    int gridSize;
    int gridSize2;
    std::vector<int> lstack,lmark;
    int imark;


    int getGridIndex(double value) const {
        return std::max(std::min(static_cast<int>(value / cellSize), gridSize - 1),0);
    }
    
    int getGridIndexGlobal(int ix,int iy,int iz)const{
        return gridSize2*ix+gridSize*iy+iz;
    }
    
    void getGridIndexLocal(int value,int&ix,int&iy,int&iz) const {
        ix=value/gridSize2;
        iy=(value-ix*gridSize2)/gridSize;
        iz=(value-ix*gridSize2-iy*gridSize);
        
        int val=getGridIndexGlobal(ix,iy,iz);
        if(val!=value){
            std::cout<<"perdu"<<std::endl;
        }
    }
    
    void getNeighbors(int ix,int iy,int iz,std::vector<int> &neighbors)const{
        neighbors.clear();
        if(ix>0){
            neighbors.push_back(getGridIndexGlobal(ix-1,iy,iz));
        }
        if(ix<gridSize-1){
            neighbors.push_back(getGridIndexGlobal(ix+1,iy,iz));
        }
        if(iy>0){
            neighbors.push_back(getGridIndexGlobal(ix,iy-1,iz));
        }
        if(iy<gridSize-1){
            neighbors.push_back(getGridIndexGlobal(ix,iy+1,iz));
        }
        if(iz>0){
            neighbors.push_back(getGridIndexGlobal(ix,iy,iz-1));
        }
        if(iz<gridSize-1){
            neighbors.push_back(getGridIndexGlobal(ix,iy,iz+1));
        }
    }

public:
    UniformGrid(const std::vector<MyPoint>& points)
         {
        gridSize=pow(points.size(),1.0/3.0)+10;
        gridSize2=gridSize*gridSize;
        cellSize = 1.0 / gridSize;
        grid.resize(gridSize * gridSize * gridSize);
        lmark.resize(grid.size(),-1);
        imark=0;

        for (const auto& point : points) {
            // convert these 3D indices into a single 1D index in order to store these cells in a 1D vector
            int idx = getGridIndexGlobal(getGridIndex(point[0]) , getGridIndex(point[1]) , getGridIndex(point[2]));
            grid[idx].push_back(point);
        }
    }

    MyPoint getFirstPoint(int xIdx,int yIdx,int zIdx){
        int idx=getGridIndexGlobal(xIdx, yIdx, zIdx);
        
        if(grid[idx].size()>0){
            return grid[idx][0];
        }
        
        int imark0=imark;
        imark++;
        lstack.clear();
        lstack.push_back(idx);
        lmark[idx]=imark0;
        int ix,iy,iz;
        std::vector<int> neighbors(6);
        for(size_t istack=0;istack<lstack.size();istack++){
            int iGrid=lstack[istack];
            getGridIndexLocal(iGrid, ix, iy, iz);
            getNeighbors(ix, iy, iz, neighbors);
            for(size_t i=0;i<neighbors.size();i++){
                int ineigh=neighbors[i];
                if(lmark[ineigh]!=imark0){
                    if(grid[ineigh].size()>0){
                        return grid[ineigh][0];
                    }
                    lmark[ineigh]=imark0;
                    lstack.push_back(ineigh);
                }
            }
        }

        return MyPoint();
    }
        
    MyPoint nearest(const MyPoint& query) {
        int xIdx = getGridIndex(query[0]);
        int yIdx = getGridIndex(query[1]);
        int zIdx = getGridIndex(query[2]);

        //Get a first close point
        MyPoint best=getFirstPoint(xIdx, yIdx, zIdx); // strategy is to get the first point fast
        double bestDist = query.distance(best);
        int nIdx=(bestDist/cellSize)+1;
        
        int nxmin=xIdx-nIdx;
        int nxmax=xIdx+nIdx+1;
        if(nxmin<0)nxmin=0;
        if(nxmax>gridSize)nxmax=gridSize;
        
        int nymin=yIdx-nIdx;
        int nymax=yIdx+nIdx+1;
        if(nymin<0)nymin=0;
        if(nymax>gridSize)nymax=gridSize;
        
        int nzmin=zIdx-nIdx;
        int nzmax=zIdx+nIdx+1;
        if(nzmin<0)nzmin=0;
        if(nzmax>gridSize)nzmax=gridSize;
        
        for (int nxIdx = nxmin; nxIdx < nxmax; ++nxIdx) {
            for (int nyIdx = nymin; nyIdx < nymax; ++nyIdx) {
                for (int nzIdx = nzmin; nzIdx <nzmax; ++nzIdx) {
                    
                    int idx = getGridIndexGlobal(nxIdx,nyIdx,nzIdx);
                    
                    for (const auto& point : grid[idx]) {
                        double dist = query.distance(point);
                        if (dist < bestDist) {
                            bestDist = dist;
                            best = point;
                        }
                    }
                }
            }
        }

        return best;
    }
};

// Main testing function
std::tuple<double, double, double> Compare(size_t npoint, bool uniform, bool error) {
    std::mt19937 mt;
    
    std::uniform_real_distribution<double> dist_uniform(0.0, 1.0);
    std::normal_distribution<double> dist_non_uniform(0.5, 0.15);  // Non-uniform distribution with mean 0.5 and std deviation 0.15


    std::vector<MyPoint> points(npoint);
        for (size_t i = 0; i < npoint; ++i) {
            if (uniform) {
                // Generate uniformly distributed points
                points[i][0] = dist_uniform(mt);
                points[i][1] = dist_uniform(mt);
                points[i][2] = dist_uniform(mt);
            } else {
                // Generate non-uniformly distributed points
                points[i][0] = std::clamp(dist_non_uniform(mt), 0.0, 1.0);  // Clamp to [0, 1] range
                points[i][1] = std::clamp(dist_non_uniform(mt), 0.0, 1.0);
                points[i][2] = std::clamp(dist_non_uniform(mt), 0.0, 1.0);
            }
        }
    
        std::vector<MyPoint> queries(npoint);
        for (size_t i = 0; i < npoint; ++i) {
            if (uniform) {
                queries[i][0] = dist_uniform(mt);
                queries[i][1] = dist_uniform(mt);
                queries[i][2] = dist_uniform(mt);
            } else {
                queries[i][0] = std::clamp(dist_non_uniform(mt), 0.0, 1.0);
                queries[i][1] = std::clamp(dist_non_uniform(mt), 0.0, 1.0);
                queries[i][2] = std::clamp(dist_non_uniform(mt), 0.0, 1.0);
            }
        }

    
    //Brute force
    std::vector<double> bestDist(npoint,0.0);
    std::vector<size_t> bestpoint(npoint,0);
    if(error){
        for (size_t i = 0; i < npoint; ++i) {
            bestDist[i]=queries[i].distance(points[0]);
            bestpoint[i]=0;
            for (size_t j = 1; j< npoint; ++j) {
                double distance=queries[i].distance(points[j]);
                if(distance<bestDist[i]){
                    bestDist[i]=distance;
                    bestpoint[i]=j;
                }
            }
        }
    }
    
    std::vector<double> KDError(npoint);
    std::vector<double> OctError(npoint);
    std::vector<double> UniformGridError(npoint);


    // Build KDTree
    auto start = std::chrono::high_resolution_clock::now();
    KDTree my_kdtree(points);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> kd_build_duration = end - start;
    std::cout << "KDTree build time: " << kd_build_duration.count() << " seconds\n";
    

    // Query KDTree
    start = std::chrono::high_resolution_clock::now();
    for (size_t i =0;i< queries.size();i++) {
        MyPoint nearest = my_kdtree.nearest(queries[i]);
//        std::cout << "KDTree nearest to (" << query[0] << ", " << query[1] << ", " << query[2]
//                  << ") is (" << nearest[0] << ", " << nearest[1] << ", " << nearest[2] << ")\n";
        double queryDistance = queries[i].distance(nearest);
        KDError[i] = abs(bestDist[i] - queryDistance);
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> kd_query_duration = end - start;
    std::cout << "KDTree query time: " << kd_query_duration.count() << " seconds\n";
    

    // Build Octree
    start = std::chrono::high_resolution_clock::now();
    Octree my_octree(points);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> octree_build_duration = end - start;
    std::cout << "Octree build time: " << octree_build_duration.count() << " seconds\n";

    // Query Octree
    start = std::chrono::high_resolution_clock::now();
    size_t count=0;
    for (size_t i =0;i< queries.size();i++) {
        count++;
        MyPoint nearest = my_octree.nearest(queries[i]);
//        std::cout << "Octree nearest to (" << query[0] << ", " << query[1] << ", " << query[2]
//                  << ") is (" << nearest[0] << ", " << nearest[1] << ", " << nearest[2] << ")\n";
        double queryDistance = queries[i].distance(nearest);
        OctError[i] = abs(bestDist[i] - queryDistance);
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> octree_query_duration = end - start;
    std::cout << "Octree query time: " << octree_query_duration.count() << " seconds\n";

    // Build Uniform Grid
    start = std::chrono::high_resolution_clock::now();
    UniformGrid my_grid(points);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> grid_build_duration = end - start;
    std::cout << "Uniform Grid build time: " << grid_build_duration.count() << " seconds\n";

    // Query Uniform Grid
    start = std::chrono::high_resolution_clock::now();
    for (size_t i =0;i< queries.size();i++) {
        MyPoint nearest = my_grid.nearest(queries[i]);
//        std::cout << "Uniform Grid nearest to (" << query[0] << ", " << query[1] << ", " << query[2]
//                  << ") is (" << nearest[0] << ", " << nearest[1] << ", " << nearest[2] << ")\n";
        double queryDistance = queries[i].distance(nearest);
        UniformGridError[i] = abs(bestDist[i] - queryDistance);
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> grid_query_duration = end - start;
    std::cout << "Uniform Grid query time: " << grid_query_duration.count() << " seconds\n";
    
    // Max errors and average
    std::cout << std::endl;
    
    if(error){
        double max1 = KDError[0];
        for (size_t i = 0; i < KDError.size(); i++)
            if (KDError[i] > max1) max1 = KDError[i];
        std::cout << "KD Tree Max Error " << max1 << std::endl;
        
        double sum1 = 0.0;
        for (size_t i = 0; i < KDError.size(); i++)
            sum1 += KDError[i];
        std::cout << "KD Tree Average Error " << sum1/KDError.size() << std::endl;
        
        double max2 = OctError[0];
        for (size_t i = 0; i < OctError.size(); i++)
            if (OctError[i] > max2) max2 = OctError[i];
        std::cout << "Oct Tree Max Error " << max2 << std::endl;
        
        double sum2 = 0.0;
        for (size_t i = 0; i < OctError.size(); i++)
            sum2 += OctError[i];
        std::cout << "Oct Tree Average Error " << sum2/OctError.size() << std::endl;
        
        double max3 = UniformGridError[0];
        size_t maxpos3=0;
        for (size_t i = 0; i < UniformGridError.size(); i++)
            if (UniformGridError[i] > max3) {
                max3 = UniformGridError[i];
                maxpos3=i;
            }
        std::cout << "Uniform Grid Max Error " << max3 << " maxpos:"<<maxpos3<<std::endl;
        
        double sum3 = 0.0;
        for (size_t i = 0; i < UniformGridError.size(); i++)
            sum3 += UniformGridError[i];
        std::cout << "Uniform Grid Average Error " << sum3/UniformGridError.size() << std::endl;
        
        std::cout << std::endl;
    }
    
    return std::make_tuple(kd_query_duration.count(), octree_query_duration.count(), grid_query_duration.count());
}

// For algorithm performance graphs
void save_results_to_file(const std::vector<size_t>& test, const std::vector<double>& kd_times,const std::vector<double>& octree_times, const std::vector<double>& grid_times, const std::string& filename) {
    
    std::ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << "Test-Size,KDTree-Time,Octree-Time,Uniform-Grid-Time\n";
        for (size_t i = 0; i < test.size(); ++i) {
            outfile << test[i] << "," << kd_times[i] << "," << octree_times[i] << "," << grid_times[i] << "\n";
        }
        outfile.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

// For plotting the uniform and non uniform distributions
void save_points_to_file(const std::vector<MyPoint>& points, const std::string& filename) {
    std::ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << "x,y,z\n";  // Write headers
        for (const auto& point : points) {
            outfile << point[0] << "," << point[1] << "," << point[2] << "\n";
        }
        outfile.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

// Function to generate and save points for uniform and non-uniform distributions
void Generate(size_t npoint, bool uniform, const std::string& filename) {
    std::mt19937 mt;
    std::uniform_real_distribution<double> dist_uniform(0.0, 1.0);
    std::normal_distribution<double> dist_non_uniform(0.5, 0.15);  // Non-uniform distribution

    std::vector<MyPoint> points(npoint);
    for (size_t i = 0; i < npoint; ++i) {
        if (uniform) {
            // Generate uniformly distributed points
            points[i][0] = dist_uniform(mt);
            points[i][1] = dist_uniform(mt);
            points[i][2] = dist_uniform(mt);
        } else {
            // Generate non-uniformly distributed points
            points[i][0] = std::clamp(dist_non_uniform(mt), 0.0, 1.0);
            points[i][1] = std::clamp(dist_non_uniform(mt), 0.0, 1.0);
            points[i][2] = std::clamp(dist_non_uniform(mt), 0.0, 1.0);
        }
    }

    // Save points to CSV file for visualization
    save_points_to_file(points, filename);
}


int main() {
    std::cout.precision(16);
    std::vector<size_t> test;
    
    // Uniform and non uniform distribution visualization
    size_t npoint = 1000; // number of points for the generate points
    Generate(npoint, true, "uniform_points.csv");
    Generate(npoint, false, "non_uniform_points.csv");
//    return 0;
    
    // Start from 1000, increment by 1.2x factor each time up to 1 million points
    for (size_t n = 1000; n <= 1000000; n = static_cast<size_t>(n * 1.2)) {
        test.push_back(n);
    }
    bool check_error=true;

    // Number of rounds for averaging
    const int num_rounds = 10;
    
    // Vectors to store averaged query times for uniform and non-uniform tests
    std::vector<double> uniform_kd_times, uniform_octree_times, uniform_grid_times;
    std::vector<double> non_uniform_kd_times, non_uniform_octree_times, non_uniform_grid_times;
    
    // Uniform test
    std::cout << "Starting uniform tests..." << std::endl;
    for (size_t size : test) {
        double total_kd_time = 0.0, total_octree_time = 0.0, total_grid_time = 0.0;
        for (int round = 0; round < num_rounds; ++round) {
            auto [kd_time, octree_time, grid_time] = Compare(size, true, check_error);
            total_kd_time += kd_time;
            total_octree_time += octree_time;
            total_grid_time += grid_time;
        }
        uniform_kd_times.push_back(total_kd_time / num_rounds);
        uniform_octree_times.push_back(total_octree_time / num_rounds);
        uniform_grid_times.push_back(total_grid_time / num_rounds);
    }
    
    // Non-uniform test
    std::cout << "Starting non-uniform tests..." << std::endl;
    for (size_t size : test) {
        double total_kd_time = 0.0, total_octree_time = 0.0, total_grid_time = 0.0;
        for (int round = 0; round < num_rounds; ++round) {
            auto [kd_time, octree_time, grid_time] = Compare(size, false, check_error);
            total_kd_time += kd_time;
            total_octree_time += octree_time;
            total_grid_time += grid_time;
        }
        non_uniform_kd_times.push_back(total_kd_time / num_rounds);
        non_uniform_octree_times.push_back(total_octree_time / num_rounds);
        non_uniform_grid_times.push_back(total_grid_time / num_rounds);
    }
    
    // Save uniform results to file
    save_results_to_file(test, uniform_kd_times, uniform_octree_times, uniform_grid_times, "uniform_test_results.csv");

    // Save non-uniform results to file
    save_results_to_file(test, non_uniform_kd_times, non_uniform_octree_times, non_uniform_grid_times, "non_uniform_test_results.csv");

    return 0;
}


