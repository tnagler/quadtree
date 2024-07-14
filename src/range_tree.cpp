#include <RcppEigen.h>
#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <unordered_map>
#include <array>

// [[Rcpp::depends(RcppEigen)]]



struct Point {
  double x, y;
  size_t index;
  bool operator==(const Point& other) const {
    return x == other.x && y == other.y && index == other.index;
  }
};

struct PointID {
  Point* ptr;
  size_t index;
};

// A set of points with O(1) insert, remove, and uniform sampling;
// inspired by https://gist.github.com/matovitch/1c262d2d870024ca0b81
class PointSet
{
public:
  size_t size() { return vector_.size(); }

  void insert(const Point& p)
  {
    if (vector_.size() == vector_.capacity()) {
      vector_.reserve(vector_.size() * 2);
      set_.clear();
      for (auto& point : vector_) {
        set_.insert(PointID{&point, point.index});
      }
    }
    vector_.push_back(p);
    set_.insert(PointID{&vector_.back(), p.index});
  }

  void remove(const Point& p)
  {
    auto el = set_.lower_bound(PointID{nullptr, p.index});
    if (el->index != p.index) {
      throw std::runtime_error("Cannot remove point that's not yet in the tree.");
    }

    auto last = set_.lower_bound(PointID{nullptr, vector_.back().index});
    if (el->index != last->index) {
      auto new_id = PointID{el->ptr, last->index};
      *el->ptr = vector_.back();
      set_.erase(last);
      set_.insert(new_id);
    }

    vector_.pop_back();
    set_.erase(el);
  }

  const Point& sample(std::mt19937* rng_ptr_)
  {
    auto n = vector_.size();
    if (n == 0) {
      throw std::runtime_error("Cannot sample from empty set.");
    }

    auto distribution = std::uniform_int_distribution<std::size_t>(0, n - 1);
    return vector_[distribution(*rng_ptr_) % vector_.size()];
  }

private:
  struct IndexCompare {
    bool operator()(const PointID& a, const PointID& b) const {
      return a.index < b.index;
    }
  };

  std::set<PointID, IndexCompare> set_;
  std::vector<Point> vector_;
};


struct BoundingBox {
  double x_min, y_min, x_mid, y_mid, x_max, y_max;

  BoundingBox(double x_min, double y_min, double x_max, double y_max)
    : x_min(x_min), y_min(y_min), x_max(x_max), y_max(y_max),
      x_mid(0.5 * (x_min + x_max)), y_mid(0.5 * (y_min + y_max))
  {}

  bool contains(const BoundingBox& other) const {
    return other.x_min >= x_min && other.x_max <= x_max &&
      other.y_min >= y_min && other.y_max <= y_max;
  }

  bool intersects(const BoundingBox& other) const {
    return !(other.x_min > x_max || other.x_max < x_min ||
             other.y_min > y_max || other.y_max < y_min);
  }


  bool contains_x(const BoundingBox& other) const {
    return other.x_min >= x_min && other.x_max <= x_max;
  }

  bool contains_y(const BoundingBox& other) const {
    return other.y_min >= y_min && other.y_max <= y_max;
  }

  BoundingBox split_x(int right) const {
    return !right ?
    BoundingBox{x_min, y_min, x_mid, y_max} :
    BoundingBox{x_mid, y_min, x_max, y_max};
  }
  BoundingBox split_y(int top) const {
    return !top ?
    BoundingBox{x_min, y_min, x_max, y_mid} :
    BoundingBox{x_min, y_mid, x_max, y_max};
  }
};

struct Node {
  Node(int node_idx, const BoundingBox& boundary)
    : node_idx(node_idx), boundary(boundary) {}

  Node() {}

  BoundingBox boundary{0, 0, 0, 0};
  Node* up_child = nullptr;
  Node* down_child = nullptr;
  Node* left_child = nullptr;
  Node* right_child = nullptr;
  int node_idx = -1;

  int y_ = 0;

  PointSet points = {};
  size_t point_count = 0;
};


class RangeTree {
public:
  RangeTree(const BoundingBox& boundary, int depth, std::vector<int> seeds = {})
    : depth_(depth)
  {
    int num_nodes = std::pow(2, depth_ + 1) - 1;

    nodes_.resize(num_nodes);
    for (auto& n : nodes_) {
      n.resize(num_nodes);
    }
    nodes_[0][0] = Node{0, boundary};
    build_x_trees();
    link_y_trees();

    auto seq = std::seed_seq(seeds.begin(), seeds.end());
    rng_ = std::mt19937(seq);
  }

  void insert(const Point& point)
  {
    Node* node = &nodes_[0][0]; // start at x-root
    node->point_count++;
    for (int d = 1; d <= depth_; ++d) {
      node = find_child_with_point(node, point, 0);
      node->point_count++;
      auto y_node = node;
      for (int k = 1; k <= depth_; ++k) {
        y_node = find_child_with_point(y_node, point, 1);
        y_node->point_count++;
      }
      y_node->points.insert(point);
    }
  }

  void remove(const Point& point) {
    Node* node = &nodes_[0][0]; // start at x-root
    node->point_count--;
    for (int d = 1; d <= depth_; ++d) {
      node = find_child_with_point(node, point, 0);
      node->point_count--;
      auto y_node = node;
      for (int k = 1; k <= depth_; ++k) {
        y_node = find_child_with_point(y_node, point, 1);
        y_node->point_count--;
      }
      y_node->points.remove(point);
    }
  }

  // Helper function to recursively calculate points within the range
  void find_terminal_nodes(const BoundingBox& range)
  {
    terminal_nodes_.clear(); // starting new search
    find_terminal_nodes_(&nodes_[0][0], range);
  };

  // Helper function to recursively calculate points within the range
  void find_terminal_nodes_(Node* node, const BoundingBox& range)
  {
    // If the node is outside the range, no points are within the range
    if (!range.intersects(node->boundary)) {
      return;
    }

    // If the node is fully contained in the range, all its points are
    if (range.contains(node->boundary) || !node->left_child) {
      terminal_nodes_.push_back(node);
      return;
    }

    if (node->left_child) {
      find_terminal_nodes_(node->left_child, range);
      find_terminal_nodes_(node->right_child, range);
    }
  };

  Point sample(const BoundingBox& range) {
    // Find nodes at highest possible level fully contained in the range,
    // starting at root
    find_terminal_nodes(range);

    int total_points_in_range = 0;
    for (auto node : terminal_nodes_) {
      total_points_in_range += node->point_count;
    }
    if (terminal_nodes_.size() == 0) {
      throw std::runtime_error("No terminal nodes in the specified range.");
    }
    if (total_points_in_range == 0) {
      throw std::runtime_error("No points in the specified range.");
    }

    std::uniform_int_distribution<int> dist(0, total_points_in_range - 1);
    int r = dist(rng_);

    Node* sampled_node;
    for (auto& node : terminal_nodes_) {
      if (r < node->point_count) {
        sampled_node = node;
        break;
      }
      r -= node->point_count;
    }

    if (sampled_node->point_count == 0) {
      throw std::runtime_error("Not a terminal node");
    }


    while (sampled_node->points.size() == 0) {
      if (!sampled_node->down_child) {
        throw std::runtime_error("No children.");
      }

      dist = std::uniform_int_distribution<int>(0, sampled_node->point_count - 1);
      if (sampled_node->down_child) {
        if (dist(rng_) < sampled_node->down_child->point_count) {
          sampled_node = sampled_node->down_child;
        }  else {
          sampled_node = sampled_node->up_child;
        }
      }
    }

    return sampled_node->points.sample(&rng_);
  };


private:

  int x_ = 0;
  int y_ = 0;

  Node* find_child_with_point(const Node* node, const Point& point, int dim)
  {
    if (dim == 0) {
      return (point.x <= node->boundary.x_mid) ? node->left_child : node->right_child;
    } else {
      return (point.y <= node->boundary.y_mid) ? node->down_child : node->up_child;
    }
  }

  void build_x_trees(int node_idx = 0, int level = 1)
  {
    build_y_tree(node_idx, 0, 1);
    // if (!nodes_[node_idx][1].down_child) {
    //   std::cout << "failed ylink" << std::endl;
    // }
    if (level > depth_) return;

    // std::cout << "splitting x tree " << node_idx << " at level " << level << std::endl;

    int left_child_idx = ++x_;
    int right_child_idx = ++x_;

    // std::cout << left_child_idx << ", " <<right_child_idx << ", " << nodes_.size()
    //           << std::endl;

    auto node = &nodes_[node_idx][0];
    node->left_child  = &nodes_[left_child_idx][0];
    node->right_child = &nodes_[right_child_idx][0];
    *node->left_child  = Node{left_child_idx, node->boundary.split_x(false)};
    *node->right_child = Node{right_child_idx, node->boundary.split_x(true)};

    build_x_trees(left_child_idx, level + 1);
    build_x_trees(right_child_idx, level + 1);
  }

  void build_y_tree(int x_node_idx, int y_node_idx, int level) {
    if (level > depth_) return;
    // std::cout << "building " << x_node_idx << " at level " << level << std::endl;

    auto x_node = &nodes_[x_node_idx][0];
    int down_child_idx = ++x_node->y_;
    int up_child_idx = ++x_node->y_;


    // std::cout << "setting children " << x_node_idx << std::endl;

    auto node = &nodes_[x_node_idx][y_node_idx];
    nodes_[x_node_idx][down_child_idx]  = Node{x_node_idx, node->boundary.split_y(false)};
    nodes_[x_node_idx][up_child_idx] = Node{x_node_idx, node->boundary.split_y(true)};

    node->down_child = &nodes_[x_node_idx][down_child_idx];
    node->up_child  = &nodes_[x_node_idx][up_child_idx];

    // if (!nodes_[x_node_idx][y_node_idx].down_child) std::cout << "failed" << std::endl;
    build_y_tree(x_node_idx, down_child_idx, level + 1);
    build_y_tree(x_node_idx, up_child_idx, level + 1);
  }


  void link_y_trees(int node_idx = 0) {
    auto root_node = &nodes_[node_idx][0];
    if (!root_node->left_child) return;

    // std::cout << "linking " << node_idx << " -> " << root_node->left_child->node_idx <<
    //   ", " << root_node->right_child->node_idx << std::endl;
    auto& root_tree = nodes_[node_idx];
    auto& left_tree = nodes_[root_node->left_child->node_idx];
    auto& right_tree = nodes_[root_node->right_child->node_idx];

    for (int i = 0; i < root_tree.size(); ++i) {
      nodes_[node_idx][i].left_child = &nodes_[root_node->left_child->node_idx][i];
      nodes_[node_idx][i].right_child = &nodes_[root_node->right_child->node_idx][i];
    }

    // std::cout << "done" << std::endl;

    link_y_trees(root_node->left_child->node_idx);
    link_y_trees(root_node->right_child->node_idx);
  }


  int depth_;
  std::vector<std::vector<Node>> nodes_;
  std::vector<Node*> terminal_nodes_;
  std::mt19937 rng_;
};



// [[Rcpp::export]]
Eigen::MatrixXd test_range(const Eigen::MatrixXd& points, const Eigen::MatrixXd& query,
                           int depth = 5) {
  RangeTree rangetree(BoundingBox{0.0, 0.0, 1.0, 1.0}, depth);
  Eigen::MatrixXd samples = points;

  // Insert points from the Eigen matrix into the quadtree
  for (size_t i = 0; i < points.rows(); ++i) {
    rangetree.insert(Point{points(i, 0), points(i, 1), i});
  }

  // Define a query bounding box
  BoundingBox query_box{
    query(0, 0), query(0, 1),
    query(1, 0), query(1, 1)
  };


  for (size_t i = 0; i < samples.rows(); ++i) {

    Point old_point{samples(i, 0), samples(i, 1), i};
    rangetree.remove(old_point);

    Point sampled_point = rangetree.sample(query_box);
    // Point sampled_point = old_point;

    samples(i, 0) = sampled_point.x;
    samples(i, 1) = sampled_point.y;

    auto new_point = Point{0.2, 0.2, i};
    rangetree.insert(new_point);
  }

  return samples;
}


