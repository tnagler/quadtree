#include <RcppEigen.h>
#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <unordered_map>
#include <functional>
#include <array>

// [[Rcpp::depends(RcppEigen)]]

struct Point
{
  double x, y;
  size_t index;
  bool operator==(const Point &other) const
  {
    return x == other.x && y == other.y && index == other.index;
  }
};

struct PointHash
{
  std::size_t operator()(const Point &p) const
  {
    return std::hash<size_t>()(p.index);
  }
  // An alternative would be to hash the position
  // std::size_t operator()(const Point &p) const
  // {
  //   std::size_t hx = std::hash<double>()(p.x);
  //   std::size_t hy = std::hash<double>()(p.y);
  //   return hx ^ (hy << 1); // Combining the two hashes
  // }
};

class PointSet
{
private:
  std::vector<Point> v;
  std::unordered_map<Point, size_t, PointHash> mp;

public:
  size_t size() { return v.size(); }

  bool search(const Point &point)
  {
    if (mp.find(point) != mp.end())
      return true;
    return false;
  }

  void insert(const Point &point)
  {

    if (search(point))
    {
      throw std::runtime_error("Cannot add point that's already in the set.");
    }

    if (v.size() == v.capacity())
    {
      v.reserve(v.size() * 2);
    }
    v.push_back(point);
    mp[point] = v.size() - 1;
  }

  void remove(const Point &point)
  {

    if (!search(point))
    {
      throw std::runtime_error("Cannot remove point that's not yet in the set.");
    }

    auto it = mp.find(point);
    v[it->second] = v.back();
    v.pop_back();
    mp[v[it->second]] = it->second;
    mp.erase(point);
  }

  const Point &sample(std::mt19937 *rng_ptr_)
  {
    auto n = v.size();
    if (n == 0)
    {
      throw std::runtime_error("Cannot sample from empty set.");
    }

    auto distribution = std::uniform_int_distribution<std::size_t>(0, n - 1);
    return v[distribution(*rng_ptr_) % v.size()];
    return v[rand() % v.size()];
  }
};

struct BoundingBox
{
  double x_min, y_min, x_mid, y_mid, x_max, y_max;

  BoundingBox(double x_min, double y_min, double x_max, double y_max)
      : x_min(x_min), y_min(y_min), x_max(x_max), y_max(y_max),
        x_mid(0.5 * (x_min + x_max)), y_mid(0.5 * (y_min + y_max))
  {
  }

  bool contains(const BoundingBox &other) const
  {
    return other.x_min >= x_min && other.x_max <= x_max &&
           other.y_min >= y_min && other.y_max <= y_max;
  }

  bool intersects(const BoundingBox &other) const
  {
    return !(other.x_min > x_max || other.x_max < x_min ||
             other.y_min > y_max || other.y_max < y_min);
  }
};

// Fixed depth quadtree class
class QuadTree2
{
private:
  struct Node
  {
    BoundingBox boundary;
    PointSet points;
    size_t point_count = 0;

    std::array<Node *, 4> children;
    Node(const BoundingBox &boundary) : boundary(boundary) {}
  };

  void construct_children(Node *node, uint16_t depth, int &node_idx)
  {
    if (depth >= depth_)
      return;

    const auto &b = node->boundary;
    nodes_.emplace_back(BoundingBox{b.x_min, b.y_min, b.x_mid, b.y_mid});
    nodes_.emplace_back(BoundingBox{b.x_mid, b.y_min, b.x_max, b.y_mid});
    nodes_.emplace_back(BoundingBox{b.x_min, b.y_mid, b.x_mid, b.y_max});
    nodes_.emplace_back(BoundingBox{b.x_mid, b.y_mid, b.x_max, b.y_max});

    for (int i = 0; i < 4; ++i)
    {
      node->children[i] = &nodes_[++node_idx];
    }

    for (auto &child : node->children)
    {
      construct_children(child, depth + 1, node_idx);
    }
  }

  Node *find_child_with_point(const Node *node, const Point &point)
  {
    if (point.x <= node->boundary.x_mid && point.y <= node->boundary.y_mid)
    {
      return node->children[0];
    }
    else if (point.y <= node->boundary.y_mid)
    {
      return node->children[1];
    }
    else if (point.x <= node->boundary.x_mid)
    {
      return node->children[2];
    }
    else
    {
      return node->children[3];
    }
  }

  std::vector<Node> nodes_;
  std::vector<Node *> terminal_nodes_;
  uint16_t depth_;
  std::mt19937 rng_;

public:
  QuadTree2(const BoundingBox &boundary, int depth, std::vector<int> seeds = {})
      : depth_(depth)
  {
    // geometric series sum formula 4^0 + 4^1 + ... + 4^depth
    size_t num_nodes = (std::pow(4, depth + 1) - 1) / 3;
    nodes_.reserve(num_nodes);
    terminal_nodes_.reserve(num_nodes);

    nodes_.emplace_back(boundary);
    int node_idx = 0;
    construct_children(&nodes_[0], 0, node_idx);

    auto seq = std::seed_seq(seeds.begin(), seeds.end());
    rng_ = std::mt19937(seq);
  }

  void insert(const Point &point)
  {
    Node *node = &nodes_[0]; // start at root
    node->point_count++;
    for (int d = 1; d <= depth_; ++d)
    {
      node = find_child_with_point(node, point);
      node->point_count++;
    }
    node->points.insert(point);
  }

  void remove(const Point &point)
  {
    Node *node = &nodes_[0]; // start at root
    node->point_count--;
    for (int d = 1; d <= depth_; ++d)
    {
      node = find_child_with_point(node, point);
      node->point_count--;
    }
    node->points.remove(point);
  }

  // Helper function to recursively calculate points within the range
  void find_terminal_nodes(const BoundingBox &range)
  {
    terminal_nodes_.clear(); // starting new search
    find_terminal_nodes_(&nodes_[0], range, 0);
  };

  // Helper function to recursively calculate points within the range
  void find_terminal_nodes_(Node *node, const BoundingBox &range, int depth)
  {
    // If the node is outside the range, no points are within the range
    if (!range.intersects(node->boundary))
    {
      return;
    }

    // If the node is fully contained in the range, all its points are
    if (depth == depth_ || range.contains(node->boundary))
    {
      terminal_nodes_.push_back(node);
      return;
    }

    // Otherwise, check each child node
    for (auto &child : node->children)
    {
      find_terminal_nodes_(child, range, depth + 1);
    }
  };

  Point sample(const BoundingBox &range)
  {
    // Find nodes at highest possible level fully contained in the range,
    // starting at root
    find_terminal_nodes(range);

    int total_points_in_range = 0;
    for (auto node : terminal_nodes_)
    {
      total_points_in_range += node->point_count;
    }
    if (total_points_in_range == 0)
    {
      throw std::runtime_error("No points in the specified range.");
    }

    std::uniform_int_distribution<int> dist(0, total_points_in_range - 1);
    int r = dist(rng_);

    Node *sampled_node;
    for (auto &node : terminal_nodes_)
    {
      if (r < node->point_count)
      {
        sampled_node = node;
        break;
      }
      r -= node->point_count;
    }

    while (sampled_node->points.size() == 0)
    {
      std::uniform_int_distribution<int> dist(0, sampled_node->point_count - 1);
      r = dist(rng_);

      for (auto &child : sampled_node->children)
      {
        if (r < child->point_count)
        {
          sampled_node = child;
          break;
        }
        r -= child->point_count;
      }
    }

    return sampled_node->points.sample(&rng_);
  };
};

// [[Rcpp::export]]
Eigen::MatrixXd test2(const Eigen::MatrixXd &points, const Eigen::MatrixXd &query,
                      int depth = 5)
{
  QuadTree2 quadtree(BoundingBox{0.0, 0.0, 1.0, 1.0}, depth);

  // Insert points from the Eigen matrix into the quadtree
  for (size_t i = 0; i < points.rows(); ++i)
  {
    quadtree.insert(Point{points(i, 0), points(i, 1), i});
  }

  // Define a query bounding box
  BoundingBox query_box{
      query(0, 0), query(0, 1),
      query(1, 0), query(1, 1)};

  Eigen::MatrixXd samples = points;

  for (size_t i = 0; i < samples.rows(); ++i)
  {

    Point old_point{samples(i, 0), samples(i, 1), i};
    quadtree.remove(old_point);

    Point sampled_point = quadtree.sample(query_box);
    // Point sampled_point = old_point;

    samples(i, 0) = sampled_point.x;
    samples(i, 1) = sampled_point.y;

    auto new_point = Point{0.2, 0.2, i};
    quadtree.insert(new_point);
  }

  return samples;
}
