#include <iostream>
#include <vector>
#include <initializer_list>
#include <cmath>
#include <limits>
#include <numbers>

const long double cEPS = 1e-9;
const long double cINF = std::numeric_limits<long double>::max();
const long double Pi = std::numbers::pi_v<long double>;


bool are_equal(long double a, long double b) {
  return std::abs(a - b) < cEPS;
}

bool less_or_equal(long double a, long double b) {
  return a < b || are_equal(a, b);
}

bool more_or_equal(long double a, long double b) {
  return a > b || are_equal(a, b);
}

struct Vector;
struct Line;

struct Point {
  long double x;
  long double y;

  Point(long double x = 0, long double y = 0);
  Point(const Point& other);

  Point& operator=(const Point& other);
  friend bool operator==(const Point& lhs, const Point& rhs);
  friend bool operator!=(const Point& lhs, const Point& rhs);
  Point& operator+=(const Vector& v);
  Point operator+(const Vector& v) const;
  Point& operator-=(const Vector& v);
  Point operator-(const Vector& v) const;

  void rotate(const Point& center, long double angle);
  void reflect(const Point& center);
  void reflect(const Line& axis);
  void scale(const Point& center, long double coefficient);
};

Point::Point(long double x, long double y)
  : x(x)
  , y(y)
  {}

Point::Point(const Point& other)
  : x(other.x)
  , y(other.y)
  {}

Point& Point::operator=(const Point& other) {
  Point tmp(other);
  std::swap(*this, tmp);
  return *this;
}

bool operator==(const Point& lhs, const Point& rhs) {
  return are_equal(lhs.x, rhs.x) && are_equal(lhs.y, rhs.y);
}

bool operator!=(const Point& lhs, const Point& rhs) {
  return !(lhs == rhs);
}


struct Vector {
  long double x;
  long double y;

  Vector();
  Vector(long double x, long double y);
  Vector(Point start, Point end);

  void reverse();
  Vector operator-() const;
  Vector& operator+=(const Vector& other);
  friend Vector operator+(const Vector& lhs, const Vector& rhs);
  Vector& operator-=(const Vector& other);
  friend Vector operator-(const Vector& lhs, const Vector& rhs);
  Vector& operator*=(const long double k);
  friend Vector operator*(const Vector& lhs, const long double k);
  friend Vector operator*(const long double k, const Vector& rhs);
  Vector& operator/=(const long double k);
  friend Vector operator/(const Vector& lhs, const long double k);
  friend long double operator*(const Vector& lhs, const Vector& rhs);
  friend long double operator,(const Vector& lhs, const Vector& rhs);
  long double length() const;
  void rotate(long double angle);
};

Vector::Vector()
  : Vector(0, 0)
  {}

Vector::Vector(long double x, long double y)
  : x(x)
  , y(y)
{}

Vector::Vector(Point start, Point end)
  : Vector(end.x - start.x, end.y - start.y)
{}

void Vector::reverse() {
  x = -x;
  y = -y;
}

Vector Vector::operator-() const {
  Vector tmp(*this);
  tmp.reverse();
  return tmp;
}

Vector& Vector::operator+=(const Vector& other) {
  x += other.x;
  y += other.y;
  return *this;
}

Vector operator+(const Vector& lhs, const Vector& rhs) {
  Vector tmp(lhs);
  tmp += rhs;
  return tmp;
}

Vector& Vector::operator-=(const Vector& other) {
  x -= other.x;
  y -= other.y;
  return *this;
}

Vector operator-(const Vector& lhs, const Vector& rhs) {
  Vector tmp(lhs);
  tmp -= rhs;
  return tmp;
}

Vector& Vector::operator*=(const long double k) {
  x *= k;
  y *= k;
  return *this;
}

Vector operator*(const Vector& lhs, const long double k) {
  Vector tmp(lhs);
  tmp *= k;
  return tmp;
}

Vector operator*(const long double k, const Vector& rhs) {
  Vector tmp(rhs);
  tmp *= k;
  return tmp;
}

Vector& Vector::operator/=(const long double k) {
  x /= k;
  y /= k;
  return *this;
}

Vector operator/(const Vector& lhs, const long double k) {
  Vector tmp(lhs);
  tmp /= k;
  return tmp;
}

long double operator*(const Vector& lhs, const Vector& rhs) {
  return lhs.x * rhs.y - rhs.x * lhs.y;
}

long double operator,(const Vector& lhs, const Vector& rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y;
};

long double Vector::length() const {
  return sqrtl(x * x + y * y);
}

void Vector::rotate(long double angle) {
  long double a_rad = angle * Pi / 180;
  long double old_x = x;
  x = old_x * cosl(a_rad) - y * sinl(a_rad);
  y = old_x * sinl(a_rad) + y * cosl(a_rad);
}


Point& Point::operator+=(const Vector& v) {
  x += v.x;
  y += v.y;
}

Point Point::operator+(const Vector& v) const {
  Point tmp(*this);
  tmp += v;
  return tmp;
}

Point& Point::operator-=(const Vector& v) {
  x -= v.x;
  y -= v.y;
  return *this;
}

Point Point::operator-(const Vector& v) const {
  Point tmp(*this);
  tmp -= v;
  return tmp;
}


class Line {
 private:
  long double a_;
  long double b_;
  long double c_;

 public:
  Line(const Point& first, const Point& second);
  Line(long double k, long double shift);
  Line(const Point& first, long double k);
  Line(const Vector& v);

  friend bool operator==(const Line& lhs, const Line& rhs);
  friend bool operator!=(const Line& lhs, const Line& rhs);
  friend struct Point;
};

Line::Line(const Point& first, const Point& second)
  : a_(second.y - first.y)
  , b_(first.x - second.x)
  , c_(second.x * first.y - first.x * second.y)
  {}

Line::Line(long double k, long double shift)
  : a_(k)
  , b_(-1)
  , c_(shift)
  {}

Line::Line(const Point& first, long double k)
  : a_(k)
  , b_(-1)
  , c_(-k * first.x + first.y)
  {}

Line::Line(const Vector& v)
  : Line(Point(0, 0), Point(v.x, v.y))
  {}

bool operator==(const Line& lhs, const Line& rhs) {
  return are_equal(lhs.a_ * rhs.b_, lhs.b_ * rhs.a_) &&
         are_equal(lhs.b_ * rhs.c_, lhs.c_ * rhs.b_) &&
         are_equal(lhs.c_ * rhs.a_, lhs.a_ * rhs.c_);
}

bool operator!=(const Line& lhs, const Line& rhs) {
  return !(lhs == rhs);
}


void Point::rotate(const Point& center, long double angle) {
  Vector v(center, *this);

  v.rotate(angle);
  *this = center + v;
}

void Point::reflect(const Point& center) {
  Vector v(*this, center);
  *this = center + v;
}

void Point::reflect(const Line& axis) {
  long double k = (axis.a_ * x + axis.b_ * y + axis.c_) / (axis.a_ * axis.a_ + axis.b_ * axis.b_);
  long double mx = x - k * axis.a_;
  long double my = y - k * axis.b_;
  Point M(mx, my);

  Vector v(*this, M);
  *this = M + v;
}

void Point::scale(const Point& center, long double coefficient) {
  Vector v(center, *this);
  *this = center + v * coefficient;
}


class Shape {
 public:
  virtual long double perimeter() const = 0;
  virtual long double area() const = 0;
  virtual bool operator==(const Shape& other) const = 0;
  virtual bool isCongruentTo(const Shape& other) const = 0;
  virtual bool isSimilarTo(const Shape& other) const = 0;
  virtual bool containsPoint(const Point& point) const = 0;

  virtual void rotate(const Point& center, long double angle) = 0;
  virtual void reflect(const Point& center) = 0;
  virtual void reflect(const Line& axis) = 0;
  virtual void scale(const Point& center, long double coefficient) = 0;
    virtual ~Shape() = default;
};


class Polygon: public Shape {
 protected:
  std::vector<Point> vertices_;

 public:
  Polygon(const std::vector<Point>& vertices);
  Polygon(std::initializer_list<Point> vertices);

  size_t verticesCount() const;
  const std::vector<Point>& getVertices() const;
  bool isConvex() const;

  long double perimeter() const;
  long double area() const;

  void rotate(const Point& center, long double angle);
  void reflect(const Point& center);
  void reflect(const Line& axis);
  void scale(const Point& center, long double coefficient);
  bool containsPoint(const Point& point) const;

  bool operator==(const Shape& other) const;
  bool isCongruentTo(const Shape& other) const;
  bool isSimilarTo(const Shape& other) const;
};

Polygon::Polygon(const std::vector<Point>& vertices)
  : vertices_(vertices)
  {}

Polygon::Polygon(std::initializer_list<Point> vertices)
  : vertices_(vertices)
{}

size_t Polygon::verticesCount() const {
  return vertices_.size();
}

const std::vector<Point>& Polygon::getVertices() const {
  return vertices_;
}

bool Polygon::isConvex() const {
  if (vertices_.size() < 3) {
    return 0;
  }

  int flag = 0;
  for (size_t i = 0; i < vertices_.size(); ++i) {
    Vector a(vertices_[i], vertices_[(i + 1) % vertices_.size()]);
    Vector b(vertices_[(i + 1) % vertices_.size()], vertices_[(i + 2) % vertices_.size()]);
    if (a * b < 0) {
      flag |= 1;
    } else {
      flag |= 2;
    }
    if (flag == 3) {
      return false;
    }
  }

  return flag != 3;
}

long double Polygon::perimeter() const {
  long double p = 0;
  for (int i = 0; i < vertices_.size(); ++i) {
    p += Vector(vertices_[i], vertices_[(i + 1) % vertices_.size()]).length();
  }

  return p;
}

long double Polygon::area() const {
  long double signed_area = 0;

  for (int i = 0; i < vertices_.size(); ++i) {
    const Point& p1 = vertices_[i];
    const Point& p2 = vertices_[(i + 1) % vertices_.size()];
    signed_area += (p1.x * p2.y - p2.x * p1.y);
  }
  signed_area /= 2.0;

  long double area = std::abs(signed_area);
  return area;
}

void Polygon::rotate(const Point& center, long double angle) {
  for (auto& p : vertices_) {
    p.rotate(center, angle);
  }
}

void Polygon::reflect(const Point& center) {
  for (auto& p : vertices_) {
    p.reflect(center);
  }
}

void Polygon::reflect(const Line& axis) {
  for (auto& p : vertices_) {
    p.reflect(axis);
  }
}

void Polygon::scale(const Point& center, long double coefficient) {
  for (auto& p : vertices_) {
    p.scale(center, coefficient);
  }
}

bool Polygon::containsPoint(const Point& point) const {
  for (size_t i = 0; i < vertices_.size(); ++i) {
    Vector edge_vec(vertices_[i], vertices_[(i + 1) % vertices_.size()]);
    Vector point_vec(vertices_[i], point);
    if (are_equal(edge_vec * point_vec, 0)) {
      long double min_x = std::min(vertices_[i].x, vertices_[(i + 1) % vertices_.size()].x);
      long double min_y = std::min(vertices_[i].y, vertices_[(i + 1) % vertices_.size()].y);
      long double max_x = std::max(vertices_[i].x, vertices_[(i + 1) % vertices_.size()].x);
      long double max_y = std::max(vertices_[i].y, vertices_[(i + 1) % vertices_.size()].y);

      if (less_or_equal(min_x, point.x) && less_or_equal(point.x, max_x) &&
          less_or_equal(min_y, point.y) && less_or_equal(point.y, max_y)) {
        return true;
      }
    }
  }

  int intersections = 0;
  for (size_t i = 0; i < vertices_.size(); ++i) {
    Point p1(vertices_[i]);
    Point p2(vertices_[(i + 1) % vertices_.size()]);

    if (are_equal(p1.y, p2.y)) {
      continue;
    }

    if ((less_or_equal(p1.y, point.y) && p2.y > point.y) || (less_or_equal(p2.y, point.y) && p1.y > point.y)) {
      long double x_intersect = p1.x + (point.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);

      if (x_intersect > point.x) {
        ++intersections;
      }
    }
  }

  return intersections % 2;
}

bool point_cmp(const Point& a, const Point& b) {
  if (!are_equal(a.x, b.x)) return a.x < b.x;
  return a.y < b.y;
}

bool Polygon::operator==(const Shape& other) const {
  if (const Polygon* other_polygon = dynamic_cast<const Polygon*>(&other)) {
    if (this->verticesCount() != other_polygon->verticesCount()) {
      return false;
    }

    auto vertices1 = this->getVertices();
    auto vertices2 = other_polygon->getVertices();

    std::sort(vertices1.begin(), vertices1.end(), point_cmp);
    std::sort(vertices2.begin(), vertices2.end(), point_cmp);

    for (size_t i = 0; i < vertices1.size(); ++i) {
      if (vertices1[i] != vertices2[i]) {
        return false;
      }
    }
    return true;
  }
  return false;
}

bool Polygon::isCongruentTo(const Shape& other) const {
  if (const Polygon* other_polygon = dynamic_cast<const Polygon*>(&other)) {
    if (this->verticesCount() != other_polygon->verticesCount()) {
      return false;
    }
    std::vector<long double> edge_lengths;
  }
  return false;
}

bool isSimilarTo(const Shape& other) const;


class Rectangle: public Polygon {
 protected:
  static std::vector<Point> calculate_vertices(const Point& p1, const Point& p3, long double k);
 public:
  Rectangle(const Point& p1, const Point& p3, long double k);

  Point center() const;
  std::pair<Line, Line> diagonals() const;
};

std::vector<Point> Rectangle::calculate_vertices(const Point& p1, const Point& p3, long double k) {
  Point c((p1.x + p3.x) / 2, (p1.y + p3.y) / 2);
  Vector diag_vec(p1, p3);

  Vector p_vec(-diag_vec.y, diag_vec.x);

  Vector side_a1 = (1 / (1 + k * k)) * (k * k * diag_vec - k * p_vec);
  Vector side_b1 = (1 / (1 + k * k)) * (diag_vec + k * p_vec);

  Vector side_a2 = (1 / (1 + k * k)) * (k * k * diag_vec + k * p_vec);
  Vector side_b2 = (1 / (1 + k * k)) * (diag_vec - k * p_vec);

  Vector short_side_cand;
  if (k < 1) {
    short_side_cand = side_a1;
  } else {
    short_side_cand = side_b1;
  }

  Point p2;
  Point p4;
  if ((diag_vec * short_side_cand) > 0) {
    p2 = p1 + side_a1;
    p4 = p1 + side_b1;
  } else {
    p2 = p1 + side_a2;
    p4 = p1 + side_b2;
  }

  return { p1, p2, p3, p4 };
}

Rectangle::Rectangle(const Point& p1, const Point& p3, long double k)
  : Polygon(calculate_vertices(p1, p3, k))
  {}

Point Rectangle::center() const {
  return Point((vertices_[0].x + vertices_[2].x) / 2, (vertices_[0].y + vertices_[2].y) / 2);
}

std::pair<Line, Line> Rectangle::diagonals() const {
  Line d1(vertices_[0], vertices_[2]);
  Line d2(vertices_[1], vertices_[3]);
  return { d1, d2 };
}


class Ellipse: public Shape {
 protected:
  Point f1_;
  Point f2_;
  long double two_a_;

 public:
  Ellipse(const Point& f1, const Point& f2, long double radius);

  std::pair<Point, Point> focuses() const;
  std::pair<Line, Line> directrices() const;
  long double eccentricity() const;
  Point center() const;
  long double focal_length() const;
  long double semimajor_axis_length() const;
  long double semiminor_axis_length() const;

  long double perimeter() const;
  long double area() const;

  void rotate(const Point& center, long double angle);
  void reflect(const Point& center);
  void reflect(const Line& axis);
  void scale(const Point& center, long double coefficient);
  bool containsPoint(const Point& point) const;

  bool operator==(const Shape& other) const;
};

Ellipse::Ellipse(const Point& f1, const Point& f2, long double two_a)
  : f1_(f1)
  , f2_(f2)
  , two_a_(two_a)
  {}

std::pair<Point, Point> Ellipse::focuses() const {
  return { f1_, f2_ };
}

std::pair<Line, Line> Ellipse::directrices() const {
  Vector major_axis_vec(f1_, f2_);
  Vector u = major_axis_vec / major_axis_vec.length();

  long double a = two_a_ / 2;
  long double e = eccentricity();
  long double d = a / e;

  Point c = center();
  Point d_point1 = c + u * d;
  Point d_point2 = c - u * d;

  Vector p = { -(u.y), u.x };

  Line d1(d_point1, d_point1 + p);
  Line d2(d_point2, d_point2 + p);

  return { d1, d2 };
}

long double Ellipse::eccentricity() const {
  long double a = two_a_ / 2;
  if (are_equal(a, 0)) {
    return 0;
  }

  long double c = focal_length();
  return c / a;
}

Point Ellipse::center() const {
  return Point((f1_.x + f2_.x) / 2, (f1_.y + f2_.y) / 2);
}

long double Ellipse::focal_length() const {
  Point c = center();
  long double x_delta = c.x - f1_.x;
  long double y_delta = c.y - f1_.y;
  return sqrtl(x_delta * x_delta + y_delta * y_delta);
}

long double Ellipse::semimajor_axis_length() const {
  long double l = two_a_ / 2.0;
  return l;
}

long double Ellipse::semiminor_axis_length() const {
  long double a = semimajor_axis_length();
  long double c = focal_length();
  long double b = sqrtl(a * a - c * c);

  return b;
}

long double Ellipse::perimeter() const {
  long double a = semimajor_axis_length();
  long double b = semiminor_axis_length();
  long double p = Pi * (3 * (a + b) - (sqrtl((3 * a + b) * (a + 3 * b))));

  return p;
}

long double Ellipse::area() const {
  long double s =  Pi * semimajor_axis_length() * semiminor_axis_length();

  return s;
}

void Ellipse::rotate(const Point& center, long double angle) {
  f1_.rotate(center, angle);
  f2_.rotate(center, angle);
}

void Ellipse::reflect(const Point& center) {
  f1_.reflect(center);
  f2_.reflect(center);
}

void Ellipse::reflect(const Line& axis) {
  f1_.reflect(axis);
  f2_.reflect(axis);
}

void Ellipse::scale(const Point& center, long double coefficient) {
  f1_.scale(center, coefficient);
  f2_.scale(center, coefficient);
  two_a_ *= coefficient;
}

bool Ellipse::containsPoint(const Point& point) const {
  Vector v1(f1_, point);
  Vector v2(f2_, point);
  long double sum = v1.length() + v2.length();

  return sum <= two_a_;
}

bool Ellipse::operator==(const Shape& other) const {
  if (const Ellipse* other_ellipse = dynamic_cast<const Ellipse*>(&other)) {
    if (((this->f1_ == other_ellipse->f1_ && this->f2_ == other_ellipse->f2_) ||
      (this->f1_ == other_ellipse->f2_ && this->f2_ == other_ellipse->f1_)) &&
      (this->two_a_ == other_ellipse->two_a_)) {
      return true;
    }
  }
  return false;
}


class Circle: public Ellipse {
 public:
  Circle(const Point& p, long double r);
  long double radius() const;

  long double perimeter() const override;
  long double area() const override;
};

Circle::Circle(const Point& p, long double r)
  : Ellipse(p, p, 2 * r)
  {}

long double Circle::radius() const {
  return two_a_ / 2;
}

long double Circle::perimeter() const {
  long double p = 2.0 * Pi * radius();

  return p;
}

long double Circle::area() const {
  long double r = radius();
  long double a = Pi * r * r;

  return a;
}


class Square: public Rectangle {
 public:
  Square(const Point& p1, const Point& p3);

  long double side() const;
  long double diagonal() const;
  Circle circumscribedCircle() const;
  Circle inscribedCircle() const;
};

Square::Square(const Point& p1, const Point& p3)
  : Rectangle(p1, p3, 1)
  {}

long double Square::side() const {
  Vector side_vec(vertices_[0], vertices_[1]);
  return side_vec.length();
}

long double Square::diagonal() const {
  Vector diag_vec(vertices_[0], vertices_[2]);
  return diag_vec.length();
}

Circle Square::circumscribedCircle() const {
  Circle tmp(center(), diagonal() / 2);
  return tmp;
}

Circle Square::inscribedCircle() const {
  Circle tmp(center(), side() / 2);
  return tmp;
}


class Triangle: public Polygon {
 public:
  Triangle(const std::vector<Point>& vertices);
  Triangle(std::initializer_list<Point> vertices);

  Circle circumscribedCircle() const;
  Circle inscribedCircle() const;
  Point centroid() const;
  Point circumcenter() const;
  Point incenter() const;
  Point orthocenter() const;
  Line EulerLine() const;
  Circle ninePointsCircle() const;
};

Triangle::Triangle(const std::vector<Point>& vertices)
  : Polygon(vertices)
  {}

Triangle::Triangle(std::initializer_list<Point> vertices)
  : Polygon(vertices)
  {}

Circle Triangle::circumscribedCircle() const {
  const Point center = circumcenter();
  const Point& A = vertices_[0];

  long double r = sqrtl((A.x - center.x) * (A.x - center.x) +
                          (A.y - center.y) * (A.y - center.y));

  return { center, r };
}

Circle Triangle::inscribedCircle() const {
  const Point center = incenter();
  const Point& A = vertices_[0];
  const Point& B = vertices_[1];
  const Point& C = vertices_[2];

  long double s = 0.5 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
  if (s < 0) {
    s = -s;
  }

  long double p = (Vector(A, B).length() +
                   Vector(B, C).length() +
                   Vector(C, A).length()) / 2.0;

  long double r = s / p;
  return { center, r };
}

Point Triangle::centroid() const {
  return { (vertices_[0].x + vertices_[1].x + vertices_[2].x) / 3.0,
           (vertices_[0].y + vertices_[1].y + vertices_[2].y) / 3.0 };
}

Point Triangle::circumcenter() const {
  const Point& A = vertices_[0];
  const Point& B = vertices_[1];
  const Point& C = vertices_[2];

  long double d = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));

  long double x = ((A.x * A.x + A.y * A.y) * (B.y - C.y) +
                   (B.x * B.x + B.y * B.y) * (C.y - A.y) +
                   (C.x * C.x + C.y * C.y) * (A.y - B.y)) / d;

  long double y = ((A.x * A.x + A.y * A.y) * (C.x - B.x) +
                   (B.x * B.x + B.y * B.y) * (A.x - C.x) +
                   (C.x * C.x + C.y * C.y) * (B.x - A.x)) / d;

  return { x, y };
}

Point Triangle::incenter() const {
  const Point& A = vertices_[0];
  const Point& B = vertices_[1];
  const Point& C = vertices_[2];

  long double a = sqrtl((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y));
  long double b = sqrtl((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y));
  long double c = sqrtl((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));

  long double x = (a * A.x + b * B.x + c * C.x) / (a + b + c);
  long double y = (a * A.y + b * B.y + c * C.y) / (a + b + c);

  return { x, y };
}

Point Triangle::orthocenter() const {
  const Point center = circumcenter();
  const Point& A = vertices_[0];
  const Point& B = vertices_[1];
  const Point& C = vertices_[2];

  Point orthocenter = center +
                      Vector(center, A) +
                      Vector(center, B) +
                      Vector(center, C);

  return orthocenter;
}

Line Triangle::EulerLine() const {
  const Point O = orthocenter();
  const Point C = circumcenter();

  return Line(O, C);
}

Circle Triangle::ninePointsCircle() const {
  Point H = orthocenter();
  Circle CC = circumscribedCircle();
  Point O = CC.center();
  long double r = CC.radius();

  Point center((H.x + O.x) / 2, (H.y + O.y) / 2);
  Circle NPC(center, r / 2);

  return NPC;
}
