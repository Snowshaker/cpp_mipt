#include <iostream>
#include <vector>
#include <initializer_list>
#include <cmath>
#include <limits>

const long double cEPS = 1e-9;
const long double cINF = std::numeric_limits<long double>::max();

bool are_equal(long double a, long double b) {
  return std::abs(a - b) < cEPS;
}

struct Vector;

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
  long double length();
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


Vector Vector::operator-() const {
  Vector tmp(*this);
  tmp.x = -x;
  tmp.y = -y;
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

long double Vector::length() {
  return sqrtl(x * x + y * y);
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


class Shape {
 public:
//  virtual long double perimeter() const = 0;
//  virtual long double area() const = 0;
//  virtual bool operator==(const Shape& other) const = 0;
//  virtual bool isCongruentTo(const Shape& other) const = 0;
//  virtual bool isSimilarTo(const Shape& other) const = 0;
//  virtual bool containsPoint(const Point& point) const = 0;
//
//  virtual void rotate(const Point& center, long double angle) = 0;
//  virtual void reflect(const Point& center) = 0;
//  virtual void reflect(const Line& axis) = 0;
//  virtual void scale(const Point& center, long double coefficient) = 0;
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


class Circle: public Ellipse {
 public:
  Circle(const Point& p, long double r);
  long double radius() const;
};

Circle::Circle(const Point& p, long double r)
  : Ellipse(p, p, 2 * r)
  {}

long double Circle::radius() const {
  return two_a_ / 2;
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
  Circle circumscribedCircle();
  Circle inscribedCircle();
  Point centroid();
  Point orthocenter();
  Line EulerLine();
  Circle ninePointsCircle();
};