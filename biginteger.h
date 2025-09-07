#include <iostream>
#include <vector>
#include <string>

class BigInteger {
 private:
  static constexpr unsigned int cBase = 1'000'000'000;
  static constexpr int cExponent = 9;
  static constexpr unsigned int cTenPow = 1'000'000'000;

  bool is_negative;
  std::vector<unsigned int> digits;

  void normalize();
  bool has_less_abs_than(const BigInteger& other) const;
  std::pair<unsigned int, BigInteger> get_digit(const BigInteger& remainder, const BigInteger& divisor) const;
  std::pair<BigInteger, BigInteger> divide_and_rem(const BigInteger& other) const;
  unsigned int as_uint() const;
  BigInteger multiply_by_uint(unsigned int n) const;
  void multiply_by_uint_inplace(unsigned int n);
  void add_uint_inplace(unsigned int n);

 public:
  BigInteger();
  BigInteger(int n);
  BigInteger(unsigned int);
  BigInteger(long long n);
  BigInteger(unsigned long n);
  BigInteger(unsigned long long n);
  BigInteger(const BigInteger& other);
  BigInteger(const std::string& str);
  BigInteger(const char* s);
  BigInteger& operator=(const BigInteger& other);
  BigInteger operator-() const;
  BigInteger& operator+=(const BigInteger& other);
  BigInteger& operator-=(const BigInteger& other);
  friend BigInteger operator+(const BigInteger& lhs, const BigInteger& rhs);
  friend BigInteger operator-(const BigInteger& lhs, const BigInteger& rhs);
  BigInteger& operator*=(const BigInteger& other);
  friend BigInteger operator*(const BigInteger& lhs, const BigInteger& rhs);
  BigInteger& operator/=(const BigInteger& other);
  friend bool operator==(const BigInteger& lhs, const BigInteger& rhs);
  friend std::strong_ordering operator<=>(const BigInteger& lhs, const BigInteger& rhs);
  friend BigInteger operator/(const BigInteger& lhs, const BigInteger& rhs);
  BigInteger& operator%=(const BigInteger& other);
  friend BigInteger operator%(const BigInteger& lhs, const BigInteger& rhs);
  BigInteger& operator++();
  BigInteger operator++(int);
  BigInteger& operator--();
  BigInteger operator--(int);
  std::string toString() const;
  friend std::istream& operator>>(std::istream& is, BigInteger& bigint);
  friend std::ostream& operator<<(std::ostream& os, const BigInteger& bigint);
  explicit operator bool() const;
  friend BigInteger biabs(const BigInteger& b);
};

void BigInteger::normalize() {
  while (!digits.empty() && digits.back() == 0) {
    digits.pop_back();
  }

  if (digits.empty()) {
    digits.push_back(0);
    is_negative = false;
  }
}

bool BigInteger::has_less_abs_than(const BigInteger& other) const {
  if (digits.size() != other.digits.size()) {
    return digits.size() < other.digits.size();
  }
  for (int i = digits.size() - 1; i >= 0; --i) {
    if (digits[i] > other.digits[i]) {
      return false;
    } else if (digits[i] < other.digits[i]) {
      return true;
    }
  }
  return false;
}

std::pair<unsigned int, BigInteger> BigInteger::get_digit(const BigInteger& remainder, const BigInteger& divisor) const {
  unsigned int l = 0;
  unsigned int r = cBase - 1;

  unsigned int ans = l;
  BigInteger ret;

  while (l <= r) {
    unsigned int mid = l + (r - l) / 2;
    if (BigInteger cur = /* divisor.multiply_by_uint(mid) */ static_cast<unsigned long long>(mid) * divisor; cur <= remainder) {
      ans = mid;
      ret = cur;
      l = mid + 1;
    } else {
      r = mid - 1;
    }
  }
  return {ans, ret};
}

std::pair<BigInteger, BigInteger> BigInteger::divide_and_rem(const BigInteger& other) const {
  if (other == 0) {
    throw std::runtime_error("BigInteger division by zero");
  }
  if (has_less_abs_than(other)) {
    return { BigInteger(0), *this };
  }
  if (other.digits.size() == 1) {
    unsigned int divisor = other.digits[0];
    if (divisor == 0) throw std::runtime_error("BigInteger division by zero");
    BigInteger quotient;
    quotient.digits.resize(digits.size());
    unsigned long long cur_rem = 0;
    for (int i = digits.size() - 1; i >= 0; --i) {
      cur_rem = cur_rem * cBase + digits[i];
      quotient.digits[i] = static_cast<unsigned int>(cur_rem / divisor);
      cur_rem %= divisor;
    }
    quotient.is_negative = is_negative ^ other.is_negative;
    quotient.normalize();
    BigInteger remainder(cur_rem);
    remainder.is_negative = is_negative;
    remainder.normalize();
    return { quotient, remainder };
  }

  BigInteger u = biabs(*this);
  BigInteger v = biabs(other);

  unsigned int d = cBase / (v.digits.back() + 1);
  if (d > 1) {
    u.multiply_by_uint_inplace(d);
    v.multiply_by_uint_inplace(d);
  }

  int n = v.digits.size();
  int m = u.digits.size() - n;

  BigInteger q;
  q.digits.resize(m + 1, 0);

  if (u.digits.size() == (size_t)n + m) {
    u.digits.push_back(0);
  }

  for (int j = m; j >= 0; --j) {
    unsigned long long u_chunk =
      static_cast<unsigned long long>(u.digits[j + n]) * cBase + u.digits[j + n - 1];
    unsigned long long q_guess = u_chunk / v.digits.back();

    if (q_guess >= cBase) {
      q_guess = cBase - 1;
    }

    BigInteger u_segment;
    u_segment.digits.assign(u.digits.begin() + j, u.digits.begin() + j + n + 1);
    u_segment.normalize();

    BigInteger product = v.multiply_by_uint(q_guess);

    while (product > u_segment) {
      q_guess--;
      product -= v;
    }

    q.digits[j] = q_guess;

    BigInteger current_rem = u_segment - product;

    for (size_t i = 0; i < current_rem.digits.size(); ++i) {
      u.digits[j + i] = current_rem.digits[i];
    }
    for (size_t i = current_rem.digits.size(); i <= (size_t)n; ++i) {
      u.digits[j + i] = 0;
    }
  }

  u.normalize();
  BigInteger r;
  if (d > 1) {
    r.digits.resize(u.digits.size());
    unsigned long long current_rem = 0;
    for (int i = u.digits.size() - 1; i >= 0; --i) {
      current_rem = current_rem * cBase + u.digits[i];
      r.digits[i] = static_cast<unsigned int>(current_rem / d);
      current_rem %= d;
    }
  } else {
    r = u;
  }

  q.is_negative = is_negative ^ other.is_negative;
  r.is_negative = is_negative;

  q.normalize();
  r.normalize();
  return { q, r };
}

unsigned int BigInteger::as_uint() const {
  return digits[0];
}

BigInteger BigInteger::multiply_by_uint(unsigned int n) const {
  if (n == 0) return BigInteger(0);
  if (n == 1) return *this;

  BigInteger result;
  result.digits.clear();
  result.is_negative = is_negative;

  unsigned long long carry = 0;
  for (unsigned int digit : digits) {
    unsigned long long current = static_cast<unsigned long long>(digit) * n + carry;
    result.digits.push_back(static_cast<unsigned int>(current % cBase));
    carry = current / cBase;
  }

  while (carry > 0) {
    result.digits.push_back(static_cast<unsigned int>(carry % cBase));
    carry /= cBase;
  }

  result.normalize();
  return result;
}

void BigInteger::multiply_by_uint_inplace(unsigned int n) {
  if (n == 0) {
    *this = BigInteger(0);
    return;
  }
  if (n == 1) return;

  unsigned long long carry = 0;
  for (size_t i = 0; i < digits.size(); ++i) {
    unsigned long long current = static_cast<unsigned long long>(digits[i]) * n + carry;
    digits[i] = static_cast<unsigned int>(current % cBase);
    carry = current / cBase;
  }

  while (carry > 0) {
    digits.push_back(static_cast<unsigned int>(carry % cBase));
    carry /= cBase;
  }

  normalize();
}

void BigInteger::add_uint_inplace(unsigned int n) {
  if (is_negative) {
    if (digits.size() == 1 && digits[0] < n) {
      is_negative = false;
      digits[0] = n - digits[0];
    } else {
      long long borrow = n;
      for (size_t i = 0; i < digits.size() && borrow > 0; ++i) {
        long long diff = static_cast<long long>(digits[i]) - borrow;
        if (diff < 0) {
          digits[i] = static_cast<unsigned int>(diff + cBase);
          borrow = 1;
        } else {
          digits[i] = static_cast<unsigned int>(diff);
          borrow = 0;
        }
      }
      normalize();
    }
    return;
  }

  unsigned long long carry = n;
  size_t i = 0;
  while (i < digits.size() && carry > 0) {
    unsigned long long sum = static_cast<unsigned long long>(digits[i]) + carry;
    digits[i] = static_cast<unsigned int>(sum % cBase);
    carry = sum / cBase;
    ++i;
  }

  while (carry > 0) {
    digits.push_back(static_cast<unsigned int>(carry % cBase));
    carry /= cBase;
  }
}

BigInteger::BigInteger()
  : is_negative(false)
  , digits{0}
{}

BigInteger::BigInteger(int n)
  : BigInteger(static_cast<long long>(n))
{}

BigInteger::BigInteger(unsigned int n)
  : BigInteger(static_cast<unsigned long long>(n))
{}

BigInteger::BigInteger(long long n) {
  if (n == 0) {
    is_negative = false;
    digits.push_back(0);
    return;
  }

  is_negative = n < 0;
  unsigned long long abs_val = (n < 0) ? -(static_cast<unsigned long long>(n)) : n;

  do {
    digits.push_back(abs_val % cBase);
    abs_val /= cBase;
  } while (abs_val > 0);
};

BigInteger::BigInteger(unsigned long n)
  : BigInteger(static_cast<unsigned long long>(n))
{}

BigInteger::BigInteger(unsigned long long n) : is_negative(false) {
  do {
    digits.push_back(n % cBase);
    n /= cBase;
  } while (n > 0);
};

BigInteger::BigInteger(const BigInteger& other)
  : is_negative(other.is_negative)
  , digits(other.digits)
{}

BigInteger::BigInteger(const std::string& str) {
  if (str.empty() || str == "-" || str == "+") {
    *this = BigInteger(0);
    return;
  }

  bool sign = (str[0] == '-');
  size_t start_pos = (str[0] == '-' || str[0] == '+') ? 1 : 0;

  digits.clear();

  for (int end_pos = str.size(); end_pos > static_cast<int>(start_pos); end_pos -= cExponent) {
    int cur_pos = end_pos - cExponent;
    int len = cExponent;

    if (cur_pos < static_cast<int>(start_pos)) {
      len = end_pos - start_pos;
      cur_pos = start_pos;
    }
    digits.push_back(stoul(str.substr(cur_pos, len)));
  }

  is_negative = sign;
  normalize();
}

BigInteger::BigInteger(const char* s)
  : BigInteger(std::string(s))
{}

BigInteger& BigInteger::operator=(const BigInteger& other) {
  is_negative = other.is_negative;
  digits = other.digits;
  return *this;
}

BigInteger BigInteger::operator-() const {
  BigInteger tmp(*this);
  if (tmp.digits.size() == 1 && tmp.digits[0] == 0) {
    return tmp;
  }
  tmp.is_negative = !is_negative;
  return tmp;
}

BigInteger& BigInteger::operator+=(const BigInteger &other) {
  if (is_negative == other.is_negative) {
    if (digits.size() < other.digits.size()) {
      digits.resize(other.digits.size(), 0);
    }

    unsigned int carry = 0;
    for (size_t i = 0; i < digits.size(); ++i) {
      unsigned long long sum =
        static_cast<unsigned long long>(digits[i]) +
          static_cast<unsigned long long>(carry) +
          (i < other.digits.size() ? static_cast<unsigned long long>(other.digits[i]) : 0ull);

      digits[i] = static_cast<unsigned int>(sum % cBase);
      carry = static_cast<unsigned int>(sum / cBase);
    }

    if (carry > 0) {
      digits.push_back(carry);
    }

    normalize();
    return *this;
  } else {
    if (this->has_less_abs_than(other)) {
      if (is_negative) {
        *this = (BigInteger(other) += *this);
      } else {
        *this = -(-BigInteger(other) += -*this);
      }
      return *this;
    }

    unsigned int borrow = 0;
    for (size_t i = 0; i < digits.size(); ++i) {

      long long dif =
        static_cast<long long>(digits[i]) -
          (i < other.digits.size() ? static_cast<long long>(other.digits[i]) : 0ll) -
          borrow;

      if (dif < 0) {
        digits[i] = static_cast<unsigned int>(dif + static_cast<long long>(cBase));
        borrow = 1;
      } else {
        digits[i] = static_cast<unsigned int>(dif);
        borrow = 0;
      }
    }

    normalize();
    return *this;
  }
}

BigInteger& BigInteger::operator-=(const BigInteger &other) {
  *this += -other;
  return *this;
}

BigInteger operator+(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp(lhs);
  tmp += rhs;
  return tmp;
}

BigInteger operator-(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp(lhs);
  tmp -= rhs;
  return tmp;
}

BigInteger& BigInteger::operator*=(const BigInteger& other) {
  BigInteger result;
  result.digits.resize(digits.size() + other.digits.size(), 0);
  result.is_negative = is_negative ^ other.is_negative;

  for (size_t i = 0; i < digits.size(); ++i) {
    unsigned int carry = 0u;
    for (size_t j = 0; j < other.digits.size(); ++j) {
      unsigned long long current_value = result.digits[i + j] +
        static_cast<unsigned long long>(digits[i]) *
          static_cast<unsigned long long>(other.digits[j]) +
        static_cast<unsigned long long>(carry);
      result.digits[i + j] = static_cast<unsigned int>(current_value % cBase);
      carry = static_cast<unsigned int>(current_value / cBase);
    }
    if (carry > 0) {
      result.digits[i + other.digits.size()] += carry;
    }
  }

  *this = result;
  normalize();
  return *this;
}

BigInteger operator*(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp(lhs);
  tmp *= rhs;
  return tmp;
}

bool operator==(const BigInteger& lhs, const BigInteger& rhs) {
  bool lhs_is_zero = (lhs.digits.size() == 1 && lhs.digits[0] == 0);
  bool rhs_is_zero = (rhs.digits.size() == 1 && rhs.digits[0] == 0);

  if (lhs_is_zero && rhs_is_zero) return true;

  if ((lhs.digits.size() != rhs.digits.size()) ||
    lhs.is_negative != rhs.is_negative) {
    return false;
  }

  for (int i = lhs.digits.size() - 1; i >= 0; --i) {
    if (lhs.digits[i] != rhs.digits[i]) {
      return false;
    }
  }
  return true;
}

std::strong_ordering operator<=>(const BigInteger& lhs, const BigInteger& rhs) {
  if ((lhs.is_negative > rhs.is_negative) ||
    (lhs.is_negative && rhs.is_negative && rhs.has_less_abs_than(lhs)) ||
    (!lhs.is_negative && !rhs.is_negative && lhs.has_less_abs_than(rhs))) {
    return std::strong_ordering::less;
  } else if ((lhs.is_negative < rhs.is_negative) ||
    (lhs.is_negative && rhs.is_negative && lhs.has_less_abs_than(rhs)) ||
    (!lhs.is_negative && !rhs.is_negative && rhs.has_less_abs_than(lhs))) {
    return std::strong_ordering::greater;
  }
  return std::strong_ordering::equal;
};

BigInteger& BigInteger::operator/=(const BigInteger& other) {
  *this = divide_and_rem(other).first;
  return *this;
}

BigInteger operator/(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp(lhs);
  tmp /= rhs;
  return tmp;
}

BigInteger& BigInteger::operator%=(const BigInteger& other) {
  *this = divide_and_rem(other).second;
  return *this;
}

BigInteger operator%(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp(lhs);
  tmp %= rhs;
  return tmp;
}

BigInteger& BigInteger::operator++() {
  *this += BigInteger(1ull);
  return *this;
}

BigInteger BigInteger::operator++(int) {
  BigInteger tmp(*this);
  ++(*this);
  return tmp;
}

BigInteger& BigInteger::operator--() {
  *this -= BigInteger(1ull);
  return *this;
}

BigInteger BigInteger::operator--(int) {
  BigInteger tmp(*this);
  --(*this);
  return tmp;
}

std::string BigInteger::toString() const {
  if (*this == BigInteger(0)) return "0";

  std::string result = (is_negative ? "-" : "");
  result += std::to_string(digits.back());
  for (int i = digits.size() - 2; i >= 0; --i) {
    std::string chunk_str = std::to_string(digits[i]);
    result += std::string(cExponent - chunk_str.size(), '0') + chunk_str;
  }
  return result;
}

std::istream& operator>>(std::istream& is, BigInteger& bigint) {
  std::string s;
  is >> s;
  if (is) {
    bigint = BigInteger(s);
  }
  return is;
}

std::ostream& operator<<(std::ostream& os, const BigInteger& bigint) {
  std::string s = bigint.toString();
  os << s;
  return os;
}

BigInteger::operator bool() const {
  return *this != BigInteger(0ull);
}

BigInteger operator""_bi(unsigned long long n) {
  return BigInteger(n);
}

BigInteger operator""_bi(const char* s, size_t len) {
  return BigInteger(std::string(s, len));
}

BigInteger biabs(const BigInteger& b) {
  BigInteger tmp(b);
  tmp.is_negative = false;
  return tmp;
}

class Rational {
 private:
  BigInteger num;
  BigInteger den;

  void normalize();

 public:
  friend BigInteger gcd(BigInteger a, BigInteger b);
  Rational();
  Rational(const BigInteger& n);
  Rational(int n);
  Rational operator-() const;
  Rational& operator+=(const Rational& other);
  friend Rational operator+(const Rational& lhs, const Rational& rhs);
  Rational& operator-=(const Rational& other);
  friend Rational operator-(const Rational& lhs, const Rational& rhs);
  Rational& operator*=(const Rational& other);
  friend Rational operator*(const Rational& lhs, const Rational& rhs);
  Rational& operator/=(const Rational& other);
  friend Rational operator/(const Rational& lhs, const Rational& rhs);
  friend bool operator==(const Rational& lhs, const Rational& rhs);
  friend std::strong_ordering operator<=>(const Rational& lhs, const Rational& rhs);
  std::string toString() const;
  std::string asDecimal(size_t precision = 0) const;
  explicit operator double() const;
};

BigInteger gcd(BigInteger a, BigInteger b) {
  a = biabs(a);
  b = biabs(b);

  while (b) {
    a %= b;
    std::swap(a, b);
  }
  return a;
}

void Rational::normalize() {
  BigInteger gcd_val = gcd(num, den);
  num /= gcd_val;
  den /= gcd_val;

  if (den < 0) {
    num = -num;
    den = -den;
  }
}

Rational::Rational()
  : num(0)
  , den(1)
{}

Rational::Rational(const BigInteger& n)
  : num(n)
  , den(1)
{}

Rational::Rational(int n)
  : Rational(BigInteger(n))
{}

Rational Rational::operator-() const {
  Rational tmp(*this);
  tmp.num = -tmp.num;
  return tmp;
}

Rational& Rational::operator+=(const Rational& other) {
  num *= other.den;
  num += other.num * den;
  den *= other.den;
  normalize();
  return *this;
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
  Rational tmp(lhs);
  tmp += rhs;
  return tmp;
}

Rational& Rational::operator-=(const Rational& other) {
  num *= other.den;
  num -= other.num * den;
  den *= other.den;
  normalize();
  return *this;
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
  Rational tmp(lhs);
  tmp -= rhs;
  return tmp;
}

Rational& Rational::operator*=(const Rational& other) {
  num *= other.num;
  den *= other.den;
  normalize();
  return *this;
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
  Rational tmp(lhs);
  tmp *= rhs;
  return tmp;
}

Rational& Rational::operator/=(const Rational& other) {
  Rational tmp(other);
  std::swap(tmp.num, tmp.den);
  *this *= tmp;
  return *this;
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
  Rational tmp(lhs);
  tmp /= rhs;
  return tmp;
}

bool operator==(const Rational& lhs, const Rational& rhs) {
  return (lhs <=> rhs) == std::strong_ordering::equal;
}

std::strong_ordering operator<=>(const Rational& lhs, const Rational& rhs) {
  BigInteger lnum = lhs.num * rhs.den;
  BigInteger rnum = rhs.num * lhs.den;
  return lnum <=> rnum;
}

std::string Rational::toString() const {
  std::string result = num.toString();
  if (den != 1) {
    result += "/" + den.toString();
  }

  return result;
}

std::string Rational::asDecimal(size_t precision) const {
  bool sign(num < 0);
  BigInteger n(biabs(num));

  BigInteger qnt = n / den;
  std::string result = (qnt).toString();

  BigInteger rem = n % den;

  if (sign && (qnt || rem)) {
    result = '-' + result;
  }

  if (rem) {
    result += '.';
    for (size_t i = 0; i < precision; ++i) {
      rem *= 10;
      result += (rem / den).toString();
      rem %= den;
    }
  }

  return result;
}

Rational::operator double() const {
  return std::stod(asDecimal(20));
}