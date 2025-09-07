#include <iostream>
#include <algorithm>
#include <cstring>

class String {
 private:
  size_t size_ = 0;
  size_t capacity_ = 0;
  char* arr_ = nullptr;

  void reallocate(size_t new_size);

 public:
  String(const char *cstyle_string);
  String(const size_t n = 0, const char c = '\0');
  String(const String& other);
  String(const char* s, size_t count);
  void swap(String& other);
  ~String();
  String& operator=(const String& other);
  friend bool operator==(const String& lhs, const String& rhs);
//  friend std::strong_ordering operator<=>(const String& lhs, const String& rhs);
  friend bool operator!=(const String& lhs, const String& rhs);
  friend bool operator<(const String& lhs, const String& rhs);
  friend bool operator>(const String& lhs, const String& rhs);
  friend bool operator<=(const String& lhs, const String& rhs);
  friend bool operator>=(const String& lhs, const String& rhs);
  char& operator[](size_t index);
  const char& operator[](size_t index) const;
  size_t length() const;
  size_t size() const;
  size_t capacity() const;
  void push_back(const char c);
  void pop_back();
  char& front();
  const char& front() const;
  char& back();
  const char& back() const;
  String& operator+=(const String& other);
  friend String operator+(const String& lhs, const String& rhs);
  String& operator+=(char c);
  friend String operator+(const String& lhs, const char rhs);
  friend String operator+(const char lhs, const String& rhs);
  friend size_t find(const String& lhs, const String& rhs);
  friend size_t rfind(const String& lhs, const String& rhs);
  size_t find(const String& substring) const;
  size_t rfind(const String& substring) const;
  String substr(size_t start, size_t count) const;
  bool empty() const;
  void clear();
  void shrink_to_fit();
  friend std::istream& operator>>(std::istream& is, String& str);
  friend std::ostream& operator<<(std::ostream& os, const String& str);
  char* data();
  const char* data() const;
};

String::String(const char* cstyle_string)
  : size_(strlen(cstyle_string))
  , capacity_(size_ + 1)
  , arr_(new char[capacity_])
{
  std::copy(cstyle_string, cstyle_string + capacity_, arr_);
}

String::String(const size_t n, const char c)
  : size_(n)
  , capacity_(n + 1)
  , arr_(new char[n + 1]) {
  std::fill(arr_, arr_ + size_, c);
  arr_[size_] = '\0';
}

String::String(const String& other)
  : size_(other.size_)
  , capacity_(other.capacity_)
  , arr_(new char[capacity_])
{
  std::copy(other.arr_, other.arr_ + capacity_, arr_);
}

String::String(const char* s, size_t count)
  : size_(count)
  , capacity_(count + 1)
  , arr_(new char[capacity_])
{
  std::copy(s, s + count, arr_);
  arr_[size_] = '\0';
}

void String::swap(String& other) {
  std::swap(size_, other.size_);
  std::swap(capacity_, other.capacity_);
  std::swap(arr_, other.arr_);
}

String::~String() {
  delete[] arr_;
}

String& String::operator=(const String& other)
{
  if (this == &other) {
    return *this;
  }

  if (other.size_ + 1 <= capacity_) {
    std::copy(other.arr_, other.arr_ + other.size_ + 1, arr_);
    size_ = other.size_;
    return *this;
  }

  String new_string(other);
  swap(new_string);
  return *this;
}

bool operator==(const String& lhs, const String& rhs) {
  if (lhs.size_ != rhs.size_) {
    return false;
  }
  return strcmp(lhs.arr_, rhs.arr_) == 0;
}

//std::strong_ordering operator<=>(const String& lhs, const String& rhs) {
//  int result = strcmp(lhs.arr_, rhs.arr_);
//  if (result < 0) {
//    return std::strong_ordering::less;
//  } else if (result > 0) {
//    return std::strong_ordering::greater;
//  }
//  return std::strong_ordering::equal;
//};

bool operator!=(const String& lhs, const String& rhs) {
  return(!(lhs == rhs));
}

bool operator<(const String& lhs, const String& rhs) {
  return strcmp(lhs.arr_, rhs.arr_) < 0;
}

bool operator>=(const String& lhs, const String& rhs) {
  return !(lhs < rhs);
}

bool operator>(const String& lhs, const String& rhs) {
  return !(lhs == rhs || lhs < rhs);
}

bool operator<=(const String& lhs, const String& rhs) {
  return rhs >= lhs;
}

char& String::operator[](size_t index) {
  return arr_[index];
}

const char& String::operator[](size_t index) const {
  return arr_[index];
}

size_t String::length() const {
  return size_;
}

size_t String::size() const {
  return size_;
}

size_t String::capacity() const {
  return capacity_ == 0 ? capacity_ : capacity_ - 1;
}

void String::reallocate(size_t new_capacity) {
  if (size_ + 1 <= new_capacity && capacity_ != new_capacity) {
    char* new_arr = new char[new_capacity];
    if (arr_ != nullptr) {
      std::copy(arr_, arr_ + size_ + 1, new_arr);
    }
    std::swap(arr_, new_arr);
    capacity_ = new_capacity;
    delete[] new_arr;
  }
}

void String::push_back(const char c) {
  if (size_ + 1 == capacity_) {
    reallocate(capacity_ == 0 ? 1 : capacity_ * 2);
  }
  arr_[size_] = c;
  arr_[size_ + 1] = '\0';
  ++size_;
}

void String::pop_back() {
  if (size_ != 0) {
    arr_[size_ - 1] = '\0';
    --size_;
  }
}

char& String::front() {
  return arr_[0];
}

const char& String::front() const {
  return arr_[0];
}

char& String::back() {
  return arr_[size_ - 1];
}

const char& String::back() const {
  return arr_[size_ - 1];
}

String& String::operator+=(const String& other) {
  size_t new_size = size_ + other.size_;
  if (new_size + 1 > capacity_) {
    reallocate(new_size + 1);
  }
  std::copy(other.arr_, other.arr_ + other.size_, arr_ + size_);
  size_ = new_size;
  arr_[size_] = '\0';
  return *this;
}

String operator+(const String& lhs, const String& rhs) {
  String new_string(lhs);
  new_string += rhs;
  return new_string;
}

String& String::operator+=(const char c) {
  push_back(c);
  return *this;
}

String operator+(const String& lhs, const char rhs) {
  String new_string(lhs.size_ + 1);
  std::copy(lhs.arr_, lhs.arr_ + lhs.size_, new_string.arr_);
  new_string.arr_[lhs.size_] = rhs;
  return new_string;
}

String operator+(const char lhs, const String& rhs) {
  String new_string(rhs.size_ + 1);
  std::copy(rhs.arr_, rhs.arr_ + rhs.size_, new_string.arr_ + 1);
  new_string.arr_[0] = lhs;
  return new_string;
}

size_t find(const String& lhs, const String& rhs) {
  size_t l_size = lhs.size_;
  size_t r_size = rhs.size_;
  for (size_t i = 0; i <= l_size - r_size; ++i) {
    if (strncmp(lhs.arr_ + i, rhs.arr_, r_size) == 0) {
      return i;
    }
  }
  return l_size;
}

size_t rfind(const String& lhs, const String& rhs) {
  if (lhs.size_ < rhs.size_) {
    return lhs.size_;
  }

  size_t l_size = lhs.size_;
  size_t r_size = rhs.size_;
  for (size_t i = l_size - r_size; ; --i) {
    if (strncmp(lhs.arr_ + i, rhs.arr_, r_size) == 0) {
      return i;
    }
    if (i == 0) {
      break;
    }
  }
  return l_size;
}

size_t String::find(const String& substring) const {
  return ::find(*this, substring);
}

size_t String::rfind(const String& substring) const {
  return ::rfind(*this, substring);
}

String String::substr(size_t start, size_t count) const {
  String new_string(arr_ + start, count);
  return new_string;
}

bool String::empty() const {
  return size_ == 0;
}

void String::clear() {
  size_ = 0;
  if (capacity_ > 0) {
    arr_[size_] = '\0';
  }
}

void String::shrink_to_fit() {
  reallocate(size_ + 1);
}

std::istream& operator>>(std::istream& is, String& str) {
  str = String();
  is >> std::ws;
  char c;
  while (is.get(c) && !std::isspace(c)) {
    str.push_back(c);
  }
  if (is && std::isspace(c)) {
    is.putback(c);
  }
  return is;
}

std::ostream& operator<<(std::ostream& os, const String& str) {
  if (str.arr_ != nullptr) {
    os << str.arr_;
  }
  return os;
}

char* String::data() {
  return arr_;
}

const char* String::data() const {
  return arr_;
}