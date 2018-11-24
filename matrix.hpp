#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <iterator>

using std::size_t;

namespace sjtu
{
	
	template <class T>
	class Matrix
	{
	private:
		size_t r, c;
		T *p;
		// your private member variables here.
		
	public:
		Matrix(size_t n = 1, size_t m = 1, T _init = T())
		{
			if (n < 0 || m < 0)
				throw std::invalid_argument("negative border!");
			r = (int) n;
			c = (int) m;
			p = new T[n * m];
			for (int i = 0; i < n * m; ++i)
				p[i] = (T) _init;
		}
		
		explicit Matrix(std::pair<size_t, size_t> sz, T _init = T())
		{
			r = (int) sz.first;
			c = (int) sz.second;
			if (r < 0 || c < 0)
				throw std::invalid_argument("negative border!");
			p = new T[r * c];
			for (int i = 0; i < r * c; ++i)
				p[i] = (T) _init;
		}
		
		Matrix(const Matrix &o)
		{
			r = o.rowLength();
			c = o.columnLength();
			p = new T[r * c];
			for (int i = 0; i < r * c; ++i)
				p[i] = o.get(i);
		}
		
		template <class U>
		Matrix(const Matrix<U> &o)
		{
			r = o.rowLength();
			c = o.columnLength();
			p = new T[r * c];
			for (int i = 0; i < r * c; ++i)
				p[i] = (T) o.get(i);
		}
		
		Matrix &operator=(const Matrix &o)
		{
			if (p == o.p) return *this;
			r = o.r;
			c = o.c;
			if (p != nullptr) delete [] p;
			p = new T[r * c];
			for (int i = 0; i < r * c; ++i)
				p[i] = o.p[i];
			return *this;
		}
		
		template <class U>
		Matrix &operator=(const Matrix<U> &o)
		{
			r = o.rowLength();
			c = o.columnLength();
			if (p != nullptr) delete [] p;
			p = new T[r * c];
			for (int i = 0; i < r * c; ++i)
				p[i] = (T) o.get(i);
			return *this;
		}
		
		Matrix(Matrix &&o) noexcept
		{
			r = o.r;
			c = o.c;
			p = o.p;
			o.p = nullptr;
		}
		
		Matrix &operator=(Matrix &&o) noexcept
		{
			r = o.r;
			c = o.c;
			if (p == o.p) return *this;
			if (p != nullptr) delete [] p;
			p = o.p;
			o.p = nullptr;
			return *this;
		}
		
		~Matrix() {if (p != nullptr) delete [] p;}
		
		Matrix(std::initializer_list<std::initializer_list<T>> il)
		{
			if (il.size() == 0)
				throw std::invalid_argument("empty Matrix!");
			r = il.size();
			auto it = il.begin();
			c = it->size();
			if (c == 0)
				throw std::invalid_argument("empty Matrix!");
			p = new T[r * c];

			for (int i = 0; i < r; ++i, ++it) {
				if (it->size() != c) {
					delete [] p;
					throw std::invalid_argument("no consistent column!");
				}

				auto ii = it->begin();
				for (int j = 0; j < c; ++j, ++ii) {
					p[i * c + j] = *ii;
				}
			}
		}
		
	public:
		size_t rowLength() const {return r;}
		
		size_t columnLength() const {return c;}
		
		T get(int t) const {
			if (t < 0 || t >= r * c)
				throw std::invalid_argument("exceed!");
			return p[t];
		}

		T &get(int t) {
			if (t < 0 || t >= r * c)
				throw std::invalid_argument("exceed!");
			return p[t];
		}

		void resize(size_t _n, size_t _m, T _init = T())
		{
			if (_n < 0 || _m < 0)
				throw std::invalid_argument("negative border!");

			if (r * c == (int) (_n * _m)) {
				r = (int) _n;
				c = (int) _m;
				return;
			}

			T *pp = new T[_n * _m];

			if (r * c > (int) (_n * _m)) {
				for (int i = 0; i < (_n * _m); ++i)
					pp[i] = p[i];
			}

			if (r * c < (int) (_n * _m)) {
				for (int i = 0; i < r * c; ++i)
					pp[i] = p[i];
				for (int i = r * c; i < (_n * _m); ++i)
					pp[i] = _init;
			}

			if (p != nullptr) delete [] p;
			p = pp;
			r = _n;
			c = _m;
		}
		
		void resize(std::pair<size_t, size_t> sz, T _init = T())
		{
			resize(sz.first, sz.second, _init);
		}
		
		std::pair<size_t, size_t> size() const
		{
			return {r, c};
		};
		
		void clear()
		{
			r = 0; c = 0;
			if (p != nullptr) delete [] p; p = nullptr;
		}
		
	public:
		const T &operator()(size_t i, size_t j) const
		{
			if (i < 0 || i >= r || j < 0 || j >= c)
				throw std::invalid_argument("exceed!!");
			return p[i * c + j];
		}
		
		T &operator()(size_t i, size_t j)
		{
			if (i < 0 || i >= r || j < 0 || j >= c)
				throw std::invalid_argument("exceed!!");
			return p[i * c + j];
		}
		
		Matrix<T> row(size_t i) const
		{
			if (i < 0 || i >= r)
				throw std::invalid_argument("exceed!!");
			Matrix<T> A(1, c);
			for (int k = 0; k < c; ++k)
				A(0, k) = p[i * c + k];
			return A;
		}
		
		Matrix<T> column(size_t i) const
		{
			if (i < 0 || i >= c)
				throw std::invalid_argument("exceed!!");
			Matrix<T> A(r, 1);
			for (int k = 0; k < r; ++k)
				A(k, 0) = p[k * c + i];
			return A;
		}
		
		
	public:
		template <class U>
		bool operator==(const Matrix<U> &o) const
		{
			if (r != o.rowLength()) return false;
			if (c != o.columnLength()) return false;
			for (int i = 0; i < r * c; ++i)
				if (p[i] != o.get(i))
					return false;
			return true;
		}
		
		template <class U>
		bool operator!=(const Matrix<U> &o) const
		{
			if (r != o.rowLength()) return true;
			if (c != o.columnLength()) return true;
			for (int i = 0; i < r * c; ++i)
				if (p[i] != o.get(i))
					return true;
			return false;
		}
		
		Matrix operator-() const
		{
			Matrix<T> A(r, c);
			for (int i = 0; i < r * c; ++i)
				A.p[i] = -p[i];
			return A;
		}
		
		template <class U>
		Matrix &operator+=(const Matrix<U> &o)
		{
			if (r != o.r || c != o.c)
				throw std::invalid_argument("not match");

			for (int i = 0; i < r * c; ++i)
				p[i] += (T) o.get(i);
			return *this;
		}
		
		template <class U>
		Matrix &operator-=(const Matrix<U> &o)
		{
			if (r != o.r || c != o.c)
				throw std::invalid_argument("not match");

			for (int i = 0; i < r * c; ++i)
				p[i] -= (T) o.get(i);
			return *this;
		}
		
		template <class U>
		Matrix &operator*=(const U &x)
		{
			for (int i = 0; i < r * c; ++i)
				p[i] *= (T) x;
			return *this;
		}
		
		Matrix tran() const
		{
			Matrix<T> A(c, r);
			for (int i = 0; i < r; ++i)
				for (int j = 0; j < c; ++j)
					A(j, i) = p[i * c + j];
			return A;
		}
		
	public: // iterator
		class iterator
		{
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type        = T;
			using pointer           = T *;
			using reference         = T &;
			using size_type         = size_t;
			using difference_type   = std::ptrdiff_t;
			
			iterator() = default;

			iterator(size_type _x, size_type _y, size_type _lx, size_type _ly, size_type rx, size_type ry, pointer _p, size_type _m, size_type _n) {
				t = _p;
				lx = _lx; ly = _ly;
				dx = rx - _lx + 1;
				dy = ry - _ly + 1;
				x = _x; y = _y;
				m = _m; n = _n;
			}
			
			iterator(const iterator &P) {
				t = P.t; lx = P.lx; ly = P.ly; dx = P.dx; dy = P.dy; x = P.x; y = P.y; m = P.m; n = P.n;
			}

			iterator &operator=(const iterator &P) {
				t = P.t; lx = P.lx; ly = P.ly; dx = P.dx; dy = P.dy; x = P.x; y = P.y; m = P.m; n = P.n;
				return (*this);
			}
			
		private:
			pointer t;
			size_type lx, ly, dx, dy, x, y, m, n;

		public:
			size_type pos() const {
				return (x - lx) * dy + y - ly;
			}
			
			difference_type operator-(const iterator &o)
			{
				return pos() - o.pos();
			}
			
			iterator &operator+=(difference_type offset)
			{
				size_type xx, yy;
				x += (xx = offset / dy);
				y += (yy = offset % dy);
				if (y >= ly + dy) y -= dy, yy -= dy, ++x, ++xx;
				if (x < lx || x >= lx + dx || y < ly || y >= ly + dy)
					throw std::invalid_argument("exceed!");
				t += (xx * n + yy);
				return *this;
			}
			
			iterator operator+(difference_type offset) const
			{
				iterator p(*this);
				p += offset;
				return p;
			}
		
			iterator &operator-=(difference_type offset)
			{
				size_type xx, yy;
				x -= (xx = offset / dy);
				y -= (yy = offset % dy);
				if (y < ly) --x, ++xx, y += dy, yy -= dy;
				if (x < lx || x >= lx + dx || y < ly || y >= ly + dy)
					throw std::invalid_argument("exceed!");
				t -= (xx * n + yy);
				return *this;
			}
			
			iterator operator-(difference_type offset) const
			{
				iterator p(*this);
				p -= offset;
				return p;
			}
			
			iterator &operator++()
			{
				++y; ++t;
				if (y == ly + dy) {
					y = ly; ++x;
					t -= dy; t += n;
				}
				// if (x < lx || x >= lx + dx || y < ly || y >= ly + dy)
					// throw std::invalid_argument("exceed!");
				return *this;
			}
			
			iterator operator++(int)
			{
				iterator p = *this;
				++(*this);
				return p;
			}
			
			iterator &operator--()
			{
				--y; --t;
				if (y < ly) {
					y = ly + dy - 1; --x;
					t += dy; t -= n;
				}
				if (x < lx || x >= lx + dx || y < ly || y >= ly + dy)
					throw std::invalid_argument("exceed!");
				return *this;
			}
			
			iterator operator--(int)
			{
				iterator p = *this;
				--(*this);
				return p;
			}
			
			reference operator*() const
			{
				return *t;
			}
			
			pointer operator->() const
			{
				return t;
			}
			
			bool operator==(const iterator &o) const
			{
				return t == o.t && x == o.x && y == o.y && m == o.m && n == o.n && lx == o.lx && ly == o.ly && dx == o.dx && dy == o.dy;
			}
			
			bool operator!=(const iterator &o) const
			{
				return t != o.t || x != o.x || y != o.y || m != o.m || n != o.n || lx != o.lx || ly != o.ly || dx != o.dx || dy != o.dy;
			}
		};
		
		iterator begin()
		{
			iterator IT(0, 0, 0, 0, r - 1, c - 1, p, r, c);
			return IT;
		}
		
		iterator end()
		{
			iterator IT(r, 0, 0, 0, r - 1, c - 1, p + r * c, r, c);
			return IT;
		}
		
		std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> ll, std::pair<size_t, size_t> rr)
		{
			iterator i1(ll.first, ll.second, ll.first, ll.second, rr.first, rr.second, p + ll.first * c + ll.second, r, c);
			iterator i2(rr.first, rr.second, ll.first, ll.second, rr.first, rr.second, p + rr.first * c + rr.second, r, c);
			return {i1, i2};
        }
	};
		
}

//
namespace sjtu
{
	template <class T, class U>
	auto operator*(const Matrix<T> &mat, const U &x)
	{
		Matrix<decltype(T() * U())> A(mat.rowLength(), mat.columnLength());
		int top = A.rowLength() * A.columnLength();
		for (int i = 0; i < top; ++i)
			A.get(i) = mat.get(i) * x;
		return A;
	}
	
	template <class T, class U>
	auto operator*(const U &x, const Matrix<T> &mat)
	{
		Matrix<decltype(T() * U())> A(mat.rowLength(), mat.columnLength());
		int top = A.rowLength() * A.columnLength();
		for (int i = 0; i < top; ++i)
			A.get(i) = mat.get(i) * x;
		return A;
	}
	
	template <class U, class V>
	auto operator*(const Matrix<U> &a, const Matrix<V> &b)
	{
		if (a.columnLength() != b.rowLength())
			throw std::invalid_argument("can't multiply!");
		int M = a.rowLength(), N = a.columnLength(), Q = b.columnLength();
		Matrix<decltype(U() * V())> A(M, Q);
		for (int i = 0; i < M; ++i)
			for (int j = 0; j < Q; ++j)
				for (int k = 0; k < N; ++k)
					A(i, j) += a(i, k) * b(k, j);
		return A;
	}
	
	template <class U, class V>
	auto operator+(const Matrix<U> &a, const Matrix<V> &b)
	{
		if (a.size() != b.size())
			throw std::invalid_argument("can't add!");
		Matrix<decltype(U() + V())> A(a.rowLength(), a.columnLength());
		int top = a.rowLength() * a.columnLength();
		for (int i = 0; i < top; ++i)
			A.get(i) = a.get(i) + b.get(i);
		return A;
	}
	
	template <class U, class V>
	auto operator-(const Matrix<U> &a, const Matrix<V> &b)
	{
		if (a.size() != b.size())
			throw std::invalid_argument("can't add!");
		Matrix<decltype(U() - V())> A(a.rowLength(), a.columnLength());
		int top = a.rowLength() * a.columnLength();
		for (int i = 0; i < top; ++i)
			A.get(i) = a.get(i) - b.get(i);
		return A;
	}
	
}

#endif //SJTU_MATRIX_HPP

