#include <immintrin.h>
#include <assert.h> 
#include <iostream>
#include <iomanip> 

const int INF = 0x3F3F3F3F;

template<class E> 
class Matrix {

	public:
    Matrix() { 
    	build(0, 0, 0); 
    }
    
    Matrix(int dim1, int dim2) { 
    	build(dim1, dim2, 0); 
    }
    
    Matrix(int dim1, int dim2, bool align) { 
    	build(dim1, dim2, align); 
    }
    
    Matrix(int dim1, int dim2, bool align, const E &val) { 
    	build(dim1, dim2, align); fill(val); 
    }

    Matrix(const Matrix<E>& m2){
    	build(m2.m_dim1, m2.m_dim2, m2.m_align);
    	if (m_pStart != 0) {
				E *pSrc  = m2.m_pStart;
				for (int i = 0; i < m_size; ++i){
					m_pStart[i] = pSrc[i];	
				}
			}
    }

    // Returns dimension 1 length.
    int dim1() const { return m_dim1; }

    // Returns dimension 2 length.
    int dim2() const { return m_dim2; }
    
    // Returns allocated (dim1*dim2) length.
    int size() const { return m_size; }
    
    void construct(int dim1, int dim2) {
    	construct(dim1, dim2, false);
    }

    void construct(int dim1, int dim2, bool align, const E &val) {
    	if (m_pStart!=NULL)
    		destroy();	
    	build(dim1, dim2, align);
    	fill(val);
    }
    
    void construct(int dim1, int dim2, bool align) {
    	if (m_pStart!=NULL)
    		destroy();	
    	build(dim1, dim2, align);
    }
    
    // Fill the matrix with val.
		void fill(const E &val) {
			for (int i = 0; i < m_size; ++i)
				m_pStart[i] = val;	
		}
	
		// Set matrix position at (i, j) with val.
		void set(int i, int j, const E &val) {
			assert(0 <= i && i <= m_dim1 && 0 <= j && j <= m_dim2);
			m_pStart[i*m_dim2+j] = val;
		}
	
		// Set matrix position at index with val.
		void set(int index, const E &val) {
			assert(0 <= index && index < m_size);
			m_pStart[index] = val;
		}
	
		int getCollumnSize() const {
			return m_dim2;
		}
	
		int getSize() const {
			return m_size;
		}
	
		E* getAddress(int index) const {
			return &m_pStart[index];
		}
	
		// Return a reference to the element with index (i, j).
		E &operator()(int i, int j) const {
			assert(0 <= i && i <= m_dim1 && 0 <= j && j <= m_dim2);
			return m_pStart[i*m_dim2+j];
		}
	
		// Return a reference to the element at index.
		E &operator()(int index) const {
			assert(0 <= index && index < m_size);
			return m_pStart[index];
		}
	
		Matrix<E> &operator=(const Matrix<E> &m2) {
			copy(m2);
			return *this;
		}

		bool operator==(const Matrix<E> &m2) {
			if (m_pStart == 0) {
				if (m2.m_pStart == 0) return true;
				else return false;
			}
			if (m2.m_pStart == 0) return false;
			if (m_size != m2.m_size) return false;

			for (int i = 0; i < m_size; i++){
				if (m_pStart[i] != m2.m_pStart[i])
					return false;
			}
			return true;
		}
	
		void copy(const Matrix<E> &array2){
			if (!(array2.m_dim1 == m_dim1 && array2.m_dim2 == m_dim2 && array2.m_align == m_align)){
				destroy();
				build(array2.m_dim1, array2.m_dim2, array2.m_align);
			}
			if (m_pStart != 0) {
				E *pSrc  = array2.m_pStart;
				for (int i = 0; i < m_size; ++i){
					m_pStart[i] = pSrc[i];	
				}
			}
		}
		
	 friend std::ostream& operator<<(std::ostream &os, const Matrix<E> &m){
      os << "{";
      for(int i = 0; i < m.m_dim1; i++){
          if(i) os << "},\n ";
          os << "{";
          for(int j=0; j< m.m_dim2; j++){
              if(j) os << ", ";
              if ((int)m(i,j) == INF) os << "∞";
              else os << m(i, j);
          }
      }
      return os << "}}";
    }

		~Matrix() { 
			destroy();
		}
    
	private:

		E   *m_pStart; 	// Start of the array (address of A[0,0]).
		int  m_dim1;
		int  m_dim2;
		int  m_size;	// m_dim1 * m_dim 2
		bool  m_align;

		// Allocates a matrix as an array of length m_size.
		void build(int dim1, int dim2, bool align) {
		    m_dim1 = dim1;
		    m_dim2 = dim2;
		    m_size = dim1*dim2;
		    m_align = align;
		    
			if (m_size == 0) {
				m_pStart = NULL;
			} else {	
				if (!align)
					m_pStart = (E*) malloc (m_size * sizeof(E));
				else
					m_pStart = (E*) _mm_malloc (m_size * sizeof(E), 32);
			}
		}
		
		void destroy() {
			if (m_pStart!=NULL)
				free(m_pStart);
		}
			   	
};

