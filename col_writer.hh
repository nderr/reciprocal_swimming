#ifndef COL_WRITER_HH
#define COL_WRITER_HH

#include <cstdio>
#include <cstring>

/**
 * Class for writing out a series of strings in an table of
 * columns with specified headers, widths, and padding
 *
 * Classes of 2-6 columns are explitly instantiated
 *
 * @tparam N number of columns
 */
template<int N> class col_writer {

	/** padding for each column */
	int padding[N];
	/** padding (left and right) inside dash underline */
	int lr_dpad[N][2];
	/** padding (left and right) for each column's header */
	int lr_hpad[N][2];
	/** maximum length of an entry in each column */
	int max[N];
	/** length of header of each column */
	int header_len[N];
	/** length of dash underline of each column */
	int dash_len[N];
	/** pointer to (explicit) C-string headers */
	const char* headers[N];
	/** C-string dashes */
	char* dashes;

	/**
	 * Returns the left header padding for a single column
	 *
	 * @param[in] i specified column index
	 * @return left header padding of specified column
	 */
	int hlpad(int i) {return padding[i]+lr_hpad[i][0]+(i>0?lr_hpad[i-1][1]:0);}

	/**
	 * Returns the left dash padding for a single column
	 *
	 * @param[in] i specified column index
	 * @return left dash padding of specified column
	 */
	int dlpad(int i) {return padding[i]-lr_dpad[i][0]-(i>0?lr_dpad[i-1][1]:0);}

	/**
	 * Returns the difference between a string and the
	 * maximum length of a specified column's entries
	 *
	 * @param[in] str string to examine length of
	 * @param[in] i specified column index
	 * @return different in lengths between column i's
	 *   maximum-length entry and str
	 */
	int diff(char* s,int i) {return max[i]-static_cast<int>(strlen(s));}

	public:

	// constructor / destructor
	col_writer(int padding,const char* headers[N],int* max);
	col_writer(int padding,const char* headers[N]);
	col_writer(int padding);
	~col_writer();

	// functions for updating column max entry lengths and headers,
	// plus re-calculating dash lengths, padding, etc
	void update_col_len(char *s,int c){int l=strlen(s);int &m=max[c];m=l>m?l:m;}
	void update_col_len(const char *s,int c){int l=strlen(s);int &m=max[c];m=l>m?l:m;}
	void set_header(const char *h,int c){headers[c]=h;header_len[c]=strlen(h);}
	void setup();

	// actual printing functions
	void print_headers(FILE*fh=stdout);
	void print_row(char* str[N],FILE *fh=stdout);
};

#endif
