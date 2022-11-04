#ifndef ARG_MANAGER_HH
#define ARG_MANAGER_HH

#include "arg_helpers.hh"
#include "col_writer.hh"
#include <vector>
#include <typeinfo>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <map>

/**
 * A class for managing command-line args and
 * input parameters of C++ programs
 */
class arg_manager {

	// output constants

	/**  max length of a string */
	static const unsigned MAX_STR_LEN=1024;
	/** indent amount for description printing */
	static const unsigned INDENT=6;
	/** columns used for description printing */
	static const unsigned COLS=56;

	// static helper functions for string ID, printing, quitting
	static int to_int(char c) {return c-'0';}
	static int to_flag(int i) {return '0'+i;}
	static bool is_digit(char c) {
		return '0' <= c && c <= '9';
	}
	static bool is_letter(char c) {
		return ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z');
	}
	void check_letter(char c) {
		if (is_letter(c)) return;
		quit("error: flag %c is not a letter\n",c);
	}
	static bool is_flag(char fl,char* str);
	void quit(const char* format,...);
	static void print_indent_wrap(const char *str,
			unsigned indent,unsigned cols);

	// arg variables

	/** number of provided args */
	int argc;
	/** number of required args */
	int n_req;
	/** number of required args read in */
	int n_read;
	/** provided arg values */
	char** argv;
	/** flag marking whether an arg has been processed */
	bool* used;
	/** whether help has been loaded */
	bool help_ready;
	/** whether the required args are good */
	bool args_good;

	// print variable
	col_writer<4> cw;

	// arg storage

	/** map (by unique flag) of all options provided */
	std::map<char,std::pair<type_ID,int> > opts;
	/** storage vector for double options */
	std::vector<option<double> > doubles;
	/** storage vector for int options */
	std::vector<option<int> > ints;
	/** storage vector for string options */
	std::vector<option<char*> > strings;
	/** storage vector for bool options (switches) */
	std::vector<option<bool> > bools;

	// vector of allocated cstrings
	std::vector<char*> allocs;

	// functions for getting pointers to specific options
	option_base* opt(type_ID typ,int ind);
	option_base* opt(char flag);

	// pointer to option storage vector
	template<class C> std::vector<option<C> >* opt_storage();

	// adding/setting/getting options
	template<class C> void option_set(option<C> &a);
	template<class C> void add_option(char flag,
			C def,const char* name,const char* descr);
	template<class C> void add_req(const char* name,const char* desc) {

		// check arg number
		if (n_read==n_req) {
			quit("error: attempt to add required arg #%d, but %d required args are"
					" specified\n",n_read+1,n_req);
		}
		add_option<C>(to_flag(n_read++),0,name,desc);
	}
	template<class C> C option_value(char fl);
	template<class C> C req_value(int i) {
		check_ind(i);
		return option_value<C>(to_flag(i));
	}

	// printing
	void print_help(bool verbose);
	bool help() {return help_ready&&(present('h')||present('H'));}
	void process_help();
	void process_reqs();
	void print_usage();

	// string helpers
	int flag_pos(char fl);
	void check_str(const char *str,char flag,bool def=false);
	void check_def(const char *str,char flag) {check_str(str,flag,true);}

	// helper for # of req args
	void check_ind(int i) {
		if (0<=i && i<n_req) return;
		quit("error: index %i is out of range of %d required args\n",
				i,n_req);
	}

	public:

	// constructor/destructor
	arg_manager(int c,char** v,int n_req=0);
	~arg_manager();

	// the function for processing all provided command-line args
	void process_args();

	// functions for adding options (just wrappers around templated add)
	typedef const char* sl;
	void add_double(char f,double v,sl n,sl d) {
		check_letter(f);
		add_option<double>(f,v,n,d);
	}
	void add_int(char f,int v,sl n,sl d) {
		check_letter(f);
		add_option<int>(f,v,n,d);
	}
	void add_switch(char f,sl n,sl d) {
		check_letter(f);
		add_option<bool>(f,false,n,d);
	}
	void add_string(char flag,char* def,sl name,sl desc);
	void add_string(char flag,sl def,sl name,sl desc);

	// functions for adding required args
	void add_double(sl name,sl desc) {add_req<double>(name,desc);}
	void add_int(sl name,sl desc) {add_req<int>(name,desc);}
	void add_string(sl name,sl desc);

	// functions for getting options / whether they've been provided
	double get_double(char flag) {return option_value<double>(flag);}
	int get_int(char flag) {return option_value<int>(flag);}
	char* get_string(char flag) {return option_value<char*>(flag);}
	bool get_switch(char flag) {return option_value<bool>(flag);}
	bool present(char fl) {return opt(fl)->present;}

	// functions for getting required args
	int get_int(int i) {return req_value<int>(i);}
	double get_double(int i) {return req_value<double>(i);}
	char* get_string(int i) {return req_value<char*>(i);}
};

#endif
