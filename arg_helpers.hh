///////////////////////////////////////////////////////////////////////////////
// arg_helpers.hh
//
// This file contains a collection of helper structs for the arg_manager class,
// used to read in and parse C++ command-line arguments.
//
// Nick Derr, Nov 12, 2022
///////////////////////////////////////////////////////////////////////////////

#ifndef ARG_HELP_HH
#define ARG_HELP_HH

#include <typeinfo>
#include <cstdio>
#include <vector>
#include <cstdlib>

/** type identification enum for comparisons */
enum type_ID {INT,DOUBLE,STRING,BOOL};

/**
 * Static class with information about the data type
 * serving as a template parameter
 *
 * @tparam T the type in question
 */
template <typename T> struct type {

	/** type identifier */
	static const type_ID ID;
	/** type name in string form */
	static const char* name;

	// I/O
	static char fmt();
	static void fmt(char *str,int n);
	static void sprint(char *str,int n,T val);
	static void fprint(FILE *fh,int n,T val);
	static bool parse(char *str,T &val);
};

/**
 * Class representing an option which can be specified,
 * with or without an argument, at the command line.
 * A type-specific subclass for any requested types
 * must be implemented
 */
struct option_base {

	/** whether the option was provided at the command line */
	bool present;
	/** 1-character flag */
	char flag;
	/** identifying name (should be unique) */
	const char* name;
	/** longer description of option's effects */
	const char* desc;
	
	/** constructs option with provided flag, name and description */
	option_base(char fl,const char* name_,const char* descr_):
		present(false),flag(fl),name(name_),desc(descr_) {}

	/** returns whether function is of type T */
	template <typename T> bool is() {return get_type_ID()==type<T>::ID;}
	/** virtual function returning unique type index */
	virtual type_ID get_type_ID()=0;
	/** virtual function returning unique type name */
	virtual const char* type_name()=0;
	/** prints default value to a C-string */
	virtual void sprint_def(char *fh,int n)=0;
	/** prints default value to a filestream */
	virtual void fprint_def(FILE *fh,int n)=0;
};

/**
 * Type-specific subclasses of the option_base abstract class
 */
template <typename T> struct option : public option_base {

	/** value to use in program execution */
	T val;
	/** default value */
	T def_val;

	/** constructor */
	option(char fl,T def,const char* name_,const char* descr_):
		option_base(fl,name_,descr_),val(def),def_val(def) {}

	/** prints default value to a C-string */
	virtual void sprint_def(char *fh,int n) {type<T>::sprint(fh,n,def_val);}
	/** prints default value to a filestream */
	virtual void fprint_def(FILE *fh,int n) {type<T>::fprint(fh,n,def_val);}
	/** returns the type index */
	type_ID get_type_ID() {return type<T>::ID;}
	/** returns the type name */
	const char* type_name() {return type<T>::name;}
	/** parses the provided string as a data member of this type */
	bool parse(char *str) {return type<T>::parse(str,val);}
};

#endif

